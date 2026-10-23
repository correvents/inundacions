#!/usr/bin/env python3
"""Passos inferiors: obre al relleu els passos per sota de terraplens que el MDT no veu.

El MDT del LiDAR és el terreny: un pas inferior, un túnel curt o un culvert sota una carretera o una via elevada hi queden
tapats pel terraplè, que fa de presa i embassa l'aigua aigües amunt. Amb OpenStreetMap:
  - Trams amb tunnel=yes/culvert (o layer negatiu) de menys de --llarg-max m: el pas és el mateix tram.
  - Trams amb bridge=*: cada via, riera o camí que hi passa per sota (sense pont i que el creua) és un pas, en un tram de
    --mig m a cada banda del creuament.
A cada pas, el fons del pas es fa anar en línia recta entre el terreny més baix just fora de cada boca (--boca m més enllà)
i les cel·les del MDT que travessa el pas que quedin per sobre d'aquesta línia es rebaixen fins a la línia. Si el pas ja és
obert al MDT (cap cel·la per sobre), no es toca res.

Sortides: el mateix --mdt modificat, data/<conca>/passos/osm.json (la consulta, per reproduir-ho) i
data/<conca>/passos/passos.geojson (cada pas, quant s'ha rebaixat i on), per veure'ls al visor.

Ús:  python3 scripts/passos.py --config config/granollers.json --mdt /tmp/relleu/mdt.tif"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import requests
import rasterio
from pyproj import Transformer
from shapely.geometry import LineString, Point, mapping
from shapely.ops import substring

OVERPASS = ["https://overpass-api.de/api/interpreter", "https://overpass.kumi.systems/api/interpreter"]
CAIXES = {"ridaura": [2.86, 41.75, 3.09, 41.89]}
TIPUS = ("highway", "railway", "waterway")


def consulta(caixa):
    w, s, e, n = caixa
    q = f"""[out:json][timeout:180];
(way["bridge"]["bridge"!="no"]({s},{w},{n},{e});
 way["tunnel"]["tunnel"!="no"]({s},{w},{n},{e});
 way["layer"~"^-"]({s},{w},{n},{e});
 way["waterway"]({s},{w},{n},{e});
 way["highway"]({s},{w},{n},{e});
 way["railway"]({s},{w},{n},{e}););
out tags geom;"""
    for url in OVERPASS:
        try:
            r = requests.post(url, data={"data": q}, timeout=300)
            r.raise_for_status()
            return r.json()
        except Exception as ex:
            print(f"{url}: {ex}")
    raise RuntimeError("Overpass no respon")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--mdt", required=True)
    ap.add_argument("--llarg-max", type=float, default=120, help="llargada màxima d'un túnel per tractar-lo com a pas inferior (m)")
    ap.add_argument("--mig", type=float, default=40, help="mitja amplada màxima d'un terraplè sota un pont (m)")
    ap.add_argument("--boca", type=float, default=20, help="distància fora de cada boca on es busca el terreny (m)")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    conca = cfg["conca"]
    d = Path("data") / conca / "passos"; d.mkdir(parents=True, exist_ok=True)
    fosm = d / "osm.json"
    if fosm.exists():
        osm = json.loads(fosm.read_text())
    else:
        osm = consulta(cfg.get("caixa") or CAIXES[conca])
        fosm.write_text(json.dumps(osm, separators=(",", ":")))

    with rasterio.open(a.mdt) as f:
        z = f.read(1); tr = f.transform; crs = f.crs; perfil = f.profile
    a_utm = Transformer.from_crs(4326, crs, always_xy=True).transform
    a_ll = Transformer.from_crs(crs, 4326, always_xy=True).transform
    vies = []
    for el in osm.get("elements", []):
        t = el.get("tags", {})
        if el.get("type") != "way" or "geometry" not in el or not any(k in t for k in TIPUS):
            continue
        xy = [a_utm(p["lon"], p["lat"]) for p in el["geometry"]]
        if len(xy) < 2:
            continue
        vies.append({"id": el["id"], "t": t, "g": LineString(xy)})
    ponts = [v for v in vies if v["t"].get("bridge", "no") != "no"]
    def es_tunel(v):
        return v["t"].get("tunnel", "no") != "no" or v["t"].get("layer", "0").startswith("-")
    print(f"OSM: {len(vies)} vies, {len(ponts)} ponts, {sum(map(es_tunel, vies))} túnels/culverts")

    # passos candidats: (geometria del pas, descripció)
    passos = []
    for v in vies:
        if es_tunel(v) and v["g"].length <= a.llarg_max:
            passos.append((v["g"], v, "túnel/culvert"))
    from shapely.strtree import STRtree
    baixes = [v for v in vies if v["t"].get("bridge", "no") == "no" and not es_tunel(v)]
    arbre = STRtree([v["g"] for v in baixes])
    for p in ponts:
        for k in arbre.query(p["g"], predicate="intersects"):
            v = baixes[k]
            if v["t"].get("layer", "0") == p["t"].get("layer", "1") or v["g"].touches(p["g"]):
                continue                                                   # es toquen pels extrems: continuació, no creuament
            x = v["g"].intersection(p["g"])
            if x.is_empty or x.geom_type not in ("Point", "MultiPoint"):
                continue
            for pt in getattr(x, "geoms", [x]):
                s = v["g"].project(pt)
                seg = substring(v["g"], max(s - a.mig, 0), min(s + a.mig, v["g"].length))
                if seg.length > 1:
                    passos.append((seg, v, f"sota pont ({p['t'].get('highway') or p['t'].get('railway') or 'pont'})"))

    nf, nc = z.shape
    def mostra(pt):
        c, r = ~tr * (pt.x, pt.y)
        return (int(r), int(c)) if 0 <= r < nf and 0 <= c < nc else None
    def terra_boca(g, extrem):
        """terreny més baix entre la boca i --boca m més enllà, en la direcció del pas"""
        p0, p1 = (g.coords[0], g.coords[1]) if extrem == 0 else (g.coords[-1], g.coords[-2])
        dx, dy = p0[0] - p1[0], p0[1] - p1[1]; L = np.hypot(dx, dy) or 1
        vals = []
        for s in np.arange(0, a.boca + 0.1, tr.a / 2):
            rc = mostra(Point(p0[0] + dx / L * s, p0[1] + dy / L * s))
            if rc and np.isfinite(z[rc]):
                vals.append(z[rc])
        return min(vals) if vals else None
    fets, fts = 0, []
    for g, v, tipus in passos:
        z0, z1 = terra_boca(g, 0), terra_boca(g, 1)
        if z0 is None or z1 is None:
            continue
        rebaix, cel = 0.0, set()
        for s in np.arange(0, g.length + 0.1, tr.a / 3):
            rc = mostra(g.interpolate(s))
            if not rc or rc in cel or not np.isfinite(z[rc]):
                continue
            cel.add(rc)
            fons = z0 + (z1 - z0) * s / g.length
            if z[rc] > fons:
                rebaix = max(rebaix, float(z[rc] - fons)); z[rc] = fons
        if rebaix > 0.3:
            fets += 1
            q = LineString([a_ll(*c) for c in g.coords])
            fts.append({"type": "Feature", "geometry": mapping(q), "properties": {
                "osm": v["id"], "tipus": tipus, "via": v["t"].get("name") or v["t"].get("highway") or v["t"].get("waterway") or v["t"].get("railway"),
                "rebaix_m": round(rebaix, 1), "fons_m": [round(float(z0), 1), round(float(z1), 1)]}})
    fts.sort(key=lambda f: -f["properties"]["rebaix_m"])
    (d / "passos.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": fts}, ensure_ascii=False))
    with rasterio.open(a.mdt, "w", **perfil) as f:
        f.write(z, 1)
    print(f"{len(passos)} passos candidats; {fets} obren un terraplè de més de 30 cm. Els més grans:")
    for f in fts[:15]:
        p = f["properties"]; c = f["geometry"]["coordinates"][len(f["geometry"]["coordinates"]) // 2]
        print(f"  {p['rebaix_m']:5.1f} m  {p['tipus']:28s} {p['via']}  ({c[1]:.5f}, {c[0]:.5f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
