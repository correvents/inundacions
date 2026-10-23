#!/usr/bin/env python3
"""Passos inferiors: obre al relleu els passos per sota de terraplens que el MDT no veu.

El MDT del LiDAR és el terreny: un pas inferior, un túnel curt o un culvert sota una carretera o una via elevada hi queden
tapats pel terraplè, que fa de presa i embassa l'aigua aigües amunt. Amb OpenStreetMap:
  - Trams amb tunnel=yes/culvert (o layer negatiu) de menys de --llarg-max m: el pas és el mateix tram.
  - Trams amb bridge=*: cada via, riera o camí que hi passa per sota (sense pont i que el creua) és un pas, en un tram de
    --mig m a cada banda del creuament.
  - Rieres i torrents d'OSM que creuen una carretera o una via de tren sense pont ni culvert dibuixat: s'hi suposa l'obra
    de drenatge transversal que sempre s'hi fa, en el mateix tram de --mig m. Es marquen com a «suposat».
A cada pas, el fons del pas es fa anar en línia recta entre el terreny més baix just fora de cada boca (--boca m més enllà)
i les cel·les del MDT que travessa el pas que quedin per sobre d'aquesta línia es rebaixen fins a la línia. Si el pas ja és
obert al MDT (cap cel·la per sobre), no es toca res.

Capacitat de les obres de drenatge (rieres sota terraplè, dibuixades a OSM o suposades): no se'n sap la mida, però es
dissenyen amb la norma de carreteres 5.2-IC perquè hi passi la crescuda de 100 anys de la seva conca. Es calcula així,
sense cap número triat a mà:
  - conca de l'obra: direccions de drenatge sobre el relleu ja obert (pysheds), àrea, longitud i desnivell;
  - pluja de 100 anys en 24 h: Gumbel regional amb els màxims anuals de pluja diària de les estacions XEMA de la conca
    (× 1,13 per passar de dia fix a 24 h mòbils);
  - intensitat per al temps de concentració de Témez amb la corba de la 5.2-IC (I1/Id = 11 al litoral català),
    coeficient d'uniformitat Kt i reducció per àrea Ka de la mateixa norma;
  - escolament amb el llindar P0 = absorció inicial mitjana del model a la conca (mediana dels episodis).
  Q = C · I · A · Kt / 3,6. La simulació 2D no deixa passar per l'obra més d'aquest cabal (--passos-cap de simula_2d.py).
  Ponts i passos de carretera no es limiten: són oberts.

Sortides: el mateix --mdt modificat, data/<conca>/passos/osm.json (la consulta, per reproduir-ho),
data/<conca>/passos/passos.geojson (cada pas, quant s'ha rebaixat, on i, si és una obra de drenatge, el cabal de disseny)
i data/<conca>/passos/capacitat.json (cel·les de cada obra i cabal màxim) per a la simulació 2D.

Ús:  python3 scripts/passos.py --config config/granollers.json --mdt /tmp/relleu/mdt.tif"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import requests
import rasterio
from pyproj import Transformer
from shapely.geometry import LineString, Point, mapping
from shapely.ops import substring

OVERPASS = ["https://overpass.private.coffee/api/interpreter", "https://overpass-api.de/api/interpreter",
            "https://overpass.kumi.systems/api/interpreter", "https://maps.mail.ru/osm/tools/overpass/api/interpreter"]
UA = {"User-Agent": "inundacions/0.1 (github.com/correvents/inundacions)"}
CAIXES = {"ridaura": [2.86, 41.75, 3.09, 41.89]}
TIPUS = ("highway", "railway", "waterway")
RIERES = ("river", "stream", "ditch", "drain", "canal")
VIALS = ("motorway", "trunk", "primary", "secondary", "tertiary", "unclassified", "residential",
         "motorway_link", "trunk_link", "primary_link", "secondary_link", "tertiary_link")


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
    for url in OVERPASS * 2:
        try:
            r = requests.post(url, data={"data": q}, headers=UA, timeout=300)
            r.raise_for_status()
            return r.json()
        except Exception as ex:
            print(f"{url}: {ex}")
    raise RuntimeError("Overpass no respon")


def p24_t100(conca):
    """Pluja de 100 anys en 24 h: Gumbel regional (índex de crescuda) amb els màxims anuals de pluja diària de la XEMA."""
    norm, mitj = [], []
    for f in sorted((Path("data") / conca / "xema").glob("*.csv")):
        if f.name.startswith("episodis"):
            continue
        s = pd.read_csv(f, parse_dates=["temps"]).set_index("temps").mm
        dia = s.resample("1D").agg(["sum", "count"])
        dia = dia[dia["count"] >= 44]                                  # dies amb ≥ 22 h de dades
        anual = dia.groupby(dia.index.year)["sum"].agg(["max", "count"])
        anual = anual[anual["count"] >= 330]["max"]                    # anys gairebé complets
        if len(anual) >= 5:
            mitj.append(anual.mean()); norm += list(anual / anual.mean())
    if len(norm) < 20:
        return None, len(norm)
    x = np.array(norm); beta = x.std(ddof=1) * np.sqrt(6) / np.pi; mu = x.mean() - 0.5772 * beta
    return 1.13 * float(np.mean(mitj)) * (mu - beta * np.log(-np.log(1 - 1 / 100))), len(norm)


def q_disseny(area_km2, l_km, j, p24, p0):
    """Cabal de disseny de la 5.2-IC (mètode racional) per a una conca petita."""
    tc = max(0.3 * (l_km / max(j, 1e-3) ** 0.25) ** 0.76, 0.25)       # Témez, h
    ka = 1 - np.log10(area_km2) / 15 if area_km2 > 1 else 1.0
    pd_ = p24 * ka; id_ = pd_ / 24
    i = id_ * 11 ** ((28 ** 0.1 - tc ** 0.1) / (28 ** 0.1 - 1))
    kt = 1 + tc ** 1.25 / (tc ** 1.25 + 14)
    r = pd_ / p0
    c = (r - 1) * (r + 23) / (r + 11) ** 2 if r > 1 else 0.0
    return c * i * area_km2 * kt / 3.6, tc, c


DIRS = {64: (-1, 0), 128: (-1, 1), 1: (0, 1), 2: (1, 1), 4: (1, 0), 8: (1, -1), 16: (0, -1), 32: (-1, -1)}   # pysheds: N, NE, E, SE, S, SO, O, NO


def recorregut(fd, mask, r0, c0):
    """Recorregut d'aigua més llarg fins a la sortida (r0, c0) dins la conca, en cel·les (D8; diagonals √2)."""
    rr, cc = np.nonzero(mask)
    dist = np.full(len(rr), np.nan); pos = {(r, c): k for k, (r, c) in enumerate(zip(rr, cc))}
    avall = np.full(len(rr), -1); pas = np.ones(len(rr))
    for k, (r, c) in enumerate(zip(rr, cc)):
        d = DIRS.get(int(fd[r, c]))
        if d and (r, c) != (r0, c0):
            avall[k] = pos.get((r + d[0], c + d[1]), -1); pas[k] = 1.4142 if d[0] and d[1] else 1.0
    dist[pos[(r0, c0)]] = 0.0
    for _ in range(len(rr)):
        nou = np.where(np.isnan(dist) & (avall >= 0), dist[np.maximum(avall, 0)] + pas, dist)
        if np.array_equal(np.isnan(nou), np.isnan(dist)):
            break
        dist = nou
    return float(np.nanmax(dist))


def capacitats(z, tr, crs, obres, conca):
    """Conca i cabal de disseny de cada obra de drenatge. obres: llista de (feature, [(fila, col) rebaixades], boca alta (x, y))."""
    p24, n = p24_t100(conca)
    idx = Path("data") / conca / "episodis" / "index.json"
    ias = [e["ia_mm"] for e in json.loads(idx.read_text())["episodis"] if e.get("ia_mm")] if idx.exists() else []
    if not p24 or not ias:
        print(f"Sense capacitat de les obres: {n} anys-estació de XEMA (cal ≥ 20) i {len(ias)} episodis amb absorció")
        return None
    p0 = float(np.median(ias))
    print(f"Pluja de 100 anys en 24 h: {p24:.0f} mm ({n} anys-estació XEMA); P0 = {p0:.0f} mm")
    import tempfile
    np.in1d = np.isin                                                  # pysheds encara fa servir np.in1d (numpy 2)
    from pysheds.grid import Grid
    tmp = Path(tempfile.mkdtemp()) / "z.tif"
    with rasterio.open(tmp, "w", driver="GTiff", height=z.shape[0], width=z.shape[1], count=1, dtype="float32", crs=crs, transform=tr, nodata=-9999) as f:
        f.write(np.where(np.isfinite(z), z, -9999).astype("float32"), 1)
    g = Grid.from_raster(str(tmp)); dem = g.read_raster(str(tmp))
    fd = g.flowdir(g.resolve_flats(g.fill_depressions(g.fill_pits(dem))))
    acc = np.asarray(g.accumulation(fd)); zz = np.asarray(dem, dtype=float)
    cel = abs(tr.a * tr.e) / 1e6
    sortida = {"p24_t100_mm": round(p24), "anys_estacio": n, "p0_mm": round(p0), "obres": []}
    for ft, cels, (x, y) in obres:
        c0, r0 = ~tr * (x, y); r0, c0 = int(r0), int(c0)
        k = int(round(30 / tr.a))                                      # la cel·la de més conca a menys de 30 m de la boca
        win = acc[max(r0 - k, 0): r0 + k + 1, max(c0 - k, 0): c0 + k + 1]
        if win.size == 0:
            continue
        dr, dc = np.unravel_index(np.argmax(win), win.shape); r, c = max(r0 - k, 0) + dr, max(c0 - k, 0) + dc
        xo, yo = tr * (c + 0.5, r + 0.5)
        mask = np.asarray(g.catchment(x=xo, y=yo, fdir=fd, xytype="coordinate"))
        a_km2 = float(mask.sum() * cel)
        if a_km2 < 0.01:
            continue
        l_km = recorregut(np.asarray(fd), mask, r, c) * tr.a / 1000
        j = (float(np.nanmax(np.where(mask, zz, np.nan))) - zz[r, c]) / max(l_km * 1000, 1)
        q, tc, cc = q_disseny(a_km2, max(l_km, 0.05), j, p24, p0)
        ft["properties"].update({"conca_km2": round(a_km2, 2), "q_disseny_m3s": round(q, 1), "tc_h": round(tc, 2)})
        xy = [tuple(round(v, 1) for v in tr * (cc_ + 0.5, rr + 0.5)) for rr, cc_ in cels]
        sortida["obres"].append({"osm": ft["properties"]["osm"], "q_m3s": round(q, 2), "xy": xy})
    return sortida


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
                    passos.append((seg, v, f"sota pont ({p['t'].get('highway') or p['t'].get('railway') or 'pont'})", p))

    # obres de drenatge suposades: riera que creua un vial o una via de tren (cap dels dos amb pont o túnel)
    rieres = [v for v in baixes if v["t"].get("waterway") in RIERES]
    vials = [v for v in baixes if v["t"].get("highway") in VIALS or v["t"].get("railway") in ("rail", "light_rail")]
    arbre_v = STRtree([v["g"] for v in vials])
    for r in rieres:
        for k in arbre_v.query(r["g"], predicate="intersects"):
            x = r["g"].intersection(vials[k]["g"])
            for pt in getattr(x, "geoms", [x]) if x.geom_type in ("Point", "MultiPoint") else []:
                s = r["g"].project(pt)
                seg = substring(r["g"], max(s - a.mig, 0), min(s + a.mig, r["g"].length))
                if seg.length > 1:
                    passos.append((seg, r, "obra de drenatge suposada", vials[k]))

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
    fets, fts, obres = 0, [], []
    for g, v, tipus, sobre in passos:
        z0, z1 = terra_boca(g, 0), terra_boca(g, 1)
        if z0 is None or z1 is None:
            continue
        rebaix, cel, baixades = 0.0, set(), []
        for s in np.arange(0, g.length + 0.1, tr.a / 3):
            rc = mostra(g.interpolate(s))
            if not rc or rc in cel or not np.isfinite(z[rc]):
                continue
            cel.add(rc)
            fons = z0 + (z1 - z0) * s / g.length
            if z[rc] > fons:
                if z[rc] - fons > 0.3: baixades.append(rc)
                rebaix = max(rebaix, float(z[rc] - fons)); z[rc] = fons
        if rebaix > 0.3:
            fets += 1
            q = LineString([a_ll(*c) for c in g.coords])
            fts.append({"type": "Feature", "geometry": mapping(q), "properties": {
                "osm": v["id"], "tipus": tipus, "via": v["t"].get("name") or v["t"].get("highway") or v["t"].get("waterway") or v["t"].get("railway"),
                "rebaix_m": round(rebaix, 1), "fons_m": [round(float(z0), 1), round(float(z1), 1)], "suposat": tipus.endswith("suposada"),
                "sobre": sobre and (sobre["t"].get("ref") or sobre["t"].get("name") or sobre["t"].get("highway") or sobre["t"].get("railway"))}})
            if v["t"].get("waterway") and not tipus.startswith("sota pont"):                                 # obra de drenatge: se'n limita el cabal
                obres.append((fts[-1], baixades, g.coords[0] if z0 >= z1 else g.coords[-1]))
    cap = capacitats(z, tr, crs, obres, conca)
    if cap:
        (d / "capacitat.json").write_text(json.dumps(cap, separators=(",", ":")))
        qs = sorted(o["q_m3s"] for o in cap["obres"])
        if qs: print(f"{len(qs)} obres de drenatge amb cabal de disseny: mediana {np.median(qs):.1f} m³/s, màxim {qs[-1]:.0f} m³/s")
    fts.sort(key=lambda f: -f["properties"]["rebaix_m"])
    (d / "passos.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": fts}, ensure_ascii=False))
    with rasterio.open(a.mdt, "w", **perfil) as f:
        f.write(z, 1)
    print(f"{len(passos)} passos candidats; {fets} obren un terraplè de més de 30 cm ({sum(f['properties']['suposat'] for f in fts)} suposats). Els més grans:")
    for f in fts[:15]:
        p = f["properties"]; c = f["geometry"]["coordinates"][len(f["geometry"]["coordinates"]) // 2]
        print(f"  {p['rebaix_m']:5.1f} m  {p['tipus']:28s} {p['via']} sota {p['sobre']}  ({c[1]:.5f}, {c[0]:.5f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
