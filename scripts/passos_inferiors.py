#!/usr/bin/env python3
"""Passos inferiors d'un poble: carrers i camins que passen per sota d'una via de tren o d'una carretera.

Són els punts baixos on l'aigua de la pluja s'acumula primer (el clavegueram no se l'empassa i no té sortida per superfície).

1. Candidats (OpenStreetMap, data/<conca>/passos/osm.json, la mateixa consulta que scripts/passos.py):
   - vial (highway) amb tunnel=* o layer negatiu que passa per sota d'una altra via (tren o carretera);
   - vial a nivell creuat per un pont (bridge=*) de tren o de carretera.
   Es descarten els vials que només passen sota un edifici (building_passage) i els que no existeixen encara.
   Els trams que són a menys de --agrupa m l'un de l'altre sota la mateixa via són un sol pas (dues calçades, vorera...).
2. Relleu (MDT del LiDAR, 2 m de l'ICGC si s'hi pot arribar, si no 5 m de l'IGN), amb --mdt: a cada pas
   - perfil del vial: cota mínima sota la via i cota de cada rampa a --rampa m; profunditat = rampa més baixa − mínim
     (el que s'ha d'omplir abans que l'aigua surti carrer avall);
   - clotada: relleu omplert − terreny (on l'aigua queda sense sortida per superfície), volum que hi cap i cota del sobreeixidor;
   - conca local: àrea que hi drena (D8 sobre el relleu omplert, fins a la primera clotada ≥ 0,3 m que troba l'aigua);
   - pluja per omplir-lo (mm) = volum / conca local (tota la pluja escolant, sense clavegueram): com més baixa, abans s'inunda.
   Si el pas és un túnel tapat al MDT (el LiDAR veu el terraplè, no el vial), el fons es fa anar en línia recta entre les boques.

Sortides: data/<conca>/passos/passos_inferiors.geojson i .csv, i la taula inundacions.punts_baixos de Supabase (si hi ha clau).

Ús:  python3 scripts/passos_inferiors.py --config config/granollers.json [--mdt /tmp/passos/mdt.tif] [--supabase]"""
import argparse, csv, json, os, sys
from pathlib import Path
import numpy as np
from pyproj import Transformer
from shapely.geometry import LineString, Point, mapping
from shapely.ops import substring, unary_union
from shapely.strtree import STRtree

VIALS_NO = ("construction", "proposed", "platform", "elevator", "corridor", "bus_stop", "raceway")
TRENS = ("rail", "light_rail", "subway", "narrow_gauge")
CARREGA = {"motorway": "autopista", "trunk": "autovia", "primary": "carretera", "secondary": "carretera", "tertiary": "carretera",
           "motorway_link": "enllaç", "trunk_link": "enllaç", "primary_link": "enllaç", "secondary_link": "enllaç", "tertiary_link": "enllaç"}
PEATONS = ("footway", "path", "pedestrian", "cycleway", "steps", "bridleway")


def classe(hs):
    if any(h not in PEATONS + ("track", "service") for h in hs):
        return "carrer o carretera"
    if any(h == "service" for h in hs):
        return "accés"
    return "camí" if any(h == "track" for h in hs) else "vianants"


def nom_via(t):
    if t.get("railway") in TRENS:
        return "tren" + (f" ({t['name']})" if t.get("name") else "")
    return t.get("ref") or t.get("name") or CARREGA.get(t.get("highway"), t.get("highway") or "via")


def candidats(osm, caixa, a_utm, agrupa):
    vies = []
    w, s, e, n = caixa
    for el in osm.get("elements", []):
        t = el.get("tags", {})
        if el.get("type") != "way" or "geometry" not in el:
            continue
        if "highway" not in t and t.get("railway") not in TRENS:
            continue
        if t.get("highway") in ("motorway", "motorway_link") and (t.get("tunnel", "no") != "no" or t.get("bridge", "no") == "no"):
            continue                                                  # l'autopista per sota no és un punt del poble
        if t.get("highway") in VIALS_NO or t.get("railway") in ("construction", "proposed", "abandoned", "razed", "disused") \
                or "construction" in t or "proposed" in t or t.get("area") == "yes":
            continue
        ll = [(p["lon"], p["lat"]) for p in el["geometry"]]
        if len(ll) < 2 or not any(w <= x <= e and s <= y <= n for x, y in ll):
            continue
        vies.append({"id": el["id"], "t": t, "g": LineString([a_utm(*p) for p in ll])})
    rang = lambda o: (o["t"].get("railway") in TRENS, o["t"].get("highway") in ("motorway", "trunk"), o["t"].get("highway") in CARREGA,
                      o["t"].get("highway") not in ("service", "track") + PEATONS)
    tunel = lambda v: v["t"].get("tunnel", "no") not in ("no", "building_passage") or v["t"].get("layer", "0").startswith("-")
    pont = lambda v: v["t"].get("bridge", "no") != "no"
    dalt = [v for v in vies if not tunel(v)]
    arbre = STRtree([v["g"] for v in dalt])
    cands = []
    # a) vial en túnel (o layer < 0) que passa per sota d'una via
    for v in vies:
        if "highway" not in v["t"] or not tunel(v):
            continue
        sobre = [dalt[k] for k in arbre.query(v["g"], predicate="crosses")]
        sobre = [o for o in sobre if o["t"].get("railway") in TRENS or o["t"].get("highway") not in PEATONS + (None,)]
        if sobre and v["g"].length <= 400:
            o = max(sobre, key=rang)
            cands.append({"g": v["g"], "v": v, "sobre": o, "pt": v["g"].intersection(o["g"]), "tipus": "túnel"})
    # b) vial a nivell creuat per un pont de tren o de carretera
    baixos = [v for v in vies if "highway" in v["t"] and not pont(v) and not tunel(v)]
    arbre_b = STRtree([v["g"] for v in baixos])
    for p in vies:
        if not pont(p) or not (p["t"].get("railway") in TRENS or p["t"].get("highway") not in PEATONS):
            continue
        for k in arbre_b.query(p["g"], predicate="intersects"):
            v = baixos[k]
            if v["g"].touches(p["g"]) or v["t"].get("layer", "0") == p["t"].get("layer", "1"):
                continue
            x = v["g"].intersection(p["g"])
            for pt in getattr(x, "geoms", [x]) if x.geom_type in ("Point", "MultiPoint") else []:
                cands.append({"g": v["g"], "v": v, "sobre": p, "pt": pt, "tipus": "sota pont"})
    # agrupa: mateix pas si són a menys de --agrupa m i sota la mateixa via (o el mateix nom de via de dalt)
    passos = []
    for c in sorted(cands, key=lambda c: (c["v"]["t"].get("highway") in PEATONS, -len(c["v"]["t"].get("name", "")))):
        for p in passos:
            if p["pt"].distance(c["pt"].centroid) < agrupa:
                p["membres"].append(c)
                if rang(c["sobre"]) > rang(p["sobre"]):
                    p["sobre"] = c["sobre"]
                break
        else:
            passos.append({"pt": c["pt"].centroid, "sobre": c["sobre"], "membres": [c]})
    return passos, vies


class Mdt:
    def __init__(self, f):
        import rasterio
        with rasterio.open(f) as r:
            self.z = r.read(1).astype("float64"); self.tr = r.transform; self.crs = r.crs
        self.z[(self.z < -100) | ~np.isfinite(self.z)] = np.nan
        self.pas = self.tr.a

    def rc(self, x, y):
        c, r = ~self.tr * (x, y)
        r, c = int(r), int(c)
        return (r, c) if 0 <= r < self.z.shape[0] and 0 <= c < self.z.shape[1] else None

    def cota(self, x, y):
        rc = self.rc(x, y)
        return self.z[rc] if rc else np.nan


def omple(z):
    """Relleu omplert (priority-flood amb skimage, reconstrucció per erosió des de la vora)."""
    from skimage.morphology import reconstruction
    zz = np.where(np.isfinite(z), z, np.nanmin(z) - 1)
    llavor = zz.copy(); llavor[1:-1, 1:-1] = zz.max()
    return reconstruction(llavor, zz, method="erosion")


def receptors(zf):
    """D8 sobre el relleu omplert amb un pendent mínim (epsilon) per travessar les parts planes."""
    nf, nc = zf.shape
    # pendent mínim: allisa les planes amb la distància a la vora més baixa (aproximació: petit pendent cap a l'est-sud)
    N = nf * nc
    best = np.full(N, -1); bdz = np.zeros(N)
    zr = zf.ravel()
    rr, cc = np.divmod(np.arange(N), nc)
    for dr, dc in ((-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1)):
        r2, c2 = rr + dr, cc + dc
        ok = (r2 >= 0) & (r2 < nf) & (c2 >= 0) & (c2 < nc)
        j = np.where(ok, r2 * nc + c2, 0)
        dz = np.where(ok, (zr - zr[j]) / np.hypot(dr, dc), -np.inf)
        m = dz > bdz
        best[m] = j[m]; bdz[m] = dz[m]
    return best


def primera_clotada(best, etiq):
    """Per a cada cel·la, la primera clotada (etiqueta ≥ 0) que troba l'aigua baixant (pointer jumping)."""
    N = best.size
    salt = np.where(best >= 0, best, np.arange(N))
    f = etiq.ravel().copy()
    for _ in range(40):
        f = np.where(f >= 0, f, f[salt])
        nou = salt[salt]
        if np.array_equal(nou, salt):
            break
        salt = nou
    return f


def analitza(passos, vies, mdt, a_ll, rampa, clot_min):
    from scipy import ndimage
    z = mdt.z.copy()
    # tapa els túnels que el LiDAR no veu: fons en línia recta entre boques
    for p in passos:
        for c in p["membres"]:
            if c["tipus"] != "túnel":
                continue
            g = c["g"]; z0, z1 = mdt.cota(*g.coords[0]), mdt.cota(*g.coords[-1])
            if not (np.isfinite(z0) and np.isfinite(z1)):
                continue
            for s in np.arange(0, g.length + 0.1, mdt.pas / 2):
                pt = g.interpolate(s); rc = mdt.rc(pt.x, pt.y)
                if rc and np.isfinite(z[rc]):
                    z[rc] = min(z[rc], z0 + (z1 - z0) * s / g.length)
    zf = omple(z)
    prof = np.where(np.isfinite(z), zf - z, 0)
    etiq, nclot = ndimage.label(prof >= clot_min, structure=np.ones((3, 3)))
    etiq = etiq - 1                                                  # -1 = fora de clotada
    best = receptors(np.where(np.isfinite(z), zf, -1e9) + 1e-4 * _distancia_sortida(zf))
    dren = primera_clotada(best, etiq)
    area_clot = np.bincount(dren[dren >= 0], minlength=nclot) * mdt.pas ** 2
    vol_clot = ndimage.sum(prof, etiq + 1, index=np.arange(1, nclot + 1)) * mdt.pas ** 2
    prof_clot = ndimage.maximum(prof, etiq + 1, index=np.arange(1, nclot + 1))
    vora_clot = ndimage.maximum(zf, etiq + 1, index=np.arange(1, nclot + 1))
    # xarxa de vials per al perfil: unió de tots els vials (per seguir el carrer més enllà del tram de l'OSM)
    vials = [v for v in vies if "highway" in v["t"]]
    arbre = STRtree([v["g"] for v in vials])
    out = []
    for i, p in enumerate(passos):
        pt = p["pt"]
        # perfil al llarg de cada membre i dels vials que hi enllacen, fins a --rampa m del creuament
        carrers = unary_union([vials[k]["g"] for k in arbre.query(pt.buffer(rampa))] + [c["g"] for c in p["membres"]])
        fons, rampes = [], []
        for c in p["membres"]:
            g = c["g"]; s0 = g.project(c["pt"].centroid)
            for s in np.arange(max(s0 - 15, 0), min(s0 + 15, g.length) + 0.1, mdt.pas / 2):
                q = g.interpolate(s); fons.append(z[mdt.rc(q.x, q.y)] if mdt.rc(q.x, q.y) else np.nan)
            for sentit in (-1, 1):
                vals = []
                for d in np.arange(rampa * 0.6, rampa + 0.1, mdt.pas):
                    s = s0 + sentit * d
                    if 0 <= s <= g.length:
                        q = g.interpolate(s); vals.append(mdt.cota(q.x, q.y))
                    else:                                             # el tram s'acaba: segueix pels carrers que hi enllacen
                        q = g.interpolate(g.length if sentit > 0 else 0)
                        anell = pt.buffer(d).exterior.intersection(carrers)
                        pts = [anell] if anell.geom_type == "Point" else list(getattr(anell, "geoms", []))
                        pts = [x for x in pts if x.geom_type == "Point" and x.distance(q) < d]
                        if pts:
                            vals.append(mdt.cota(*min(pts, key=lambda x: x.distance(q)).coords[0]))
                vals = [v for v in vals if np.isfinite(v)]
                if vals:
                    rampes.append(float(np.median(vals)))
        fons = [f for f in fons if np.isfinite(f)]
        if not fons:
            continue
        zmin = float(min(fons))
        # clotada a menys de 15 m del punt més baix del pas
        rc = mdt.rc(pt.x, pt.y); k = max(1, int(15 / mdt.pas))
        win = (slice(max(rc[0] - k, 0), rc[0] + k + 1), slice(max(rc[1] - k, 0), rc[1] + k + 1)) if rc else None
        e = etiq[win]; pr = prof[win]
        if win is not None and (e >= 0).any():
            j = int(e.ravel()[np.argmax(np.where(e >= 0, pr, -1))])
            clot = {"clotada_m": round(float(prof_clot[j]), 2), "volum_m3": round(float(vol_clot[j])),
                    "conca_local_ha": round(float(area_clot[j]) / 1e4, 2), "cota_sobreeix_m": round(float(vora_clot[j]), 2)}
            clot["pluja_omplir_mm"] = round(1000 * clot["volum_m3"] / max(area_clot[j], 1), 1)
        else:
            clot = {"clotada_m": 0.0, "volum_m3": 0, "conca_local_ha": None, "cota_sobreeix_m": None, "pluja_omplir_mm": None}
        m0 = p["membres"][0]; t = m0["v"]["t"]
        lon, lat = a_ll(pt.x, pt.y)
        noms = sorted({c["v"]["t"].get("name") for c in p["membres"] if c["v"]["t"].get("name")})
        out.append({
            "id": f"pi{i:03d}", "nom": (", ".join(noms) or t.get("highway")) + " sota " + nom_via(p["sobre"]["t"]),
            "vial": ", ".join(noms) or t.get("highway"),
            "classe": classe([c["v"]["t"].get("highway") for c in p["membres"]]),
            "sota": nom_via(p["sobre"]["t"]), "tipus": m0["tipus"], "osm": sorted({c["v"]["id"] for c in p["membres"]}),
            "lon": round(lon, 6), "lat": round(lat, 6), "cota_fons_m": round(zmin, 2),
            "profunditat_m": round(max(min(rampes) - zmin, 0), 2) if len(rampes) >= 2 else None, **clot})
    return out


def _distancia_sortida(zf):
    """Distància (en cel·les) dins de cada plana fins a la vora per on pot sortir: petit pendent perquè les planes drenin."""
    from scipy import ndimage
    pla = np.zeros(zf.shape, bool)
    pla[1:-1, 1:-1] = True
    for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)):
        pla[1:-1, 1:-1] &= zf[1 + dr: zf.shape[0] - 1 + dr, 1 + dc: zf.shape[1] - 1 + dc] >= zf[1:-1, 1:-1] - 1e-6
    return ndimage.distance_transform_edt(pla)


WCS = [  # (url, descripció); el primer que respon amb una cobertura de 2 m o, si no, de 5 m
    ("https://geoserveis.icgc.cat/servei/catalunya/mdt2m/wcs", "ICGC"),
    ("https://geoserveis.icgc.cat/servei/catalunya/mdt/wcs", "ICGC"),
    ("https://geoserveis.icgc.cat/icgc_mdt2m/wcs/service", "ICGC"),
    ("https://servicios.idee.es/wcs-inspire/mdt", "IGN"),
]


def cobertures(url):
    import re, requests
    try:
        r = requests.get(url, params={"SERVICE": "WCS", "REQUEST": "GetCapabilities"}, timeout=(15, 120))
        ids = re.findall(r"<(?:wcs:)?CoverageId>([^<]+)<", r.text)
        print(f"{url}: HTTP {r.status_code}, {len(ids)} cobertures {ids[:40]}")
        return ids
    except Exception as ex:
        print(f"{url}: {ex}"); return []


def baixa_mdt(caixa, sortida, marge=1500):
    """MDT del LiDAR a la caixa (+ marge), a la resolució més fina que es trobi, en UTM 31N."""
    import re, time, requests, rasterio
    from rasterio.merge import merge
    from rasterio.warp import reproject, Resampling, transform_bounds
    from rasterio.transform import from_origin
    tria = None
    for url, font in WCS:
        ids = cobertures(url)
        for pref in ("2", "5"):
            m = [i for i in ids if re.search(rf"(^|[^0-9]){pref}(m|x{pref}|_?m?$)", i, re.I) and ("25831" in i or font == "ICGC")] \
                or [i for i in ids if re.search(rf"_{pref}$", i) and "25830" in i]
            if m:
                tria = (url, font, m[0], float(pref)); break
        if tria and tria[3] == 2:
            break
    if not tria:
        tria = (WCS[-1][0], "IGN", "Elevacion25830_5", 5.0)
    url, font, cov, pas = tria
    epsg = 25830 if "25830" in cov else 25831
    print(f"MDT: {font} {cov} ({pas:g} m)")
    w, s, e, n = transform_bounds("EPSG:4326", f"EPSG:{epsg}", *caixa)
    w, s, e, n = w - marge, s - marge, e + marge, n + marge
    tmp = Path(sortida).parent; tmp.mkdir(parents=True, exist_ok=True)
    mida = 1000 * pas; trossos = []
    for xa in np.arange(w, e, mida):
        for ya in np.arange(s, n, mida):
            p = [("SERVICE", "WCS"), ("VERSION", "2.0.1"), ("REQUEST", "GetCoverage"), ("COVERAGEID", cov),
                 ("SUBSET", f"x({xa:.0f},{min(xa + mida, e):.0f})"), ("SUBSET", f"y({ya:.0f},{min(ya + mida, n):.0f})"),
                 ("FORMAT", "image/tiff")]
            for intent in range(5):
                try:
                    r = requests.get(url, params=p, timeout=180)
                    if r.content[:2] in (b"II", b"MM"):
                        break
                except requests.RequestException:
                    pass
                time.sleep(10 * (intent + 1))
            if r.content[:2] not in (b"II", b"MM"):
                sys.exit(f"WCS {xa:.0f},{ya:.0f}: HTTP {r.status_code} {r.text[:300]}")
            f = tmp / f"tros_{len(trossos)}.tif"; f.write_bytes(r.content); trossos.append(rasterio.open(f))
    mos, tr = merge(trossos)
    v = mos[0][np.isfinite(mos[0]) & (mos[0] > -100)]
    print(f"Trossos: tipus {trossos[0].dtypes[0]}, {np.mean(np.abs(v - np.round(v)) < 1e-3):.0%} de cotes enteres")
    W, S, E, N = transform_bounds("EPSG:4326", "EPSG:25831", *caixa)
    W, S, E, N = W - marge, S - marge, E + marge, N + marge
    nc, nf = int((E - W) // pas), int((N - S) // pas)
    out = np.full((nf, nc), np.nan, "float32"); dtr = from_origin(W, N, pas, pas)
    reproject(mos[0].astype("float32"), out, src_transform=tr, src_crs=trossos[0].crs, src_nodata=trossos[0].nodata,
              dst_transform=dtr, dst_crs="EPSG:25831", dst_nodata=np.nan, resampling=Resampling.bilinear)
    with rasterio.open(sortida, "w", driver="GTiff", height=nf, width=nc, count=1, dtype="float32", crs="EPSG:25831",
                       transform=dtr, nodata=np.nan, compress="deflate") as f:
        f.write(out, 1)
    print(f"MDT {font} {pas:g} m: {nc}×{nf} cel·les → {sortida}")
    return font


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--mdt", help="MDT en UTM (GeoTIFF); sense, només llista els candidats")
    ap.add_argument("--agrupa", type=float, default=60)
    ap.add_argument("--rampa", type=float, default=60, help="distància del creuament on es mira la cota de cada rampa (m)")
    ap.add_argument("--clotada-min", type=float, default=0.3)
    ap.add_argument("--supabase", action="store_true")
    ap.add_argument("--baixa-mdt", action="store_true", help="baixa l'MDT a --mdt abans d'analitzar")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    conca = cfg["conca"]; caixa = cfg["passos_inferiors"]["caixa"]
    d = Path("data") / conca / "passos"
    osm = json.loads((d / "osm.json").read_text())
    crs = "EPSG:25831"
    a_utm = Transformer.from_crs(4326, crs, always_xy=True).transform
    a_ll = Transformer.from_crs(crs, 4326, always_xy=True).transform
    passos, vies = candidats(osm, caixa, a_utm, a.agrupa)
    print(f"{len(passos)} passos inferiors a la caixa {caixa}")
    if not a.mdt:
        for p in passos:
            lon, lat = a_ll(p["pt"].x, p["pt"].y)
            ns = sorted({c["v"]["t"].get("name") or c["v"]["t"].get("highway") for c in p["membres"]})
            print(f"  {lat:.5f},{lon:.5f}  {', '.join(ns)} sota {nom_via(p['sobre']['t'])} [{p['membres'][0]['tipus']}]")
        return 0
    font = baixa_mdt(caixa, a.mdt) if a.baixa_mdt else "?"
    mdt = Mdt(a.mdt)
    rs = analitza(passos, vies, mdt, a_ll, a.rampa, a.clotada_min)
    rs.sort(key=lambda r: (r["pluja_omplir_mm"] is None, r["pluja_omplir_mm"] or 0))
    for k, r in enumerate(rs):
        r["ordre"] = k + 1
    camps = ["ordre", "id", "nom", "vial", "classe", "sota", "tipus", "lat", "lon", "cota_fons_m", "profunditat_m", "clotada_m",
             "volum_m3", "conca_local_ha", "pluja_omplir_mm", "cota_sobreeix_m", "osm"]
    with open(d / "passos_inferiors.csv", "w", newline="") as f:
        w = csv.DictWriter(f, camps); w.writeheader()
        for r in rs:
            w.writerow({**r, "osm": " ".join(map(str, r["osm"]))})
    fts = [{"type": "Feature", "geometry": {"type": "Point", "coordinates": [r["lon"], r["lat"]]},
            "properties": {k: r[k] for k in camps if k not in ("lat", "lon")}} for r in rs]
    (d / "passos_inferiors.geojson").write_text(json.dumps({"type": "FeatureCollection", "mdt": f"{font} {mdt.pas:g} m", "features": fts}, ensure_ascii=False))
    print(f"MDT {mdt.pas:g} m. {len(rs)} passos; ordenats per la pluja que cal per omplir-los:")
    for r in rs:
        print(f"  {r['ordre']:3d} {r['nom'][:60]:60s} prof {r['profunditat_m']} m clotada {r['clotada_m']} m "
              f"vol {r['volum_m3']} m³ conca {r['conca_local_ha']} ha → {r['pluja_omplir_mm']} mm")
    if a.supabase and os.environ.get("SUPABASE_SECRET_KEY"):
        sys.path.insert(0, str(Path(__file__).parent))
        import sb
        sb.esborra("punts_baixos", {"conca": f"eq.{conca}", "tipus": "in.(pas inferior)"})
        n = sb.insereix("punts_baixos", [{"conca": conca, "id": r["id"], "tipus": "pas inferior", "nom": r["nom"],
                                          "lon": r["lon"], "lat": r["lat"], "dades": {k: r[k] for k in camps if k not in ("id", "nom", "lat", "lon")}}
                                         for r in rs], actualitza=True)
        print(f"Supabase: {n} punts a inundacions.punts_baixos")
    return 0


if __name__ == "__main__":
    sys.exit(main())
