#!/usr/bin/env python3
"""MDT de 5 m de tot Catalunya, conca per conca, i les capes de relleu de cada conca.

Fa servir les conques del tall de 100 km² (data/catalunya/conques/conques_100km2.geojson, totes, també les
petites de costa, perquè no quedin forats). Per a cada conca:

1. Baixa l'MDT de 5 m de l'IGN (WCS servicios.idee.es, cobertura Elevacion25830_5, LiDAR del PNOA) a trossos de
   5×5 km dins la caixa de la conca + 500 m, i deixa sense dades el que queda fora de la conca + 500 m.
   El desa tal com el dona l'IGN (ETRS89 UTM 30N, 5 m), en cm, GeoTIFF comprimit:  /tmp/mdt_cat/mdt5_cNNN.tif
   (el workflow el puja a la Release `mdt5-cat`). Si ja el tenim (--mdt), no el torna a baixar.
2. Mateix mètode que scripts/relleu.py (mapa de Platja d'Aro): reixa de --pas m en UTM 31N, direccions de
   drenatge i acumulació (pysheds) → rieretes a partir de --llindar-km2 de conca, amb superfície drenada i
   ordre de Strahler; HAND sobre aquesta xarxa i profunditat de les depressions → zones baixes.
   Diferència: la superfície drenada és només la de dins la conca (no suma el que entra d'aigües amunt).
3. Retalla els resultats a la conca i els escriu lleugers per a la web:
     data/catalunya/relleu/cNNN/rieretes.json     trams amb coordenades enteres (1e-5 °) i en diferències
     data/catalunya/relleu/cNNN/zones_baixes.png  PNG de paleta (5 colors), en lat/lon
     data/catalunya/relleu/cNNN/resum.json        límits del PNG, classes, km² de cada classe, trams

Ús:  python3 scripts/mdt_cat.py lots --lots 30 [--conques 1,2,3 | totes | pendents]       → matriu JSON per al workflow
     python3 scripts/mdt_cat.py conca --id 124 [--mdt fitxer.tif] [--pas 10] [--llindar-km2 0.05]
     python3 scripts/mdt_cat.py index                                     → data/catalunya/relleu/index.json"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np, requests
np.in1d = np.isin                      # pysheds encara fa servir np.in1d (numpy 2)
import rasterio
from rasterio.merge import merge
from rasterio.warp import reproject, Resampling, transform_bounds
from rasterio.transform import from_origin
from rasterio import features
from pyproj import Transformer
from shapely.geometry import shape, mapping, LineString, Point
from shapely.ops import transform as stransform
from shapely.prepared import prep

WCS = "https://servicios.idee.es/wcs-inspire/mdt"
CONQUES = Path("data/catalunya/conques/conques_100km2.geojson")
SORT = Path("data/catalunya/relleu")
TMP = Path("/tmp/mdt_cat")
MARGE = 500          # m al voltant de la conca
TROS = 5000          # m de costat de cada petició al WCS
CLASSES = {"1": "a menys de 2 m del llit", "2": "a menys d'1 m del llit", "3": "a menys de 0,5 m del llit",
           "4": "depressió sense sortida (> 0,3 m)"}
PALETA = {1: (198, 219, 239, 110), 2: (107, 174, 214, 150), 3: (33, 113, 181, 180), 4: (128, 0, 128, 200)}


def conques():
    fc = json.loads(CONQUES.read_text())
    return {f["properties"]["id"]: f for f in fc["features"]}


def peca(xa, ya, xb, yb):
    p = [("SERVICE", "WCS"), ("VERSION", "2.0.1"), ("REQUEST", "GetCoverage"), ("COVERAGEID", "Elevacion25830_5"),
         ("SUBSET", f"x({xa:.0f},{xb:.0f})"), ("SUBSET", f"y({ya:.0f},{yb:.0f})"), ("FORMAT", "image/tiff")]
    r = None
    for intent in range(6):              # el WCS de l'IGN torna 500 de tant en tant
        try:
            r = requests.get(WCS, params=p, timeout=300)
            if r.content[:2] in (b"II", b"MM"):
                return r.content
            if r.status_code == 404 or b"empty" in r.content[:500].lower():
                return None              # fora de la cobertura (mar, França)
        except requests.RequestException:
            pass
        time.sleep(10 * (intent + 1))
    raise RuntimeError(f"WCS {xa:.0f},{ya:.0f}: HTTP {getattr(r, 'status_code', '?')} {getattr(r, 'text', '')[:300]}")


def baixa_mdt(geom, nom):
    """MDT 5 m (EPSG:25830) de la conca + MARGE, sense dades fora. Torna el camí del GeoTIFF."""
    t = Transformer.from_crs(4326, 25830, always_xy=True)
    g30 = stransform(t.transform, geom).buffer(MARGE)
    x0, y0, x1, y1 = g30.bounds
    x0, y0 = np.floor(x0 / 5) * 5, np.floor(y0 / 5) * 5
    x1, y1 = np.ceil(x1 / 5) * 5, np.ceil(y1 / 5) * 5
    pg = prep(g30)
    trossos = []
    for xa in np.arange(x0, x1, TROS):
        for ya in np.arange(y0, y1, TROS):
            xb, yb = min(xa + TROS, x1), min(ya + TROS, y1)
            if xb - xa < 50 or yb - ya < 50:     # tros massa estret: l'amplia (el WCS en torna de buits)
                xa, ya = min(xa, xb - 50), min(ya, yb - 50)
            from shapely.geometry import box
            if not pg.intersects(box(xa, ya, xb, yb)):
                continue
            c = peca(xa, ya, xb, yb)
            if c is None:
                continue
            f = TMP / f"{nom}_{len(trossos)}.tif"
            f.write_bytes(c)
            r = rasterio.open(f)
            if r.width < 2 or r.height < 2 or not all(np.isfinite(r.bounds)):
                print(f"  tros {xa:.0f},{ya:.0f} buit o estrany ({r.width}×{r.height}), el salto")
                r.close(); continue
            trossos.append(r)
    if not trossos:
        raise RuntimeError("cap tros amb dades")
    nd = trossos[0].nodata
    z, tr = merge(trossos, bounds=(x0, y0, x1, y1), res=5, nodata=nd if nd is not None else -9999)
    z = z[0].astype("float32")
    if nd is not None:
        z[z == nd] = np.nan
    z[z < -50] = np.nan
    fora = features.geometry_mask([mapping(g30)], z.shape, tr, invert=False)
    z[fora] = np.nan
    for s in trossos:
        s.close()
        Path(s.name).unlink()
    cm = np.where(np.isnan(z), -32768 * 100, np.round(z * 100)).astype("int32")
    out = TMP / f"mdt5_{nom}.tif"
    with rasterio.open(out, "w", driver="GTiff", height=z.shape[0], width=z.shape[1], count=1, dtype="int32",
                       crs="EPSG:25830", transform=tr, nodata=-32768 * 100, compress="deflate", predictor=2, zlevel=9,
                       tiled=True, blockxsize=512, blockysize=512) as f:
        f.write(cm, 1)
        f.update_tags(unitats="cm", font="IGN MDT05 (WCS Elevacion25830_5)", conca=nom)
    print(f"MDT 5 m: {z.shape[1]}×{z.shape[0]} cel·les, {len(trossos)} trossos, "
          f"{np.isfinite(z).sum() * 25 / 1e6:.1f} km² amb dades, {out.stat().st_size / 1e6:.1f} MB")
    return out


def llegeix_mdt(cami):
    with rasterio.open(cami) as f:
        z = f.read(1).astype("float32")
        nd = f.nodata
        z[z == nd] = np.nan
        return z / 100, f.transform, f.crs


def a_utm31(z, tr, crs, pas, caixa):
    dst = "EPSG:25831"
    w, s, e, n = transform_bounds(crs, dst, *caixa)
    w, s = np.floor(w / pas) * pas, np.floor(s / pas) * pas
    ncol, nfil = int(np.ceil((e - w) / pas)), int(np.ceil((n - s) / pas))
    out = np.full((nfil, ncol), np.nan, dtype="float32")
    dtr = from_origin(w, s + nfil * pas, pas, pas)
    reproject(z, out, src_transform=tr, src_crs=crs, src_nodata=np.nan, dst_transform=dtr, dst_crs=dst,
              dst_nodata=np.nan, resampling=Resampling.average)
    return out, dtr, dst


def codifica(coords):
    """[[lon,lat],...] → [x0,y0,dx1,dy1,...] enters a 1e-5 °."""
    q = np.round(np.asarray(coords) * 1e5).astype(np.int64)
    d = np.vstack([q[:1], np.diff(q, axis=0)])
    d = d[(d != 0).any(axis=1) | (np.arange(len(d)) == 0)]
    return d.ravel().tolist()


def relleu(cami_mdt, geom, nom, pas, llindar):
    z5, tr5, crs5 = llegeix_mdt(cami_mdt)
    h, w = z5.shape
    caixa = (tr5.c, tr5.f + h * tr5.e, tr5.c + w * tr5.a, tr5.f)
    z, tr, crs = a_utm31(z5, tr5, crs5, pas, caixa)
    del z5
    mar = np.isnan(z) | (z <= 0.2)
    perfil = dict(driver="GTiff", height=z.shape[0], width=z.shape[1], count=1, dtype="float32", crs=crs, transform=tr,
                  nodata=np.nan, compress="deflate")
    f10 = TMP / f"mdt{int(pas)}_{nom}.tif"
    with rasterio.open(f10, "w", **perfil) as f:
        f.write(np.where(mar, np.nan, z), 1)
    print(f"Reixa de {pas:.0f} m: {z.shape[1]}×{z.shape[0]}")

    from pysheds.grid import Grid
    g = Grid.from_raster(str(f10))
    dem = g.read_raster(str(f10))
    omplert = g.fill_depressions(g.fill_pits(dem))
    plana = g.resolve_flats(omplert)
    fd = g.flowdir(plana)
    acc = g.accumulation(fd)
    cel_km2 = pas ** 2 / 1e6
    xarxa = acc > llindar / cel_km2
    ordre = g.stream_order(fd, xarxa)
    rius = g.extract_river_network(fd, xarxa)

    t31 = Transformer.from_crs(4326, crs, always_xy=True)
    g31 = stransform(t31.transform, geom)
    dins = prep(g31)
    a4326 = Transformer.from_crs(crs, 4326, always_xy=True)
    trams = []
    for f in rius["features"]:
        coords = f["geometry"]["coordinates"]
        if len(coords) < 2:
            continue
        ln = LineString(coords)
        if not dins.contains(ln.interpolate(0.5, normalized=True)):
            continue                      # el tram és d'una conca veïna
        cel = [(int((tr.f - y) // pas), int((x - tr.c) // pas)) for x, y in coords]
        cel = [(r, c) for r, c in cel if 0 <= r < acc.shape[0] and 0 <= c < acc.shape[1]]
        area = max(float(acc[r, c]) for r, c in cel) * cel_km2
        o = int(max(ordre[r, c] for r, c in cel))
        sim = ln.simplify(pas / 2)
        ll = [a4326.transform(x, y) for x, y in sim.coords]
        trams.append([o, int(round(area * 100))] + codifica(ll))
    trams.sort(key=lambda t: -t[1])

    hand = np.asarray(g.compute_hand(fd, plana, xarxa), dtype="float32")
    prof = np.asarray(omplert, dtype="float32") - np.asarray(dem, dtype="float32")
    classe = np.zeros(z.shape, dtype=np.uint8)
    classe[(hand < 2)] = 1
    classe[(hand < 1)] = 2
    classe[(hand < 0.5)] = 3
    classe[prof > 0.3] = 4                      # depressió: l'aigua s'hi queda
    classe[mar] = 0
    classe[features.geometry_mask([mapping(g31)], z.shape, tr, invert=False)] = 0
    km2 = {k: round(float((classe == int(k)).sum() * cel_km2), 2) for k in CLASSES}

    # PNG de paleta en lat/lon per superposar-lo amb Leaflet (cel·les ~quadrades)
    w, s, e, n = geom.bounds
    lat = (s + n) / 2
    pas_y = pas / 111320
    pas_x = pas / (111320 * np.cos(np.radians(lat)))
    ncol, nfil = int(np.ceil((e - w) / pas_x)), int(np.ceil((n - s) / pas_y))
    cl = np.zeros((nfil, ncol), dtype=np.uint8)
    reproject(classe, cl, src_transform=tr, src_crs=crs, dst_transform=rasterio.Affine(pas_x, 0, w, 0, -pas_y, n),
              dst_crs="EPSG:4326", resampling=Resampling.nearest)
    from PIL import Image
    d = SORT / nom
    d.mkdir(parents=True, exist_ok=True)
    im = Image.fromarray(cl, "P")
    pal = [0, 0, 0] * 256
    alfa = [0] * 256
    for k, c in PALETA.items():
        pal[3 * k:3 * k + 3] = c[:3]
        alfa[k] = c[3]
    im.putpalette(pal)
    im.save(d / "zones_baixes.png", optimize=True, transparency=bytes(alfa[:5]))
    (d / "rieretes.json").write_text(json.dumps({"escala": 1e5, "camps": ["ordre", "area_km2*100", "x0", "y0", "dx", "dy", "..."],
                                                 "trams": trams}, separators=(",", ":")))
    resum = {"conca": nom, "limits": [[n - nfil * pas_y, w], [n, w + ncol * pas_x]], "pas_m": pas,
             "llindar_km2": llindar, "trams": len(trams), "classes": CLASSES, "km2": km2,
             "font": "MDT 5 m IGN (PNOA LiDAR)", "calculat": time.strftime("%Y-%m-%d")}
    (d / "resum.json").write_text(json.dumps(resum, ensure_ascii=False, separators=(",", ":")))
    print(f"rieretes: {len(trams)} trams (conca > {llindar} km²), {(d / 'rieretes.json').stat().st_size / 1e3:.0f} kB; "
          f"zones baixes {(d / 'zones_baixes.png').stat().st_size / 1e3:.0f} kB: "
          + ", ".join(f"{CLASSES[k].split(' (')[0]} {v} km²" for k, v in km2.items()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="ordre", required=True)
    a1 = sub.add_parser("lots"); a1.add_argument("--lots", type=int, default=30); a1.add_argument("--conques", default="totes")
    a2 = sub.add_parser("conca"); a2.add_argument("--id", type=int, required=True); a2.add_argument("--mdt")
    a2.add_argument("--pas", type=float, default=10); a2.add_argument("--llindar-km2", type=float, default=0.05)
    a2.add_argument("--nomes-mdt", action="store_true")
    sub.add_parser("index")
    a = ap.parse_args()
    C = conques()

    if a.ordre == "lots":
        if a.conques == "pendents":       # les que encara no tenen relleu desat
            a.conques = ",".join(str(i) for i in sorted(C) if not (SORT / f"c{i:03d}" / "resum.json").exists()) or "totes"
        ids = sorted(C) if a.conques == "totes" else [int(x) for x in a.conques.split(",") if x.strip()]
        # repartiment per superfície (les grans tarden més): la més gran al lot amb menys feina
        lots = [[] for _ in range(min(a.lots, len(ids)))]
        pes = [0.0] * len(lots)
        for i in sorted(ids, key=lambda i: -C[i]["properties"]["km2"]):
            k = int(np.argmin(pes)); lots[k].append(i); pes[k] += C[i]["properties"]["km2"] + 20
        print(json.dumps({"lot": [" ".join(map(str, l)) for l in lots]}))
        return 0

    if a.ordre == "index":
        idx = {}
        for r in sorted(SORT.glob("c*/resum.json")):
            j = json.loads(r.read_text())
            idx[j["conca"]] = {"limits": j["limits"], "trams": j["trams"], "km2": j["km2"]}
        (SORT / "index.json").write_text(json.dumps({"conques": idx, "classes": CLASSES}, ensure_ascii=False, separators=(",", ":")))
        print(f"index: {len(idx)} conques de {len(C)}")
        return 0

    TMP.mkdir(parents=True, exist_ok=True)
    nom = f"c{a.id:03d}"
    geom = shape(C[a.id]["geometry"])
    print(f"== Conca {a.id}: {C[a.id]['properties']['km2']} km²")
    cami = Path(a.mdt) if a.mdt and Path(a.mdt).exists() else baixa_mdt(geom, nom)
    if not a.nomes_mdt:
        relleu(cami, geom, nom, a.pas, a.llindar_km2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
