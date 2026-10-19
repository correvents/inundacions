#!/usr/bin/env python3
"""Relleu d'alta resolució: totes les rieretes i les zones baixes que drenen malament.

1. Baixa l'MDT de 5 m de l'IGN (WCS servicios.idee.es, cobertura Elevacion25830_5, del LiDAR del PNOA) a trossos,
   el mosaica i el passa a una reixa de --pas metres en ETRS89 UTM 31N.
2. Direccions de drenatge i acumulació (pysheds) → xarxa de rieretes a partir de --llindar-km2 de conca,
   amb la superfície drenada i l'ordre de Strahler de cada tram.
3. HAND (alçada sobre el drenatge més proper) i profunditat de les depressions: on l'aigua no té sortida
   o queda a menys d'1 m del llit, és on s'embassa amb tempestes fortes encara que no hi passi cap riera.

Escriu:
  data/<conca>/rieretes.geojson
  data/<conca>/zones_baixes.png i zones_baixes.json (límits per superposar-la al mapa)
  /tmp/relleu/mdt.tif (la reixa, per a la simulació 2D; no va al repo)

Ús:  python3 scripts/relleu.py --config config/ridaura.json [--pas 10] [--llindar-km2 0.05]"""
import argparse, json, sys
from pathlib import Path
import numpy as np, requests
np.in1d = np.isin                      # pysheds encara fa servir np.in1d (numpy 2)
import rasterio
from rasterio.merge import merge
from rasterio.warp import reproject, Resampling, transform_bounds
from rasterio.transform import from_origin
from pyproj import Transformer
from shapely.geometry import shape, mapping
from shapely.ops import unary_union

WCS = "https://servicios.idee.es/wcs-inspire/mdt"
CAIXA = [2.86, 41.75, 3.09, 41.89]      # conca del Ridaura + plana de Platja d'Aro
TMP = Path("/tmp/relleu")


def baixa_mdt():
    t = Transformer.from_crs(4326, 25830, always_xy=True)
    x0, y0 = t.transform(CAIXA[0], CAIXA[1]); x1, y1 = t.transform(CAIXA[2], CAIXA[3])
    x0, x1 = min(x0, x1) - 500, max(x0, x1) + 500
    y0, y1 = min(y0, y1) - 500, max(y0, y1) + 500
    trossos, mida = [], 5000
    for xa in np.arange(x0, x1, mida):
        for ya in np.arange(y0, y1, mida):
            p = [("SERVICE", "WCS"), ("VERSION", "2.0.1"), ("REQUEST", "GetCoverage"), ("COVERAGEID", "Elevacion25830_5"),
                 ("SUBSET", f"x({xa:.0f},{min(xa + mida, x1):.0f})"), ("SUBSET", f"y({ya:.0f},{min(ya + mida, y1):.0f})"),
                 ("FORMAT", "image/tiff")]
            r = requests.get(WCS, params=p, timeout=180)
            if r.content[:2] not in (b"II", b"MM"):
                sys.exit(f"WCS {xa:.0f},{ya:.0f}: HTTP {r.status_code} {r.text[:300]}")
            fitxer = TMP / f"mdt5_{len(trossos)}.tif"
            fitxer.write_bytes(r.content)
            trossos.append(rasterio.open(fitxer))
    mosaic, tr = merge(trossos)
    print(f"MDT 5 m: {mosaic.shape[1]}×{mosaic.shape[2]} cel·les, {len(trossos)} trossos")
    return mosaic[0].astype("float32"), tr, trossos[0].crs, trossos[0].nodata


def a_utm31(z, tr, crs, nodata, pas):
    dst = "EPSG:25831"
    w, s, e, n = transform_bounds("EPSG:4326", dst, *CAIXA)
    ncol, nfil = int((e - w) // pas), int((n - s) // pas)
    out = np.full((nfil, ncol), np.nan, dtype="float32")
    dtr = from_origin(w, n, pas, pas)
    reproject(z, out, src_transform=tr, src_crs=crs, src_nodata=nodata, dst_transform=dtr, dst_crs=dst,
              dst_nodata=np.nan, resampling=Resampling.average)
    return out, dtr, dst


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--pas", type=float, default=10)
    ap.add_argument("--llindar-km2", type=float, default=0.05)
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    d = Path("data") / cfg["conca"]
    TMP.mkdir(parents=True, exist_ok=True)

    z5, tr5, crs5, nd5 = baixa_mdt()
    z, tr, crs = a_utm31(z5, tr5, crs5, nd5, a.pas)
    z[z < -50] = np.nan
    mar = np.isnan(z) | (z <= 0.2)
    perfil = dict(driver="GTiff", height=z.shape[0], width=z.shape[1], count=1, dtype="float32", crs=crs, transform=tr,
                  nodata=np.nan, compress="deflate")
    with rasterio.open(TMP / "mdt.tif", "w", **perfil) as f:
        f.write(np.where(mar, np.nan, z), 1)
    print(f"Reixa de {a.pas:.0f} m: {z.shape[1]}×{z.shape[0]}")

    from pysheds.grid import Grid
    g = Grid.from_raster(str(TMP / "mdt.tif"))
    dem = g.read_raster(str(TMP / "mdt.tif"))
    omplert = g.fill_depressions(g.fill_pits(dem))
    plana = g.resolve_flats(omplert)
    fd = g.flowdir(plana)
    acc = g.accumulation(fd)
    cel_km2 = a.pas ** 2 / 1e6
    xarxa = acc > a.llindar_km2 / cel_km2
    ordre = g.stream_order(fd, xarxa)
    rius = g.extract_river_network(fd, xarxa)
    t = Transformer.from_crs(crs, 4326, always_xy=True)
    trams = []
    for f in rius["features"]:
        coords = f["geometry"]["coordinates"]
        if len(coords) < 2:
            continue
        cx, cy = coords[-1]
        fila, col = int((tr.f - cy) // a.pas), int((cx - tr.c) // a.pas)     # cel·la de sortida del tram
        fila = min(max(fila, 0), acc.shape[0] - 1); col = min(max(col, 0), acc.shape[1] - 1)
        ll = [t.transform(x, y) for x, y in coords[:: max(1, len(coords) // 60)] + [coords[-1]]]
        trams.append({"type": "Feature", "properties": {"area_km2": round(float(acc[fila, col]) * cel_km2, 3),
                                                       "ordre": int(ordre[fila, col])},
                      "geometry": {"type": "LineString", "coordinates": [[round(x, 6), round(y, 6)] for x, y in ll]}})
    (d / "rieretes.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": trams}, separators=(",", ":")))
    print(f"rieretes: {len(trams)} trams (conca > {a.llindar_km2} km²)")

    # zones baixes: HAND sobre la xarxa i depressions sense sortida
    hand = np.asarray(g.compute_hand(fd, plana, xarxa), dtype="float32")
    prof = np.asarray(omplert, dtype="float32") - np.asarray(dem, dtype="float32")
    with rasterio.open(TMP / "hand.tif", "w", **perfil) as f:
        f.write(hand, 1)
    classe = np.zeros(z.shape, dtype=np.uint8)
    classe[(hand < 2)] = 1
    classe[(hand < 1)] = 2
    classe[(hand < 0.5)] = 3
    classe[prof > 0.3] = 4                      # depressió: l'aigua s'hi queda
    classe[mar] = 0
    # PNG en lat/lon per superposar-lo amb Leaflet
    w, s, e, n = CAIXA
    pas_g = a.pas / 111320
    ll_tr = from_origin(w, n, pas_g, pas_g)
    ncol, nfil = int((e - w) / pas_g), int((n - s) / pas_g)
    cl = np.zeros((nfil, ncol), dtype=np.uint8)
    reproject(classe, cl, src_transform=tr, src_crs=crs, dst_transform=ll_tr, dst_crs="EPSG:4326", resampling=Resampling.nearest)
    paleta = {1: (198, 219, 239, 110), 2: (107, 174, 214, 150), 3: (33, 113, 181, 180), 4: (128, 0, 128, 200)}
    rgba = np.zeros((nfil, ncol, 4), dtype=np.uint8)
    for k, c in paleta.items():
        rgba[cl == k] = c
    from PIL import Image
    Image.fromarray(rgba, "RGBA").save(d / "zones_baixes.png", optimize=True)
    (d / "zones_baixes.json").write_text(json.dumps({"limits": [[s, w], [n, e]], "pas_m": a.pas,
        "classes": {"1": "a menys de 2 m del llit", "2": "a menys d'1 m del llit", "3": "a menys de 0,5 m del llit",
                    "4": "depressió sense sortida (> 0,3 m)"}}, ensure_ascii=False))
    for k, v in {1: "<2 m", 2: "<1 m", 3: "<0,5 m", 4: "depressió"}.items():
        print(f"  {v}: {(classe == k).sum() * cel_km2:.2f} km²")
    return 0


if __name__ == "__main__":
    sys.exit(main())
