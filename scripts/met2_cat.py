#!/usr/bin/env python3
"""MET de 2 m de l'ICGC a les zones urbanes de Catalunya, conca per conca (opció de relleu fi per a la simulació 2D).

Fonts (--font), models d'elevacions del terreny de Catalunya de l'ICGC, cadascun un sol GeoTIFF optimitzat per al núvol
(COG) a datacloud.icgc.cat, EPSG:25831, que es llegeix per finestres (/vsicurl):
  2008: MET 2×2 m, LiDAR 2008-2011 (30 GB), tal qual                       → met2_cNNN.tif   (Release met2-cat)
  2021: MET 1×1 m, LiDAR 2021-2023 (92 GB), mitjana de cada 2×2 cel·les → 2 m → met21_cNNN.tif  (Release met21-cat)
datacloud no respon des de GitHub Actions: aquest script s'executa des d'un ordinador normal i els fitxers es pugen
a la Release.

Per a cada conca del tall de 100 km² (data/catalunya/conques/conques_100km2.geojson):
  1. Zones urbanes: les cel·les construïdes del mapa de cobertes de l'ICGC (Release cobertes-cat, tipus 1 a
     config/cobertes_parametres.json), en blocs de 64 m, eixamplades --marge m i retallades a la conca + 500 m.
  2. Llegeix el MET 2 m només als blocs de 512 m que en toquen i deixa sense dades la resta.
  3. Desa el fitxer de la conca: reixa de 2 m alineada amb la de l'ICGC (EPSG:25831), alçada en cm (int32), sense
     dades -3276800 (com mdt5_cNNN.tif). Les conques sense zones urbanes no tenen fitxer.

Ús:  python3 scripts/met2_cat.py conca --id 124 [--font 2008|2021] [--sortida DIR]
     python3 scripts/met2_cat.py totes [--font 2008|2021] [--conques 1,2,3] [--fils 4] [--sortida DIR]   (salta les fetes)"""
import argparse, json, os, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("GDAL_HTTP_MULTIRANGE", "YES")
os.environ.setdefault("GDAL_HTTP_MAX_RETRY", "5")
os.environ.setdefault("GDAL_HTTP_RETRY_DELAY", "5")
import rasterio
from rasterio import features
from rasterio.transform import Affine
from rasterio.warp import reproject, Resampling
from rasterio.windows import Window
from pyproj import Transformer
from shapely.geometry import shape, mapping
from shapely.ops import transform as stransform

sys.path.insert(0, str(Path(__file__).parent))
import cobertes as C

DC = "/vsicurl/https://datacloud.icgc.cat/datacloud/model-elevacions-terreny/tif_unzip/"
FONTS = {   # camí, factor fins a 2 m, prefix del fitxer, descripció
    "2008": (DC + "model-elevacions-terreny-lidar-catalunya-2m-2008-2011.tif", 1, "met2",
             "ICGC, Model d'elevacions del terreny de Catalunya 2x2 m, LiDAR 2008-2011 (CC BY 4.0)"),
    "2021": (DC + "model-elevacions-terreny-lidar-catalunya-1m-2021-2023.tif", 2, "met21",
             "ICGC, Model d'elevacions del terreny de Catalunya 1x1 m, LiDAR 2021-2023 (CC BY 4.0), mitjana a 2 m"),
}
CONQUES = Path("data/catalunya/conques/conques_100km2.geojson")
ND = -32768 * 100
BLOC = 32            # píxels de 2 m: blocs de 64 m per a la màscara urbana
TROS = 256           # píxels de 2 m: lectures de 512 m (els blocs interns del COG)
MARGE_CONCA = 500


def conques():
    return {f["properties"]["id"]: f for f in json.loads(CONQUES.read_text())["features"]}


def eixampla(m, n):
    for _ in range(n):
        p = np.pad(m, 1)
        m = (p[1:-1, 1:-1] | p[:-2, 1:-1] | p[2:, 1:-1] | p[1:-1, :-2] | p[1:-1, 2:]
             | p[:-2, :-2] | p[:-2, 2:] | p[2:, :-2] | p[2:, 2:])
    return m


def llegeix(src, fac, c, f, w, h):
    """Finestra de la reixa de 2 m (columna c, fila f, w×h); amb fac 2, mitjana de 2×2 cel·les d'1 m (sense dades si
    en falta alguna)."""
    v = src.read(1, window=Window(c * fac, f * fac, w * fac, h * fac), boundless=True, fill_value=src.nodata)
    v = v.astype("float32")
    v[(v == src.nodata) | ~np.isfinite(v) | (v < -100)] = np.nan
    if fac > 1:
        v = v.reshape(h, fac, w, fac).mean(axis=(1, 3))
    return v


def conca(i, feat, src, sortida, marge, font="2008"):
    _, agr, prefix, desc = FONTS[font]
    nom = f"c{i:03d}"
    out = sortida / f"{prefix}_{nom}.tif"
    t = Transformer.from_crs(4326, 25831, always_xy=True)
    g = stransform(t.transform, shape(feat["geometry"])).buffer(MARGE_CONCA)
    T = src.transform * Affine.scale(agr)
    W2, H2 = src.width // agr, src.height // agr
    x0, y0, x1, y1 = g.bounds
    c0 = int(np.floor((x0 - T.c) / T.a / TROS)) * TROS
    f0 = int(np.floor((T.f - y1) / -T.e / TROS)) * TROS
    c1 = int(np.ceil((x1 - T.c) / T.a / TROS)) * TROS
    f1 = int(np.ceil((T.f - y0) / -T.e / TROS)) * TROS
    c0, f0 = max(c0, 0), max(f0, 0)
    c1, f1 = min(c1, W2), min(f1, H2)
    nf, nc = f1 - f0, c1 - c0
    tr = T * Affine.translation(c0, f0)

    # màscara urbana en blocs de 64 m (cobertes de 5 m → reixa de 16 m, qualsevol cel·la urbana marca el bloc)
    fac = BLOC // 8
    trc = tr * Affine.scale(8)
    cod = np.zeros((nf // 8, nc // 8), np.uint8)
    with rasterio.open(C.COB_URL.format(nom)) as cs:
        reproject(rasterio.band(cs, 1), cod, dst_transform=trc, dst_crs="EPSG:25831", src_nodata=0, dst_nodata=0,
                  resampling=Resampling.nearest)
    _, _, tipus = C.parametres(cod)
    urb = (tipus == 1) & (cod > 0)
    bf, bc = nf // BLOC, nc // BLOC
    urb = urb[: bf * fac, : bc * fac].reshape(bf, fac, bc, fac).any(axis=(1, 3))
    urb = eixampla(urb, int(np.ceil(marge / (BLOC * 2))))
    dins = features.geometry_mask([mapping(g)], (bf, bc), tr * Affine.scale(BLOC), invert=True, all_touched=True)
    urb &= dins
    if not urb.any():
        print(f"{nom}: sense zones urbanes"); return None

    # retalla als trossos de 512 m que tenen zona urbana
    k = TROS // BLOC
    tf, tc = -(-bf // k), -(-bc // k)
    pad = np.zeros((tf * k, tc * k), bool); pad[:bf, :bc] = urb
    tros = pad.reshape(tf, k, tc, k).any(axis=(1, 3))
    fs, cs_ = np.where(tros.any(axis=1))[0], np.where(tros.any(axis=0))[0]
    tf0, tf1, tc0, tc1 = fs[0], fs[-1] + 1, cs_[0], cs_[-1] + 1
    H, W = min((tf1 - tf0) * TROS, nf - tf0 * TROS), min((tc1 - tc0) * TROS, nc - tc0 * TROS)
    z = np.full((H, W), ND, np.int32)
    llegits = 0
    for a in range(tf0, tf1):
        for b in range(tc0, tc1):
            if not tros[a, b]:
                continue
            fa, ca = a * TROS, b * TROS
            h, w = min(TROS, nf - fa), min(TROS, nc - ca)
            if h <= 0 or w <= 0:
                continue
            for intent in range(6):
                try:
                    v = llegeix(src, agr, c0 + ca, f0 + fa, w, h); break
                except rasterio.errors.RasterioIOError as e:
                    if intent == 5: raise
                    time.sleep(10 * (intent + 1))
            m = np.repeat(np.repeat(pad[a * k:(a + 1) * k, b * k:(b + 1) * k], BLOC, 0), BLOC, 1)[:h, :w]
            ok = m & np.isfinite(v)
            bloc = np.where(ok, np.round(v * 100), ND).astype(np.int32)
            z[fa - tf0 * TROS: fa - tf0 * TROS + h, ca - tc0 * TROS: ca - tc0 * TROS + w] = bloc
            llegits += 1
    if not (z != ND).any():
        print(f"{nom}: el MET 2 m no té dades a les zones urbanes"); return None
    trz = tr * Affine.translation(tc0 * TROS, tf0 * TROS)
    tmp = out.with_suffix(".part")
    with rasterio.open(tmp, "w", driver="GTiff", height=H, width=W, count=1, dtype="int32", crs="EPSG:25831",
                       transform=trz, nodata=ND, compress="deflate", predictor=2, zlevel=9, tiled=True,
                       blockxsize=256, blockysize=256) as f:
        f.write(z, 1)
        f.update_tags(unitats="cm", conca=nom, marge_urba_m=str(marge), font=desc)
    tmp.rename(out)
    km2 = (z != ND).sum() * 4 / 1e6
    print(f"{prefix}_{nom}: {km2:.1f} km² amb dades, {llegits} trossos, {out.stat().st_size / 1e6:.1f} MB")
    return km2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="ordre", required=True)
    for nom in ("conca", "totes"):
        s = sub.add_parser(nom)
        s.add_argument("--sortida", default="/tmp/met2_cat")
        s.add_argument("--marge", type=float, default=250, help="m al voltant de les zones construïdes")
        s.add_argument("--font", choices=sorted(FONTS), default="2008")
    sub.choices["conca"].add_argument("--id", type=int, required=True)
    sub.choices["totes"].add_argument("--conques", default="totes")
    sub.choices["totes"].add_argument("--fils", type=int, default=4)
    a = ap.parse_args()
    sortida = Path(a.sortida); sortida.mkdir(parents=True, exist_ok=True)
    CQ = conques()
    if a.ordre == "conca":
        with rasterio.open(FONTS[a.font][0]) as src:
            conca(a.id, CQ[a.id], src, sortida, a.marge, a.font)
        return
    ids = sorted(CQ) if a.conques == "totes" else [int(x) for x in a.conques.split(",")]
    fet = sortida / ("fetes.txt" if a.font == "2008" else f"fetes_{a.font}.txt")
    ja = set(fet.read_text().split()) if fet.exists() else set()
    ids = [i for i in ids if str(i) not in ja]
    print(f"{len(ids)} conques per fer")

    def una(i):
        try:
            with rasterio.open(FONTS[a.font][0]) as src:   # un lector per fil (GDAL no és segur entre fils)
                conca(i, CQ[i], src, sortida, a.marge, a.font)
            with open(fet, "a") as f:
                f.write(f"{i}\n")
        except Exception as e:
            print(f"c{i:03d}: ERROR {type(e).__name__} {e}")
    with ThreadPoolExecutor(a.fils) as ex:
        list(ex.map(una, ids))


if __name__ == "__main__":
    main()
