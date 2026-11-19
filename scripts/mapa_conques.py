#!/usr/bin/env python3
"""Mapa de conques de Catalunya fet amb el nostre propi MDT (IGN, 25 m → reixa de --pas m, ETRS89 UTM 31N).

Les conques es parteixen com un arbre: a cada confluència on arriben dos o més afluents de més de
--llindar km², cada afluent passa a ser una conca pròpia; el tram de sota és una altra conca.
Les rieres que desemboquen al mar amb més de --costa km² també són conques pròpies.
Serveix per decidir com partir el MDT (Releases) i com vigilar tot Catalunya amb el radar.

Escriu, per a cada llindar:
  data/catalunya/conques/conques_LLkm2.geojson   (id, km², punt de sortida)
  data/catalunya/conques/conques_LLkm2.png
  data/catalunya/conques/rius.geojson             (rius > 50 km²)
  data/catalunya/conques/resum.json

Ús:  python3 scripts/mapa_conques.py [--pas 100] [--llindar 25 100] [--costa 2]"""
import argparse, json, time
from pathlib import Path
import numpy as np, requests
np.in1d = np.isin
import rasterio, geopandas as gpd
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from rasterio import features
from pyproj import Transformer
from shapely.geometry import shape
from numba import njit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WCS = "https://servicios.idee.es/wcs-inspire/mdt"
COBERTURES = ["Elevacion25830_25", "Elevacion25830_5"]
NE = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_admin_1_states_provinces.geojson"
CAIXA = (0.10, 40.45, 3.40, 42.90)
SORT = Path("data/catalunya/conques")
TMP = Path("/tmp/conques")
POBLES = {"Platja d'Aro": (3.067, 41.817), "Granollers": (2.288, 41.608), "Cabrera de Mar": (2.396, 41.527)}
# direccions de pysheds (N, NE, E, SE, S, SW, W, NW)
DIRS = np.array([64, 128, 1, 2, 4, 8, 16, 32]); DR = np.array([-1, -1, 0, 1, 1, 1, 0, -1]); DC = np.array([0, 1, 1, 1, 0, -1, -1, -1])


def peca(cob, xa, ya, xb, yb):
    p = [("SERVICE", "WCS"), ("VERSION", "2.0.1"), ("REQUEST", "GetCoverage"), ("COVERAGEID", cob),
         ("SUBSET", f"x({xa:.0f},{xb:.0f})"), ("SUBSET", f"y({ya:.0f},{yb:.0f})"), ("FORMAT", "image/tiff")]
    r = None
    for intent in range(6):
        try:
            r = requests.get(WCS, params=p, timeout=300)
            if r.content[:2] in (b"II", b"MM"):
                return r.content
        except requests.RequestException:
            pass
        time.sleep(10 * (intent + 1))
    raise RuntimeError(f"WCS {cob} {xa:.0f},{ya:.0f}: {getattr(r, 'status_code', '?')} {getattr(r, 'text', '')[:200]}")


def mdt(pas):
    t30 = Transformer.from_crs(4326, 25830, always_xy=True); t31 = Transformer.from_crs(4326, 25831, always_xy=True)
    xs, ys = zip(*[t31.transform(x, y) for x in CAIXA[::2] for y in CAIXA[1::2]])
    w, s, e, n = min(xs), min(ys), max(xs), max(ys)
    dtr = from_origin(w, n, pas, pas); out = np.full((int((n - s) // pas), int((e - w) // pas)), np.nan, "float32")
    xs, ys = zip(*[t30.transform(x, y) for x in CAIXA[::2] for y in CAIXA[1::2]])
    x0, y0, x1, y1 = min(xs) - 2000, min(ys) - 2000, max(xs) + 2000, max(ys) + 2000
    cob, mida = COBERTURES[0], 50000
    TMP.mkdir(parents=True, exist_ok=True)
    k = 0
    for xa in np.arange(x0, x1, mida):
        for ya in np.arange(y0, y1, mida):
            f = TMP / f"t{k}.tif"; k += 1
            f.write_bytes(peca(cob, xa, ya, min(xa + mida, x1), min(ya + mida, y1)))
            with rasterio.open(f) as src:
                z = src.read(1).astype("float32"); nd = src.nodata
                if nd is not None: z[z == nd] = np.nan
                z[z < -50] = np.nan
                tmp = np.full_like(out, np.nan)
                reproject(z, tmp, src_transform=src.transform, src_crs=src.crs, src_nodata=np.nan, dst_transform=dtr,
                          dst_crs="EPSG:25831", dst_nodata=np.nan, resampling=Resampling.average)
            ok = np.isfinite(tmp) & np.isnan(out); out[ok] = tmp[ok]
            print(f"tros {k}: {np.isfinite(out).mean():.0%} ple", flush=True)
            f.unlink()
    out[out <= 0.2] = np.nan                                 # mar
    return out, dtr


@njit(cache=True)
def etiqueta(acc, avall, ordre, sortida):
    lab = np.full(acc.size, -1, np.int64)
    for c in ordre:
        if sortida[c]:
            lab[c] = c
        elif avall[c] >= 0:
            lab[c] = lab[avall[c]]
    return lab


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pas", type=float, default=100)
    ap.add_argument("--llindar", type=float, nargs="+", default=[25, 100])
    ap.add_argument("--costa", type=float, default=2)
    a = ap.parse_args()
    SORT.mkdir(parents=True, exist_ok=True)
    z, tr = mdt(a.pas)
    perfil = dict(driver="GTiff", height=z.shape[0], width=z.shape[1], count=1, dtype="float32", crs="EPSG:25831",
                  transform=tr, nodata=np.nan, compress="deflate")
    with rasterio.open(TMP / "mdt.tif", "w", **perfil) as f:
        f.write(z, 1)
    print(f"MDT {a.pas:.0f} m: {z.shape[1]}×{z.shape[0]}", flush=True)

    from pysheds.grid import Grid
    g = Grid.from_raster(str(TMP / "mdt.tif")); dem = g.read_raster(str(TMP / "mdt.tif"))
    fd_r = g.flowdir(g.resolve_flats(g.fill_depressions(g.fill_pits(dem))))
    fd = np.asarray(fd_r); acc = np.asarray(g.accumulation(fd_r), dtype="float64")
    print("drenatge fet", flush=True)
    nf, nc = z.shape; terra = np.isfinite(z)
    cel = a.pas ** 2 / 1e6
    # cel·la de sota
    avall = np.full(z.size, -1, np.int64)
    rr, cc = np.indices(z.shape)
    for d, dr, dc in zip(DIRS, DR, DC):
        m = (fd == d) & terra
        r2, c2 = rr[m] + dr, cc[m] + dc
        dins = (r2 >= 0) & (r2 < nf) & (c2 >= 0) & (c2 < nc)
        idx = np.flatnonzero(m.ravel())[dins]
        r2, c2 = r2[dins], c2[dins]
        bo = terra[r2, c2]
        avall[idx[bo]] = (r2 * nc + c2)[bo]
    accf = acc.ravel(); terraf = terra.ravel()
    ordre = np.argsort(-accf, kind="stable"); ordre = ordre[terraf[ordre]]

    prov = gpd.read_file(NE)
    cat = prov[prov["iso_3166_2"].isin(["ES-B", "ES-GI", "ES-L", "ES-T"])].dissolve().to_crs(25831)
    cat_geom = cat.geometry.iloc[0]
    t_ll = Transformer.from_crs(25831, 4326, always_xy=True)

    # rius > 50 km² per al mapa
    rius = g.extract_river_network(fd_r, acc > 50 / cel)
    rg = gpd.GeoDataFrame.from_features(rius["features"], crs=25831)
    rg = rg[rg.intersects(cat_geom)]
    rg.to_crs(4326).to_file(SORT / "rius.geojson", driver="GeoJSON", COORDINATE_PRECISION=5)

    resum = {"pas_m": a.pas, "costa_km2": a.costa, "font": "MDT IGN 25 m (servicios.idee.es)"}
    for llindar in a.llindar:
        T = llindar / cel
        gran = (accf >= T) & terraf
        # afluents grans que conflueixen: cada un és sortida de conca
        n_gran = np.bincount(avall[gran & (avall >= 0)], minlength=z.size)
        sortida = gran & (avall >= 0) & (n_gran[np.maximum(avall, 0)] >= 2)
        sortida |= terraf & (avall < 0) & (accf >= a.costa / cel)          # desembocadures i vores
        lab = etiqueta(accf, avall, ordre, sortida).reshape(z.shape)
        dins_cat = features.rasterize([cat_geom], out_shape=z.shape, transform=tr, fill=0, default_value=1).astype(bool)
        lab[~dins_cat] = -1
        ids, inv = np.unique(lab, return_inverse=True); lab2 = inv.reshape(z.shape).astype("int32")
        fora = ids == -1
        polis = [(shape(gm), int(ids[v])) for gm, v in features.shapes(lab2, mask=~fora[lab2], transform=tr)]
        gdf = gpd.GeoDataFrame({"sortida": [s for _, s in polis]}, geometry=[p for p, _ in polis], crs=25831).dissolve("sortida").reset_index()
        gdf["km2"] = (gdf.area / 1e6).round(1)
        gdf = gdf[gdf.km2 >= 0.5].copy()
        r_o, c_o = np.divmod(gdf.sortida.values, nc)
        xo, yo = rasterio.transform.xy(tr, r_o, c_o)
        lon, lat = t_ll.transform(xo, yo)
        gdf["sortida_lon"] = np.round(lon, 5); gdf["sortida_lat"] = np.round(lat, 5)
        gdf["km2_amunt"] = np.round(accf[gdf.sortida.values] * cel, 1)
        gdf["id"] = np.arange(1, len(gdf) + 1)
        gdf = gdf[["id", "km2", "km2_amunt", "sortida_lon", "sortida_lat", "geometry"]]
        g2 = gdf.copy(); g2["geometry"] = g2.geometry.simplify(a.pas / 2)
        g2.to_crs(4326).to_file(SORT / f"conques_{llindar:g}km2.geojson", driver="GeoJSON", COORDINATE_PRECISION=5)
        petites = gdf.km2 < 5
        resum[f"{llindar:g}km2"] = {"conques": int(len(gdf)), "conques_>=5km2": int((~petites).sum()),
                                    "km2_mediana_>=5": float(gdf.km2[~petites].median()), "km2_max": float(gdf.km2.max())}
        fig, ax = plt.subplots(figsize=(11, 11), dpi=120)
        gdf[~petites].plot(ax=ax, column="id", cmap="tab20", alpha=0.55, edgecolor="#444", linewidth=0.35)
        gdf[petites].plot(ax=ax, color="#ddd", edgecolor="#bbb", linewidth=0.2)
        rg.plot(ax=ax, color="#1f5fa8", linewidth=0.5)
        cat.boundary.plot(ax=ax, color="black", linewidth=1.1)
        for nom, (lo, la) in POBLES.items():
            x, y = Transformer.from_crs(4326, 25831, always_xy=True).transform(lo, la)
            ax.plot(x, y, "o", color="red", ms=6, mec="white")
            ax.annotate(nom, (x, y), xytext=(5, 5), textcoords="offset points", fontsize=9, weight="bold")
        r = resum[f"{llindar:g}km2"]
        ax.set_title(f"Conques de Catalunya (tall a {llindar:g} km²): {r['conques_>=5km2']} conques ≥5 km², "
                     f"mediana {r['km2_mediana_>=5']:.0f} km²  ·  gris: franja costanera < 5 km²", fontsize=10)
        ax.set_axis_off(); fig.tight_layout(); fig.savefig(SORT / f"conques_{llindar:g}km2.png"); plt.close(fig)
        print(llindar, r, flush=True)
    (SORT / "resum.json").write_text(json.dumps(resum, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
