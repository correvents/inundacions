#!/usr/bin/env python3
"""Pluja del radar OPERA (EUMETNET Open Radar Data, CC BY 4.0) sobre cada subconca, cada 5 minuts.

Per a cada compost de reflectivitat (DBZH, 5 min, ~1 km) que encara no sigui a Supabase:
  - mostreja la reflectivitat en una reixa de punts (~250 m) dins de cada subconca,
  - la passa a mm/h amb Marshall-Palmer (Z = 200 R^1.6, màx. 125 mm/h),
  - desa a inundacions.pluja_conca: mm caiguts en els 5 min (mitjana de la subconca), intensitat màxima i cobertura.
L'acumulat (1/3/6/24 h) el calcula la vista inundacions.pluja_acumulada.
El bucket públic només guarda les últimes 24 h.

Ús:  python3 scripts/radar_conca.py --config config/ridaura.json [--hores 1]"""
import argparse, json, re, sys, tempfile, urllib.parse
from datetime import datetime, timedelta, timezone
from pathlib import Path
import numpy as np, h5py, requests
from pyproj import CRS, Transformer
from shapely.geometry import shape, Point
from shapely.prepared import prep
from sb import insereix, llegeix

S3 = "https://s3.waw3-1.cloudferro.com/openradar-24h"
PAS = 0.0025   # graus (~250 m): prou per tenir diversos punts per píxel de radar


def claus(ini, fi):
    """{instant (UTC): clau S3} dels compostos DBZH entre ini i fi."""
    out, h = {}, ini.replace(minute=0)
    while h <= fi:
        xml = requests.get(S3, params={"list-type": 2, "prefix": h.strftime("%Y/%m/%d/OPERA/COMP/OPERA@%Y%m%dT%H")}, timeout=30).text
        for k in re.findall(r"<Key>([^<]+)</Key>", xml):
            m = re.search(r"OPERA@(\d{8}T\d{4})@0@DBZH\.h5$", k)
            if m:
                t = datetime.strptime(m.group(1), "%Y%m%dT%H%M").replace(tzinfo=timezone.utc)
                if ini < t <= fi:
                    out[t] = k
        h += timedelta(hours=1)
    return out


def punts(geom):
    w, s, e, n = geom.bounds
    g = prep(geom)
    lon, lat = np.meshgrid(np.arange(w, e, PAS) + PAS / 2, np.arange(s, n, PAS) + PAS / 2)
    dins = np.array([g.contains(Point(x, y)) for x, y in zip(lon.ravel(), lat.ravel())])
    return lon.ravel()[dins], lat.ravel()[dins]


def mostreja(h5, lon, lat):
    """dBZ a cada punt; NaN = sense cobertura; -inf = cobertura sense eco."""
    with h5py.File(h5, "r") as h:
        ds = h["dataset1/data1/data"]
        wt, wh = h["dataset1/data1/what"].attrs, h["where"].attrs
        gain, off, nodata, undet = (float(wt[k]) for k in ("gain", "offset", "nodata", "undetect"))
        pd = wh["projdef"].decode() if isinstance(wh["projdef"], bytes) else str(wh["projdef"])
        tr = Transformer.from_crs("EPSG:4326", CRS.from_proj4(pd), always_xy=True)
        xul, yul = tr.transform(float(wh["UL_lon"]), float(wh["UL_lat"]))
        X, Y = tr.transform(lon, lat)
        col = np.rint((X - xul) / float(wh["xscale"])).astype(int)
        row = np.rint((yul - Y) / float(wh["yscale"])).astype(int)
        r0, c0 = row.min(), col.min()
        tros = ds[r0:row.max() + 1, c0:col.max() + 1].astype(np.float64)
    v = tros[row - r0, col - c0]
    z = v * gain + off
    z[v == undet] = -np.inf
    z[v == nodata] = np.nan
    return z


def mmh(z):
    r = np.where(np.isfinite(z), (10 ** (np.where(np.isfinite(z), z, 0) / 10) / 200) ** (1 / 1.6), 0.0)
    r[np.isnan(z)] = np.nan
    return np.minimum(r, 125.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--hores", type=float, default=1)
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    fc = json.loads(Path(cfg["geometria"]).read_text())
    subs = {}
    for f in fc["features"]:
        if f["geometry"]["type"].endswith("Polygon"):
            subs[f["properties"]["subconca"]] = punts(shape(f["geometry"]))
    fi = datetime.now(timezone.utc); ini = fi - timedelta(hours=a.hores)
    disponibles = claus(ini, fi)
    fets = {r["temps"][:16] for r in llegeix("pluja_conca", {"select": "temps", "conca": f"eq.{cfg['conca']}", "font": "eq.radar",
                                                               "temps": f"gt.{ini.isoformat()}"})}
    pendents = {t: k for t, k in disponibles.items() if t.strftime("%Y-%m-%dT%H:%M") not in fets}
    print(f"OPERA: {len(disponibles)} compostos a la finestra, {len(pendents)} per processar")
    files, desades = [], 0
    with tempfile.TemporaryDirectory() as tmp:
        for t, k in sorted(pendents.items()):
            f = Path(tmp) / "c.h5"
            try:
                r = requests.get(S3 + "/" + urllib.parse.quote(k, safe="/@"), timeout=120)
                r.raise_for_status(); f.write_bytes(r.content)
                for nom, (lon, lat) in subs.items():
                    i = mmh(mostreja(f, lon, lat))
                    ok = ~np.isnan(i)
                    files.append({"conca": cfg["conca"], "subconca": nom, "font": "radar", "temps": t.isoformat(), "minuts": 5,
                                  "mm": round(float(i[ok].mean()) * 5 / 60, 3) if ok.any() else None,
                                  "mmh_max": round(float(i[ok].max()), 2) if ok.any() else None,
                                  "cobertura": round(float(ok.mean()), 3)})
            except Exception as e:
                print(f"OPERA {t:%Y-%m-%d %H:%M}: {type(e).__name__} {e}")
            if len(files) >= 24:                 # desa a mesura que avança (la finestra de 24 h és llarga)
                desades += insereix("pluja_conca", files); files = []
    if files:
        desades += insereix("pluja_conca", files)
    print(f"OPERA: {desades} files desades")
    return 0


if __name__ == "__main__":
    sys.exit(main())
