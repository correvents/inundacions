#!/usr/bin/env python3
"""Pluja del radar OPERA (EUMETNET Open Radar Data, CC BY 4.0) sobre cada subconca, cada 5 minuts.

Per a cada compost de reflectivitat (DBZH, 5 min, ~1 km) que encara no sigui a Supabase:
  - mostreja la reflectivitat en una reixa de punts (~250 m) dins de cada subconca,
  - la passa a mm/h amb la relació Z-R ajustada amb les estacions (config/<conca>_zr.json, calibra_zr.py) o, als
    episodis passats (--des), amb Marshall-Palmer (Z = 200 R^1.6); màx. 125 mm/h,
  - desa a inundacions.pluja_conca: mm caiguts en els 5 min (mitjana de la subconca), intensitat màxima i cobertura.
L'acumulat (1/3/6/24 h) el calcula la vista inundacions.pluja_acumulada.
El bucket de temps real només guarda les últimes 24 h; per a episodis passats es fa servir l'arxiu
(mateixos noms de fitxer, DBZH cada 5 min des de juliol de 2024).

Ús:  python3 scripts/radar_conca.py --config config/ridaura.json [--hores 1]
     python3 scripts/radar_conca.py --config config/ridaura.json --des 2026-10-03T00:00 --fins 2026-10-05T00:00"""
import argparse, json, re, sys, tempfile, urllib.parse
from datetime import datetime, timedelta, timezone
from pathlib import Path
import numpy as np, h5py, requests
from pyproj import CRS, Transformer
from shapely.geometry import shape, Point
from shapely.prepared import prep
from sb import insereix, llegeix

S3_24H = "https://s3.waw3-1.cloudferro.com/openradar-24h"
S3_ARXIU = "https://s3.waw3-1.cloudferro.com/openradar-archive"
PAS = 0.0025   # graus (~250 m): prou per tenir diversos punts per píxel de radar


def claus(s3, ini, fi):
    """{instant (UTC): clau S3} dels compostos DBZH entre ini i fi."""
    out, h = {}, ini.replace(minute=0)
    while h <= fi:
        xml = requests.get(s3, params={"list-type": 2, "prefix": h.strftime("%Y/%m/%d/OPERA/COMP/OPERA@%Y%m%dT%H")}, timeout=30).text
        for k in re.findall(r"<Key>([^<]+)</Key>", xml):
            m = re.search(r"OPERA@(\d{8}T\d{4})@0@DBZH(_QIND)?\.h5$", k)      # els anys antics porten _QIND
            if m:
                t = datetime.strptime(m.group(1), "%Y%m%dT%H%M").replace(tzinfo=timezone.utc)
                if ini < t <= fi:
                    out[t] = k
        h += timedelta(hours=1)
    return out


def durades(ts, defecte=5):
    """Minuts que representa cada compost: l'interval fins al següent (5 min ara, 15 min en anys antics)."""
    ts = sorted(ts)
    d = {t: min(max(int((b - t).total_seconds() // 60), 1), 15) for t, b in zip(ts, ts[1:])}
    if ts:
        d[ts[-1]] = d[ts[-2]] if len(ts) > 1 else defecte
    return d


def punts(geom):
    w, s, e, n = geom.bounds
    g = prep(geom)
    lon, lat = np.meshgrid(np.arange(w, e, PAS) + PAS / 2, np.arange(s, n, PAS) + PAS / 2)
    dins = np.array([g.contains(Point(x, y)) for x, y in zip(lon.ravel(), lat.ravel())])
    return lon.ravel()[dins], lat.ravel()[dins]


def mostreja(h5, lon, lat):
    """dBZ a cada punt; NaN = sense cobertura; -inf = cobertura sense eco."""
    with h5py.File(h5, "r") as h:
        # ODIM: els atributs 'what' poden ser a datasetN/data1 (ara) o a datasetN (arxiu antic, amb la capa QIND a dataset2)
        grup, wt = None, None
        for ds_nom in [k for k in h if k.startswith("dataset")]:
            for cami in (f"{ds_nom}/data1/what", f"{ds_nom}/what"):
                if cami in h:
                    q = h[cami].attrs.get("quantity")
                    if q is not None and (q.decode() if isinstance(q, bytes) else str(q)) == "DBZH":
                        grup, wt = f"{ds_nom}/data1", h[cami].attrs
                        break
            if grup:
                break
        if grup is None:
            grup, wt = "dataset1/data1", h["dataset1/data1/what"].attrs
        ds = h[grup + "/data"]
        wh = h["where"].attrs
        gain, off, nodata, undet = (float(wt[k]) for k in ("gain", "offset", "nodata", "undetect"))
        pd = wh["projdef"].decode() if isinstance(wh["projdef"], bytes) else str(wh["projdef"])
        tr = Transformer.from_crs("EPSG:4326", CRS.from_proj4(pd), always_xy=True)
        xul, yul = tr.transform(float(wh["UL_lon"]), float(wh["UL_lat"]))
        X, Y = tr.transform(lon, lat)
        col = np.floor((X - xul) / float(wh["xscale"])).astype(int)     # UL és la cantonada del primer píxel
        row = np.floor((yul - Y) / float(wh["yscale"])).astype(int)
        r0, c0 = row.min(), col.min()
        tros = ds[r0:row.max() + 1, c0:col.max() + 1].astype(np.float64)
    v = tros[row - r0, col - c0]
    z = v * gain + off
    z[v == undet] = -np.inf
    z[v == nodata] = np.nan
    return z


def relacio_zr(conca):
    """Relació Z-R ajustada amb les estacions (calibra_zr.py) si n'hi ha; si no, Marshall-Palmer."""
    f = Path("config") / f"{conca}_zr.json"
    return json.loads(f.read_text()) if f.exists() else None


def mmh(z, zr=None):
    """mm/h a partir de la reflectivitat (dBZ). Sense zr, Marshall-Palmer (Z = 200 R^1.6); amb zr, una relació per a la
    pluja estratiforme i una altra per a la convectiva, separades per zr["llindar_dbz"]."""
    zz = np.where(np.isfinite(z), z, 0)
    lin = 10 ** (zz / 10)
    if zr is None:
        r = (lin / 200) ** (1 / 1.6)
    else:
        e, c = zr["estratiforme"], zr["convectiva"]
        r = np.where(zz >= zr["llindar_dbz"], (lin / c["a"]) ** (1 / c["b"]), (lin / e["a"]) ** (1 / e["b"]))
    r = np.where(np.isfinite(z), r, 0.0)
    r[np.isnan(z)] = np.nan
    return np.minimum(r, 125.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--hores", type=float, default=1)
    ap.add_argument("--des", help="inici UTC (YYYY-MM-DDTHH:MM): llegeix de l'arxiu en lloc de les últimes hores")
    ap.add_argument("--fins", help="final UTC (per defecte, 24 h després de --des)")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    ZR = None if a.des else relacio_zr(cfg["conca"])     # els episodis passats es desen amb Marshall-Palmer (calibra_zr.py en parteix)
    fc = json.loads(Path(cfg["geometria"]).read_text())
    subs = {}
    for f in fc["features"]:
        if f["geometry"]["type"].endswith("Polygon"):
            subs[f["properties"]["subconca"]] = punts(shape(f["geometry"]))
    if a.des:
        s3 = S3_ARXIU
        ini = datetime.fromisoformat(a.des).replace(tzinfo=timezone.utc)
        fi = datetime.fromisoformat(a.fins).replace(tzinfo=timezone.utc) if a.fins else ini + timedelta(hours=24)
    else:
        s3 = S3_24H
        fi = datetime.now(timezone.utc); ini = fi - timedelta(hours=a.hores)
    disponibles = claus(s3, ini, fi)
    fets = {r["temps"][:16] for r in llegeix("pluja_conca", {"select": "temps", "conca": f"eq.{cfg['conca']}", "font": "eq.radar",
                                                               "temps": f"gt.{ini.isoformat()}",
                                                               "and": f"(temps.lte.{fi.isoformat()})"})}
    pendents = {t: k for t, k in disponibles.items() if t.strftime("%Y-%m-%dT%H:%M") not in fets}
    minuts = durades(disponibles)
    print(f"OPERA: {len(disponibles)} compostos a la finestra, {len(pendents)} per processar")
    files, desades = [], 0
    with tempfile.TemporaryDirectory() as tmp:
        for t, k in sorted(pendents.items()):
            f = Path(tmp) / "c.h5"
            try:
                r = requests.get(s3 + "/" + urllib.parse.quote(k, safe="/@"), timeout=120)
                r.raise_for_status(); f.write_bytes(r.content)
                for nom, (lon, lat) in subs.items():
                    i = mmh(mostreja(f, lon, lat), ZR)
                    ok = ~np.isnan(i)
                    files.append({"conca": cfg["conca"], "subconca": nom, "font": "radar", "temps": t.isoformat(), "minuts": minuts[t],
                                  "mm": round(float(i[ok].mean()) * minuts[t] / 60, 3) if ok.any() else None,
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
