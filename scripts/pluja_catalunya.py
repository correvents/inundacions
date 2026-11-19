#!/usr/bin/env python3
"""Pluja horària de tot Catalunya a partir del radar OPERA (EUMETNET Open Radar Data, CC BY 4.0), per a la pestanya catalunya.html.

Per a cada hora del mes llegeix el compost d'acumulat d'1 h (ACRR, Marshall-Palmer, l'hora que acaba a HH:00) i el
mostreja en una reixa regular lat/lon de 0,02° (~2 km) sobre Catalunya. Sense correcció amb estacions (és la pluja del radar tal qual).
Escriu data/pluja_cat/AAAA-MM.bin.gz (uint8, [hores][files de N a S][columnes d'O a E]; mm = (v/20)^2, 255 = sense dada)
i actualitza data/pluja_cat/index.json.

Ús:  python3 scripts/pluja_catalunya.py --mes 2026-10"""
import argparse, gzip, json, re, tempfile, calendar
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from pathlib import Path
import numpy as np, h5py, requests
from pyproj import CRS, Transformer

S3 = "https://s3.waw3-1.cloudferro.com/openradar-archive"
OEST, EST, SUD, NORD, PAS = 0.15, 3.35, 40.5, 42.9, 0.02
LON = np.arange(OEST, EST, PAS) + PAS / 2
LAT = np.arange(NORD, SUD, -PAS) - PAS / 2
LO, LA = np.meshgrid(LON, LAT)
_proj = {}


def claus_dia(dia):
    xml = requests.get(S3, params={"list-type": 2, "prefix": dia.strftime("%Y/%m/%d/OPERA/COMP/"), "max-keys": 1000}, timeout=60).text
    out = {}
    for k in re.findall(r"<Key>([^<]+)</Key>", xml):
        m = re.search(r"OPERA@(\d{8}T\d{2})00@0@ACRR\.h5$", k)
        if m:
            out[datetime.strptime(m.group(1), "%Y%m%dT%H").replace(tzinfo=timezone.utc)] = k
    return out


def llegeix(clau):
    for i in range(3):
        try:
            r = requests.get(f"{S3}/{clau}", timeout=120); r.raise_for_status(); break
        except Exception as e:
            if i == 2: print("falla", clau, e); return None
    with tempfile.NamedTemporaryFile(suffix=".h5") as f:
        f.write(r.content); f.flush()
        with h5py.File(f.name, "r") as h:
            grup = wt = None
            for ds in [k for k in h if k.startswith("dataset")]:
                for cami in (f"{ds}/data1/what", f"{ds}/what"):
                    if cami in h:
                        q = h[cami].attrs.get("quantity"); q = q.decode() if isinstance(q, bytes) else str(q)
                        if q == "ACRR":
                            grup, wt = f"{ds}/data1", h[cami].attrs; break
                if grup: break
            if grup is None:
                print("sense ACRR", clau); return None
            wh = h["where"].attrs
            pd = wh["projdef"].decode() if isinstance(wh["projdef"], bytes) else str(wh["projdef"])
            if pd not in _proj:
                _proj[pd] = Transformer.from_crs("EPSG:4326", CRS.from_proj4(pd), always_xy=True)
            tr = _proj[pd]
            xul, yul = tr.transform(float(wh["UL_lon"]), float(wh["UL_lat"]))
            X, Y = tr.transform(LO, LA)
            col = np.floor((X - xul) / float(wh["xscale"])).astype(int)
            row = np.floor((yul - Y) / float(wh["yscale"])).astype(int)
            r0, c0 = row.min(), col.min()
            v = h[grup + "/data"][r0:row.max() + 1, c0:col.max() + 1].astype(np.float64)[row - r0, col - c0]
            gain, off, nodata, undet = (float(wt[k]) for k in ("gain", "offset", "nodata", "undetect"))
    mm = v * gain + off
    mm[v == undet] = 0.0
    mm[v == nodata] = np.nan
    return np.clip(mm, 0, 160)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mes", required=True, help="AAAA-MM")
    a = ap.parse_args()
    any_, mes = map(int, a.mes.split("-"))
    ini = datetime(any_, mes, 1, tzinfo=timezone.utc)
    ndies = calendar.monthrange(any_, mes)[1]
    hores = [ini + timedelta(hours=i + 1) for i in range(ndies * 24)]     # hora que acaba a t
    claus = {}
    for d in range(ndies + 1):
        claus.update(claus_dia(ini + timedelta(days=d)))
    feina = [(i, claus.get(t)) for i, t in enumerate(hores)]
    cub = np.full((len(hores), len(LAT), len(LON)), 255, np.uint8)
    falten = []
    def fes(it):
        i, k = it
        return i, (llegeix(k) if k else None)
    with ThreadPoolExecutor(8) as ex:
        for i, mm in ex.map(fes, feina):
            if mm is None:
                falten.append(hores[i].strftime("%Y-%m-%dT%H:%M")); continue
            q = np.round(np.sqrt(mm) * 20)
            q[np.isnan(mm)] = 255
            cub[i] = np.clip(q, 0, 254).astype(np.uint8) * ~np.isnan(mm) + 255 * np.isnan(mm)
            if i % 48 == 0: print(hores[i], "fet")
    d = Path("data/pluja_cat"); d.mkdir(parents=True, exist_ok=True)
    (d / f"{a.mes}.bin.gz").write_bytes(gzip.compress(cub.tobytes(), 9))
    idx_f = d / "index.json"
    idx = json.loads(idx_f.read_text()) if idx_f.exists() else {"mesos": {}}
    idx.update({"oest": OEST, "est": EST, "sud": SUD, "nord": NORD, "pas": PAS, "files": len(LAT), "columnes": len(LON),
                "codificacio": "uint8 [hora][fila N→S][columna O→E]; mm = (v/20)^2; 255 = sense dada",
                "font": "EUMETNET OPERA Open Radar Data (CC BY 4.0), acumulat d'1 h ACRR, sense correcció amb estacions"})
    idx["mesos"][a.mes] = {"primera": hores[0].strftime("%Y-%m-%dT%H:%M"), "hores": len(hores), "falten": sorted(falten),
                           "generat": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M")}
    idx_f.write_text(json.dumps(idx, ensure_ascii=False, indent=1))
    tot = np.where(cub < 255, (cub.astype(float) / 20) ** 2, 0).sum(0)
    print(f"{a.mes}: {len(hores) - len(falten)}/{len(hores)} hores; màxim mensual {tot.max():.0f} mm; fitxer {(d / f'{a.mes}.bin.gz').stat().st_size/1e6:.1f} MB")


if __name__ == "__main__":
    main()
