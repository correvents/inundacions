#!/usr/bin/env python3
"""Graella de pluja del radar OPERA d'un episodi passat, per animar-la al visor.

Mostreja els compostos DBZH (5 min, ~1 km) de l'arxiu en una reixa regular sobre la zona i desa
l'acumulat de cada interval (per defecte 15 min) en mm, a data/<conca>/episodis/<data>_radar.json.

Ús:  python3 scripts/radar_graella.py --config config/ridaura.json --des 2026-10-03T18:00 --fins 2026-10-04T14:00"""
import argparse, json, sys, tempfile, urllib.parse
from datetime import datetime, timedelta, timezone
from pathlib import Path
import numpy as np, requests
from radar_conca import S3_ARXIU, claus, durades, mostreja, mmh

PAS = 0.01      # graus (~1 km, la resolució del radar)
CAIXA = [2.84, 41.74, 3.12, 41.90]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--des", required=True)
    ap.add_argument("--fins", required=True)
    ap.add_argument("--interval", type=int, default=15, help="minuts de cada fotograma")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    ini = datetime.fromisoformat(a.des).replace(tzinfo=timezone.utc)
    fi = datetime.fromisoformat(a.fins).replace(tzinfo=timezone.utc)
    lons = np.round(np.arange(CAIXA[0], CAIXA[2] + 1e-9, PAS), 4)
    lats = np.round(np.arange(CAIXA[1], CAIXA[3] + 1e-9, PAS), 4)
    LON, LAT = np.meshgrid(lons, lats)
    disponibles = claus(S3_ARXIU, ini, fi)
    minuts = durades(disponibles)
    print(f"OPERA: {len(disponibles)} compostos (cada {sorted(set(minuts.values()))} min)")
    n = int((fi - ini).total_seconds() // (a.interval * 60))
    acum = np.zeros((n, LAT.size))
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "c.h5"
        for t, k in sorted(disponibles.items()):
            i = int((t - ini).total_seconds() // (a.interval * 60) - (1 if (t - ini).total_seconds() % (a.interval * 60) == 0 else 0))
            if not 0 <= i < n:
                continue
            try:
                r = requests.get(S3_ARXIU + "/" + urllib.parse.quote(k, safe="/@"), timeout=120); r.raise_for_status()
                f.write_bytes(r.content)
                acum[i] += np.nan_to_num(mmh(mostreja(f, LON.ravel(), LAT.ravel()))) * minuts[t] / 60
            except Exception as e:
                print(f"OPERA {t:%Y-%m-%d %H:%M}: {type(e).__name__} {e}")
    out = {"font": "EUMETNET OPERA (CC BY 4.0), DBZH, Z=200R^1.6", "pas_graus": PAS, "lons": lons.tolist(), "lats": lats.tolist(),
           "interval_min": a.interval,
           "temps": [(ini + timedelta(minutes=a.interval * (i + 1))).strftime("%Y-%m-%dT%H:%M") for i in range(n)],
           "mm": [[round(float(v), 1) for v in fila] for fila in acum]}
    dest = Path("data") / cfg["conca"] / "episodis"; dest.mkdir(parents=True, exist_ok=True)
    (dest / f"{a.des[:10]}_radar.json").write_text(json.dumps(out, separators=(",", ":")))
    print(f"{n} fotogrames, màx. {acum.max():.1f} mm en {a.interval} min")
    return 0


if __name__ == "__main__":
    sys.exit(main())
