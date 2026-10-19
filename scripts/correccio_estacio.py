#!/usr/bin/env python3
"""Sèrie per corregir el radar d'un episodi amb l'estació XEMA de Castell d'Aro (DO), per a simula_2d.py --episodi.

El radar OPERA sense corregir s'equivoca molt en tots dos sentits (gener 2026: 29 mm de radar per 162 mm a l'estació;
setembre 2023: més radar que estació). simula_2d.py corregeix cada fotograma amb el quocient pluviòmetre/radar de l'hora
anterior; aquí el «pluviòmetre» és l'estació i el «radar» la mitjana de les cel·les a menys de 3 km de l'estació.
Escriu data/<conca>/episodis/<id>_estacio.json (temps_utc cada 30 min, pluja_pluviometre_mm, pluja_radar_conca_mm).

Ús:  python3 scripts/correccio_estacio.py --config config/ridaura.json --radar 2019-09-20_radar.json --id 2019-09-21"""
import argparse, json, math, sys
from pathlib import Path
import numpy as np, pandas as pd

ESTACIO, LAT, LON, RADI_KM = "DO", 41.80855, 3.03241, 3.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--radar", required=True)
    ap.add_argument("--id", required=True)
    a = ap.parse_args()
    d = Path("data") / json.loads(Path(a.config).read_text())["conca"]
    R = json.loads((d / "episodis" / a.radar).read_text())
    LO, LA = np.meshgrid(R["lons"], R["lats"])
    km = np.hypot((LA - LAT) * 111, (LO - LON) * 111 * math.cos(math.radians(LAT))).ravel()
    rad = pd.Series(np.array(R["mm"])[:, km <= RADI_KM].mean(axis=1), index=pd.to_datetime(R["temps"]))
    # les lectures de la XEMA porten l'hora d'inici de la mitja hora; el radar, la de final del fotograma
    est = pd.read_csv(d / "xema" / f"{ESTACIO}.csv", parse_dates=["temps"]).set_index("temps").mm
    est.index = est.index + pd.Timedelta(minutes=30)
    rad30 = rad.resample("30min", label="right", closed="right").sum()
    t = rad30.index
    out = {"titol": f"Correcció del radar amb l'estació XEMA {ESTACIO} (Castell d'Aro)", "temps_utc": [x.strftime("%Y-%m-%dT%H:%M") for x in t],
           "pas_min": 30, "pluja_pluviometre_mm": [None if pd.isna(v) else round(float(v), 2) for v in est.reindex(t)],
           "pluja_radar_conca_mm": [round(float(v), 2) for v in rad30]}
    (d / "episodis" / f"{a.id}_estacio.json").write_text(json.dumps(out))
    print(f"{a.id}: estació {np.nansum([v or 0 for v in out['pluja_pluviometre_mm']]):.0f} mm, radar {rad30.sum():.0f} mm")
    return 0


if __name__ == "__main__":
    sys.exit(main())
