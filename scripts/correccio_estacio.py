#!/usr/bin/env python3
"""Sèrie per corregir el radar d'un episodi amb les estacions XEMA properes, per a simula_2d.py --episodi.

El radar OPERA sense corregir s'equivoca molt en tots dos sentits (gener 2026: 29 mm de radar per 162 mm a l'estació;
setembre 2023: més radar que estació). simula_2d.py corregeix cada fotograma amb el quocient pluviòmetre/radar de l'hora
anterior. Aquí el «pluviòmetre» és la suma de les estacions XEMA a menys de 20 km amb dades (Castell d'Aro, Cassà,
Palafrugell, la Bisbal, Monells) i el «radar» la suma del radar sobre cadascuna (cel·les a menys de 3 km): és el biaix
mitjà del camp (mean field bias), més estable que una sola estació. Amb --nomes-do es fa servir només Castell d'Aro.
Escriu data/<conca>/episodis/<id>_estacio.json (temps_utc cada 30 min, pluja_pluviometre_mm, pluja_radar_conca_mm).

Ús:  python3 scripts/correccio_estacio.py --config config/ridaura.json --radar 2019-09-20_radar.json --id 2019-09-21"""
import argparse, json, math, sys
from pathlib import Path
import numpy as np, pandas as pd

RADI_KM = 3.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--radar", required=True)
    ap.add_argument("--id", required=True)
    ap.add_argument("--nomes-do", action="store_true")
    a = ap.parse_args()
    d = Path("data") / json.loads(Path(a.config).read_text())["conca"]
    R = json.loads((d / "episodis" / a.radar).read_text())
    LO, LA = np.meshgrid(R["lons"], R["lats"]); mm = np.array(R["mm"]); t_rad = pd.to_datetime(R["temps"])
    estacions = [e for e in json.loads((d / "xema" / "estacions.json").read_text()) if (d / "xema" / f"{e['codi_estacio']}.csv").exists()]
    if a.nomes_do:
        estacions = [e for e in estacions if e["codi_estacio"] == "DO"]
    t = None; plv = rdr = None; usades = []
    for e in estacions:
        lat, lon = float(e["latitud"]), float(e["longitud"])
        km = np.hypot((LA - lat) * 111, (LO - lon) * 111 * math.cos(math.radians(lat))).ravel()
        if not (km <= RADI_KM).any():
            continue
        rad30 = pd.Series(mm[:, km <= RADI_KM].mean(axis=1), index=t_rad).resample("30min", label="right", closed="right").sum()
        # les lectures de la XEMA porten l'hora d'inici de la mitja hora; el radar, la de final del fotograma
        est = pd.read_csv(d / "xema" / f"{e['codi_estacio']}.csv", parse_dates=["temps"]).set_index("temps").mm
        est.index = est.index + pd.Timedelta(minutes=30)
        est = est.reindex(rad30.index)
        if est.notna().mean() < 0.8:
            continue                                                   # l'estació no funcionava durant l'episodi
        ok = est.notna()
        t = rad30.index if t is None else t
        plv = est.fillna(0) if plv is None else plv + est.fillna(0)
        rdr = rad30.where(ok, 0) if rdr is None else rdr + rad30.where(ok, 0)
        usades.append(f"{e['codi_estacio']} {e['nom_estacio']} {est.sum():.0f}/{rad30.sum():.0f} mm")
    if t is None:
        sys.exit("Cap estació amb dades durant l'episodi")
    est, rad30 = plv, rdr
    out = {"titol": "Correcció del radar amb estacions XEMA: " + "; ".join(usades), "estacions": usades, "temps_utc": [x.strftime("%Y-%m-%dT%H:%M") for x in t],
           "pas_min": 30, "pluja_pluviometre_mm": [None if pd.isna(v) else round(float(v), 2) for v in est.reindex(t)],
           "pluja_radar_conca_mm": [round(float(v), 2) for v in rad30]}
    (d / "episodis" / f"{a.id}_estacio.json").write_text(json.dumps(out))
    print(f"{a.id}: {'; '.join(usades)} (estació/radar)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
