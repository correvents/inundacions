#!/usr/bin/env python3
"""Estat del mar davant de Platja d'Aro (onatge i nivell del mar) per a cada episodi de pluja.

Font: Open-Meteo Marine API (reanàlisi i models d'onatge; nivell del mar amb marea i sobreelevació quan n'hi ha).
Hipòtesi a comprovar: amb mar alta o llevantada forta l'aigua de pluja no surt pel clavegueram i les rieres, i hi ha
més inundació amb la mateixa pluja.
Escriu data/<conca>/model/mar_horari.csv (temps UTC, onatge m, direcció, període, nivell del mar m).

Ús:  python3 scripts/mar.py --config config/ridaura.json"""
import argparse, json, sys, time
from pathlib import Path
import requests

LAT, LON = 41.805, 3.075
VARS = "wave_height,wave_direction,wave_period,sea_level_height_msl"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--des", default="2009-01-01")
    ap.add_argument("--fins", default="2026-10-18")
    a = ap.parse_args()
    d = Path("data") / json.loads(Path(a.config).read_text())["conca"] / "model"
    files = []
    for any_ in range(int(a.des[:4]), int(a.fins[:4]) + 1):
        ini, fi = max(f"{any_}-01-01", a.des), min(f"{any_}-12-31", a.fins)
        for i in range(4):
            try:
                r = requests.get("https://marine-api.open-meteo.com/v1/marine", timeout=120, params={
                    "latitude": LAT, "longitude": LON, "hourly": VARS, "start_date": ini, "end_date": fi, "timezone": "GMT"})
                r.raise_for_status(); break
            except Exception as e:
                print(any_, "reintent", e, getattr(e, "response", None) and e.response.text[:300]); time.sleep(10 * (i + 1))
        else:
            continue
        h = r.json()["hourly"]
        n = sum(v is not None for v in h["wave_height"]); m = sum(v is not None for v in h["sea_level_height_msl"])
        print(any_, len(h["time"]), "hores; amb onatge", n, "amb nivell", m)
        files += list(zip(h["time"], h["wave_height"], h["wave_direction"], h["wave_period"], h["sea_level_height_msl"]))
    with open(d / "mar_horari.csv", "w") as f:
        f.write("temps,onatge_m,direccio,periode_s,nivell_m\n")
        for t, *v in files:
            f.write(t + "," + ",".join("" if x is None else str(x) for x in v) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
