#!/usr/bin/env python3
"""Extreu d'un arxiu de l'ACA (aca_arxiu.py) la sèrie d'un sensor entre dues dates i la desa com a JSON per a evolucio.html.

Sortida: {"font", "sensor", "temps_utc": [...], "q_aforador_obs": [...]} (valors tal com els dona l'ACA, normalment m³/s).

Ús:  python3 scripts/extreu_aforador.py --arxiu /tmp/aca --provider AFORAMENT-EST --sensor CALC001212 \\
        --des 2026-10-03T18:00 --fins 2026-10-06T00:00 --sortida data/granollers/episodis/2026-10-04_aforador.json"""
import argparse, csv, gzip, json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arxiu", required=True); ap.add_argument("--provider", required=True); ap.add_argument("--sensor", required=True)
    ap.add_argument("--des", required=True); ap.add_argument("--fins", required=True); ap.add_argument("--sortida", required=True)
    a = ap.parse_args()
    t, v = [], []
    with gzip.open(Path(a.arxiu) / a.provider / f"{a.sensor}.csv.gz", "rt") as f:
        for r in csv.DictReader(f):
            if a.des <= r["temps_utc"][:16] <= a.fins:
                try: x = float(r["valor"])
                except ValueError: continue
                t.append(r["temps_utc"][:16]); v.append(x)
    Path(a.sortida).write_text(json.dumps({"font": f"ACA {a.provider}/{a.sensor}", "sensor": a.sensor, "temps_utc": t, "q_aforador_obs": v}))
    print(f"{len(t)} mesures, màx. {max(v) if v else None} el {t[v.index(max(v))] if v else ''}")
    return 0


if __name__ == "__main__":
    main()
