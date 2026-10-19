#!/usr/bin/env python3
"""Actuacions dels Bombers de la Generalitat (dades obertes g2ay-3vnj) als municipis de la plana del Ridaura.

Serveixen d'etiqueta independent de la premsa: un dia amb actuacions per aigua (inundacions, sanejaments, salvaments
en riera...) a Castell-Platja d'Aro és un dia que es va inundar alguna cosa, encara que no surti als diaris.
Escriu data/<conca>/model/bombers.csv (totes les actuacions) i bombers_aigua_dia.csv (nombre d'actuacions per aigua per dia).

Ús:  python3 scripts/bombers.py --config config/ridaura.json"""
import argparse, csv, json, re, sys, time
from collections import Counter
from pathlib import Path
import requests

API = "https://analisi.transparenciacatalunya.cat/resource/g2ay-3vnj.json"
MUNICIPIS = ["Castell-Platja d'Aro", "Santa Cristina d'Aro"]
AIGUA = re.compile(r"aigua|inund|pluj|tempest|riera|riu|achic|esgot|desguas|clavegu|filtrac|meteor", re.I)


def baixa(params):
    for i in range(5):
        try:
            r = requests.get(API, params=params, timeout=120); r.raise_for_status(); return r.json()
        except Exception as e:
            print("reintent", e); time.sleep(10 * (i + 1))
    raise SystemExit("No respon l'API")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    d = Path("data") / json.loads(Path(a.config).read_text())["conca"] / "model"; d.mkdir(parents=True, exist_ok=True)
    mostra = baixa({"$limit": 3})
    print("camps:", list(mostra[0].keys()) if mostra else "cap fila")
    files = []
    for m in MUNICIPIS:
        off = 0
        while True:
            r = baixa({"$where": f"upper(municipi) like upper('%{m.replace(chr(39), chr(39)*2)}%')", "$limit": 50000, "$offset": off,
                       "$order": "act_dat_actuacio"})
            files += r; off += len(r)
            if len(r) < 50000: break
        print(m, off, "actuacions")
    if not files:
        print("Cap actuació: revisa el nom del municipi"); return 1
    claus = sorted({k for f in files for k in f})
    with open(d / "bombers.csv", "w", newline="") as f:
        w = csv.DictWriter(f, claus); w.writeheader(); w.writerows(files)
    tipus = Counter(f.get("tal_nom_alarma", "") for f in files)
    print("tipus d'alarma (els 60 més freqüents):"); [print(f"  {n:6d}  {t}") for t, n in tipus.most_common(60)]
    print("grups:", Counter(f.get("tga_nom_grupo", "") for f in files).most_common(20))
    dies = Counter((f.get("act_dat_actuacio", "")[:10], f.get("municipi", "")) for f in files
                   if AIGUA.search(f.get("tal_nom_alarma", "") + " " + f.get("tga_nom_grupo", "")))
    with open(d / "bombers_aigua_dia.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["dia", "municipi", "actuacions"])
        for (dia, mun), n in sorted(dies.items()): w.writerow([dia, mun, n])
    print(len(dies), "dies amb actuacions per aigua; els 30 amb més:")
    [print(" ", k, n) for k, n in dies.most_common(30)]
    return 0


if __name__ == "__main__":
    sys.exit(main())
