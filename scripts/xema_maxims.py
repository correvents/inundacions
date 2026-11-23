#!/usr/bin/env python3
"""Màxims anuals de pluja diària de totes les estacions de la XEMA (Meteocat, dades obertes), per al cabal de disseny
dels passos de tot Catalunya (scripts/passos_cat.py, norma 5.2-IC).

Dades: analisi.transparenciacatalunya.cat, estacions yqwd-vj5e, lectures nzvn-apee (variable 35 = pluja en 30 min).
La suma de cada dia la fa el servidor (SoQL, $group per estació i dia). Un dia compta si té ≥ 44 lectures (≥ 22 h)
i un any si té ≥ 330 dies vàlids: el mateix criteri que p24_t100() a passos_cat.py.

Sortida: data/catalunya/xema/maxims_diaris.json  {codi: {nom, lat, lon, anys: {AAAA: mm}}}

Ús:  python3 scripts/xema_maxims.py --des 2009 --fins 2025"""
import argparse, json, time
from pathlib import Path
import requests

SOC = "https://analisi.transparenciacatalunya.cat/resource/"
SORT = Path("data/catalunya/xema/maxims_diaris.json")


def socrata(ds, params):
    for i in range(6):
        try:
            r = requests.get(SOC + ds + ".json", params=params, timeout=300); r.raise_for_status(); return r.json()
        except Exception as e:
            print("  reintent", ds, type(e).__name__, str(e)[:200]); time.sleep(15 * (i + 1))
    raise SystemExit(f"No respon {ds}")


def any_(y):
    """Pluja de cada dia i estació de l'any y: {codi: [(dia, mm, n)]}."""
    out, off = {}, 0
    while True:
        r = socrata("nzvn-apee", {
            "$select": "codi_estacio, date_trunc_ymd(data_lectura) AS dia, sum(valor_lectura) AS mm, count(*) AS n",
            "$where": f"codi_variable = '35' AND data_lectura >= '{y}-01-01T00:00:00' AND data_lectura < '{y + 1}-01-01T00:00:00'",
            "$group": "codi_estacio, dia", "$order": "codi_estacio, dia", "$limit": 50000, "$offset": off})
        for x in r:
            out.setdefault(x["codi_estacio"], []).append((x["dia"][:10], float(x.get("mm") or 0), int(x.get("n") or 0)))
        off += len(r)
        if len(r) < 50000:
            return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--des", type=int, default=2009)
    ap.add_argument("--fins", type=int, default=2025)
    a = ap.parse_args()
    est = {e["codi_estacio"]: e for e in socrata("yqwd-vj5e", {"$limit": 5000})}
    res = json.loads(SORT.read_text()) if SORT.exists() else {}
    for y in range(a.des, a.fins + 1):
        dies = any_(y)
        n_ok = 0
        for c, ll in dies.items():
            bons = [mm for _, mm, n in ll if n >= 44]
            if len(bons) < 330 or c not in est:
                continue
            e = est[c]
            r = res.setdefault(c, {"nom": e.get("nom_estacio", c), "lat": float(e["latitud"]), "lon": float(e["longitud"]), "anys": {}})
            r["anys"][str(y)] = round(max(bons), 1)
            n_ok += 1
        print(f"{y}: {len(dies)} estacions amb dades, {n_ok} amb l'any complet")
        SORT.parent.mkdir(parents=True, exist_ok=True)
        SORT.write_text(json.dumps(res, ensure_ascii=False, separators=(",", ":")))   # desa a cada any (si es talla, no es perd)
    print(f"{len(res)} estacions, {sum(len(r['anys']) for r in res.values())} anys-estació → {SORT}")


if __name__ == "__main__":
    main()
