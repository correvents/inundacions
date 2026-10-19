#!/usr/bin/env python3
"""Sèrie llarga de pluja de les estacions de la XEMA (Meteocat, dades obertes) al voltant de la conca, i catàleg d'episodis de pluja forta.

Serveix per saber totes les vegades que ha plogut fort (també les que no van inundar), per creuar-les amb les notícies.
Dades: analisi.transparenciacatalunya.cat, estacions yqwd-vj5e, lectures nzvn-apee (variable 35 = precipitació en 30 min).
Escriu data/<conca>/xema/<codi>.csv (temps UTC, mm en 30 min), estacions.json i episodis_pluja.csv:
un episodi és un període amb pluja separat de l'anterior per 12 h seques i amb ≥ 30 mm o ≥ 20 mm/h.

Ús:  python3 scripts/pluja_xema.py --config config/ridaura.json [--radi-km 20]"""
import argparse, json, math, sys, time, csv
from pathlib import Path
from datetime import datetime, timedelta
import requests

API = "https://analisi.transparenciacatalunya.cat/resource/"


def socrata(ds, params):
    for i in range(5):
        try:
            r = requests.get(API + ds + ".json", params=params, timeout=120); r.raise_for_status(); return r.json()
        except Exception as e:
            print("reintent", e); time.sleep(10 * (i + 1))
    raise SystemExit(f"No respon {ds}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--radi-km", type=float, default=20)
    ap.add_argument("--des", help="només les lectures des d'aquesta data (AAAA-MM-DD); per a conques noves en zones amb moltes estacions")
    ap.add_argument("--max-estacions", type=int, default=0, help="només les N més properes (0 = totes)")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    d = Path("data") / cfg["conca"] / "xema"; d.mkdir(parents=True, exist_ok=True)
    lat0, lon0 = cfg.get("centre", [41.81, 3.03])       # per defecte, la plana de Platja d'Aro
    est = socrata("yqwd-vj5e", {"$limit": 5000})
    km = lambda e: math.hypot((float(e["latitud"]) - lat0) * 111, (float(e["longitud"]) - lon0) * 111 * math.cos(math.radians(lat0)))
    prop = sorted([dict(e, km=round(km(e), 1)) for e in est if km(e) <= a.radi_km], key=lambda e: e["km"])
    if a.max_estacions:
        prop = prop[:a.max_estacions]
    (d / "estacions.json").write_text(json.dumps(prop, ensure_ascii=False, indent=1))
    print("estacions:", [(e["codi_estacio"], e["nom_estacio"], e["km"], e.get("data_inici"), e.get("data_fi")) for e in prop])
    episodis = []
    for e in prop:
        codi, files, off = e["codi_estacio"], [], 0
        while True:
            q = {"codi_estacio": codi, "codi_variable": "35", "$select": "data_lectura,valor_lectura", "$order": "data_lectura", "$limit": 50000, "$offset": off}
            if a.des:
                q["$where"] = f"data_lectura >= '{a.des}T00:00:00.000'"
            r = socrata("nzvn-apee", q)
            files += r; off += len(r)
            if len(r) < 50000: break
        if not files: continue
        with open(d / f"{codi}.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["temps", "mm"])
            for x in files: w.writerow([x["data_lectura"][:16], x["valor_lectura"]])
        print(codi, e["nom_estacio"], len(files), "lectures", files[0]["data_lectura"][:10], "→", files[-1]["data_lectura"][:10])
        # episodis
        cur = None
        for x in files:
            t, v = datetime.fromisoformat(x["data_lectura"][:16]), float(x["valor_lectura"])
            if v <= 0: continue
            if cur and t - cur["fi"] > timedelta(hours=12):
                episodis.append(cur); cur = None
            if not cur: cur = {"estacio": codi, "nom": e["nom_estacio"], "km": e["km"], "inici": t, "fi": t, "mm": 0, "max30": 0, "serie": []}
            cur["fi"] = t; cur["mm"] += v; cur["max30"] = max(cur["max30"], v); cur["serie"].append((t, v))
        if cur: episodis.append(cur)
    fort = []
    for p in episodis:
        s = p.pop("serie")
        def maxh(h):
            best = 0
            for i, (t, _) in enumerate(s):
                best = max(best, sum(v for tt, v in s[i:] if tt - t < timedelta(hours=h)))
            return best
        p["max1h"], p["max3h"], p["max6h"], p["max24h"] = maxh(1), maxh(3), maxh(6), maxh(24)
        if p["mm"] >= 30 or p["max1h"] >= 20:
            fort.append(p)
    with open(d / "episodis_pluja.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["estacio", "nom", "km", "inici_utc", "fi_utc", "mm", "max30min", "max1h", "max3h", "max6h", "max24h"])
        for p in sorted(fort, key=lambda p: p["inici"]):
            w.writerow([p["estacio"], p["nom"], p["km"], p["inici"].isoformat(), p["fi"].isoformat(),
                        round(p["mm"], 1), round(p["max30"], 1), round(p["max1h"], 1), round(p["max3h"], 1), round(p["max6h"], 1), round(p["max24h"], 1)])
    print(len(fort), "episodis de pluja forta")
    return 0


if __name__ == "__main__":
    sys.exit(main())
