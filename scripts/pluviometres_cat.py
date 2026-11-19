#!/usr/bin/env python3
"""Pluja horària de tots els pluviòmetres de Catalunya d'un mes, alineada amb la reixa de data/pluja_cat (hora que acaba a t, UTC).

Fonts (base): XEMA del Meteocat (dades obertes, pluja en 30 min, sense clau) i pluviòmetres de l'ACA (API Sentilo, sense clau,
només els ~3 últims mesos). AEMET s'hi afegirà quan hi hagi la clau (secret AEMET_API_KEY).
Escriu data/pluja_cat/AAAA-MM_estacions.json: [{id, font, nom, lat, lon, mm: [H valors o null]}].

Ús:  python3 scripts/pluviometres_cat.py --mes 2026-10 [--sense-aca]"""
import argparse, calendar, json, re, sys, time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from pathlib import Path
import numpy as np, requests
sys.path.insert(0, str(Path(__file__).parent))
import aca_arxiu as aca

SOC = "https://analisi.transparenciacatalunya.cat/resource/"
OEST, EST, SUD, NORD = 0.15, 3.35, 40.5, 42.9


def socrata(ds, params):
    for i in range(6):
        try:
            r = requests.get(SOC + ds + ".json", params=params, timeout=180); r.raise_for_status(); return r.json()
        except Exception as e:
            print("  reintent", ds, type(e).__name__); time.sleep(10 * (i + 1))
    raise SystemExit(f"No respon {ds}")


def xema(ini, H):
    est = {e["codi_estacio"]: e for e in socrata("yqwd-vj5e", {"$limit": 5000})}
    sums = {}
    fi = ini + timedelta(hours=H)
    # la lectura de les T cobreix [T, T+30 min): va a l'hora que acaba a T+30 arrodonida amunt
    dia = ini - timedelta(hours=1)
    while dia < fi:
        seg = min(dia + timedelta(days=4), fi)
        off = 0
        while True:
            r = socrata("nzvn-apee", {"codi_variable": "35", "$select": "codi_estacio,data_lectura,valor_lectura",
                                      "$where": f"data_lectura >= '{dia:%Y-%m-%dT%H:%M:%S}' AND data_lectura < '{seg:%Y-%m-%dT%H:%M:%S}'",
                                      "$order": "data_lectura,codi_estacio", "$limit": 50000, "$offset": off})
            for x in r:
                t = datetime.fromisoformat(x["data_lectura"][:19]).replace(tzinfo=timezone.utc) + timedelta(minutes=30)
                h = int(np.ceil((t - ini).total_seconds() / 3600)) - 1        # índex de l'hora que acaba a ini+(h+1)h
                if 0 <= h < H:
                    s = sums.setdefault(x["codi_estacio"], [[0.0, 0] for _ in range(H)])
                    s[h][0] += float(x["valor_lectura"]); s[h][1] += 1
            off += len(r)
            if len(r) < 50000: break
        print(f"  XEMA fins {seg:%Y-%m-%d}: {len(sums)} estacions")
        dia = seg
    out = []
    for c, s in sums.items():
        e = est.get(c)
        if not e: continue
        out.append({"id": "XEMA-" + c, "font": "XEMA", "nom": e.get("nom_estacio", c), "lat": float(e["latitud"]), "lon": float(e["longitud"]),
                    "mm": [round(v, 1) if n == 2 else None for v, n in s]})
    return out


def horaria(m, ini, H):
    """{ms: valor} → mm per hora. Detecta si el sensor dona pluja acumulada (creixent amb reinicis) o per interval."""
    ts = sorted(m)
    v = []
    for t in ts:
        try: v.append(float(m[t]))
        except ValueError: v.append(np.nan)
    t = np.array(ts) / 1000.0; v = np.array(v)
    ok = np.isfinite(v); t, v = t[ok], v[ok]
    if len(v) < 10: return None
    d = np.diff(v)
    acumulat = (d >= -1e-6).mean() > 0.97 and v.max() > 20 * max(np.median(v[v > 0]) if (v > 0).any() else 1, 0.1)
    if acumulat:
        inc = np.where(d >= 0, d, 0.0); t = t[1:]
    else:
        inc = v
    out = np.full(H, np.nan); n = np.zeros(H)
    h = np.ceil((t - ini.timestamp()) / 3600).astype(int) - 1
    s = np.zeros(H)
    for hh, x in zip(h, inc):
        if 0 <= hh < H: s[hh] += x; n[hh] += 1
    # una hora val si té prou mesures (la majoria de sensors donen cada 5 o 10 min)
    tip = np.median(n[n > 0]) if (n > 0).any() else 0
    out[n >= max(1, 0.75 * tip)] = s[n >= max(1, 0.75 * tip)]
    out[(out < 0) | (out > 200)] = np.nan
    return out, acumulat


def aca_pluvios(ini, H):
    s = aca.sessio()
    cat = aca.get(s, "/catalog")
    if not cat: print("ACA: sense catàleg"); return []
    sens = []
    for p in cat.get("providers", []):
        for se in p.get("sensors", []):
            txt = " ".join(str(se.get(k, "")) for k in ("type", "sensorType", "unit", "description", "componentType", "componentDesc"))
            if re.search(r"plu|precip|rain", txt, re.I) and not re.search(r"temp|humit|vent|wind", txt, re.I):
                loc = str(se.get("location") or "").split()
                if len(loc) > 1:
                    la, lo = float(loc[0]), float(loc[1])
                    if SUD <= la <= NORD and OEST <= lo <= EST:
                        sens.append((p["provider"], se, la, lo, txt))
    print(f"ACA: {len(sens)} pluviòmetres al catàleg dins de la finestra")
    fi = ini + timedelta(hours=H)
    def baixa(x):
        p, se, la, lo, txt = x
        m = aca.llegeix(aca.sessio(), p, se["sensor"], ini - timedelta(hours=1), fi)
        r = horaria(m, ini, H) if m else None
        if r is None: return None
        mm, acum = r
        return {"id": f"ACA-{se['sensor']}", "font": "ACA", "nom": se.get("componentDesc") or se.get("description") or se["sensor"],
                "lat": la, "lon": lo, "acumulat": bool(acum), "mm": [None if not np.isfinite(v) else round(float(v), 1) for v in mm]}
    out = []
    with ThreadPoolExecutor(8) as ex:
        for k, r in enumerate(ex.map(baixa, sens)):
            if r and sum(v is not None for v in r["mm"]) > 24: out.append(r)
    print(f"ACA: {len(out)} pluviòmetres amb dades ({sum(r['acumulat'] for r in out)} en format acumulat)")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mes", required=True)
    ap.add_argument("--sense-aca", action="store_true")
    a = ap.parse_args()
    y, m = map(int, a.mes.split("-"))
    ini = datetime(y, m, 1, tzinfo=timezone.utc); H = calendar.monthrange(y, m)[1] * 24
    est = xema(ini, H)
    if not a.sense_aca and (datetime.now(timezone.utc) - ini).days < 95:
        est += aca_pluvios(ini, H)
    Path("data/pluja_cat").mkdir(parents=True, exist_ok=True)
    Path(f"data/pluja_cat/{a.mes}_estacions.json").write_text(json.dumps(est, ensure_ascii=False, separators=(",", ":")))
    print(f"{a.mes}: {len(est)} pluviòmetres (" + ", ".join(f"{f} {sum(e['font'] == f for e in est)}" for f in ("XEMA", "ACA")) + ")")


if __name__ == "__main__":
    main()
