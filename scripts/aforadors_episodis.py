#!/usr/bin/env python3
"""Pic de cabal de cada aforador de l'ACA en els episodis de pluja de 2026 (per separar desbordaments de riu).

Llegeix l'arxiu de l'ACA desat al bucket privat de Supabase (aca/<dia>/<provider>.tar, fet per aca_arxiu.py) i, per als
episodis més recents que l'arxiu, l'API de l'ACA en directe (només sensors de cabal). Per a cada sensor de cabal i episodi
escriu el màxim, l'hora del màxim, el cabal de base (mediana de tota la sèrie) i el quocient pic/base.
Sortida: data/catalunya/aforadors/episodis_2026.json (només dades públiques de l'ACA).

Ús:  SUPABASE_URL=... SUPABASE_SECRET_KEY=... python3 scripts/aforadors_episodis.py"""
import csv, gzip, io, json, os, sys, tarfile, statistics
from datetime import datetime, timedelta, timezone
from pathlib import Path
import requests

sys.path.insert(0, str(Path(__file__).parent))
import aca_arxiu as A

EPISODIS = [  # inici i fi (dies UTC) dels episodis de pluja amb afectacions el 2026
    ("2026-08-24", "2026-08-25"), ("2026-09-09", "2026-09-09"), ("2026-09-16", "2026-09-17"), ("2026-09-29", "2026-09-30"),
    ("2026-10-01", "2026-10-02"), ("2026-10-03", "2026-10-07"), ("2026-10-17", "2026-10-17"), ("2026-10-19", "2026-10-20"),
    ("2026-11-11", "2026-11-14"),
]
ARXIUS = ["2026-10-19", "2026-11-01"]
DEST = Path("data/catalunya/aforadors/episodis_2026.json")


def es_cabal(r):
    t = " ".join(str(r.get(k, "")) for k in ("tipus", "unitat", "descripcio")).lower()
    return ("m3/s" in t or "m³/s" in t or "cabal" in t or "caudal" in t or "flow" in t) and "embass" not in t


def storage(nom):
    url = os.environ["SUPABASE_URL"].rstrip("/") + "/storage/v1/object/arxiu-inundacions/" + nom
    k = os.environ["SUPABASE_SECRET_KEY"]
    r = requests.get(url, headers={"Authorization": f"Bearer {k}", "apikey": k}, timeout=300)
    r.raise_for_status(); return r.content


def main():
    series, meta = {}, {}
    for dia in ARXIUS:
        txt = storage(f"aca/{dia}/sensors.csv").decode()
        files = [r for r in csv.DictReader(io.StringIO(txt)) if es_cabal(r)]
        for r in files: meta[(r["provider"], r["sensor"])] = r
        provs = sorted({r["provider"] for r in files})
        for p in provs:
            try:
                tf = tarfile.open(fileobj=io.BytesIO(storage(f"aca/{dia}/{p}.tar")))
            except Exception as e:
                print("no hi ha", dia, p, e); continue
            for m in tf.getmembers():
                s = Path(m.name).name.replace(".csv.gz", "")
                if (p, s) not in meta: continue
                d = series.setdefault((p, s), {})
                for lin in gzip.decompress(tf.extractfile(m).read()).decode().splitlines()[1:]:
                    t, v = lin.split(",", 1)
                    try: d[t] = float(v)
                    except ValueError: pass
        print(dia, len(files), "sensors de cabal")
    # episodis posteriors a l'últim arxiu: API en directe, només els sensors de cabal
    ini = datetime(2026, 10, 30, tzinfo=timezone.utc); fi = datetime.now(timezone.utc)
    s = A.sessio()
    for (p, se) in list(meta):
        m = A.llegeix(s, p, se, ini, fi)
        d = series.setdefault((p, se), {})
        for ms, v in m.items():
            try: d[datetime.fromtimestamp(ms / 1000, timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")] = float(v)
            except ValueError: pass
    out = []
    for k, d in series.items():
        if len(d) < 100: continue
        r = meta[k]; vals = [v for v in d.values() if v >= 0]
        base = statistics.median(vals) if vals else None
        eps = []
        for a, b in EPISODIS:
            t0 = (datetime.fromisoformat(a) - timedelta(days=1)).strftime("%Y-%m-%d"); t1 = (datetime.fromisoformat(b) + timedelta(days=2)).strftime("%Y-%m-%d")
            w = [(t, v) for t, v in d.items() if t0 <= t[:10] <= t1 and v >= 0]
            if not w: eps.append(None); continue
            t, v = max(w, key=lambda x: x[1])
            eps.append({"pic": round(v, 2), "hora": t[:16], "n": len(w)})
        out.append({"provider": k[0], "sensor": k[1], "nom": r.get("descripcio"), "unitat": r.get("unitat"),
                    "lat": float(r["lat"]) if r.get("lat") else None, "lon": float(r["lon"]) if r.get("lon") else None,
                    "base_mediana": round(base, 3) if base is not None else None, "max_serie": round(max(vals), 2) if vals else None,
                    "primera": min(d)[:16], "ultima": max(d)[:16], "episodis": eps})
    DEST.parent.mkdir(parents=True, exist_ok=True)
    DEST.write_text(json.dumps({"episodis": EPISODIS, "font": "ACA (Sentilo), arxiu propi + API", "aforadors": out}, ensure_ascii=False, indent=0))
    print(len(out), "aforadors")
    for o in sorted(out, key=lambda o: -max((e or {}).get("pic", 0) / max(o["base_mediana"] or 0.1, 0.1) for e in o["episodis"]))[:25]:
        print(o["nom"], o["base_mediana"], [(e or {}).get("pic") for e in o["episodis"]])


if __name__ == "__main__":
    main()
