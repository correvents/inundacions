#!/usr/bin/env python3
"""Arxiva totes les mesures de l'ACA a Catalunya (aforadors, pluviòmetres, embassaments...) abans que caduquin.

L'API Sentilo de l'ACA només guarda uns 3 mesos. Aquest script llegeix el catàleg sencer, filtra els sensors
hidrològics i en descarrega tota la finestra disponible. Sortida (--dest):
  cataleg.json                       el catàleg tal com el dona l'API
  sensors.csv                        provider, sensor, tipus, unitat, component, descripció, lat, lon, mesures, primera, última
  <provider>/<sensor>.csv.gz         temps_utc, valor

Ús:  python3 scripts/aca_arxiu.py --dest /tmp/aca [--dies 100] [--fils 8]"""
import argparse, csv, gzip, json, re, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
from pathlib import Path
import requests

API = "https://aplicacions.aca.gencat.cat/sdim2/apirest"
FMT = "%d/%m/%YT%H:%M:%S"
# sensors hidrològics: nivell, cabal, pluja, embassament, piezòmetre...
INTERES = re.compile(r"niv|cabal|caudal|flow|level|plu|precip|rain|emb|volum|piez|aforament|calc", re.I)


def sessio():
    s = requests.Session()
    s.headers.update({"Accept": "application/json", "User-Agent": "inundacions/0.1 (github.com/correvents/inundacions)"})
    return s


def get(s, ruta, params=None, intents=4):
    for i in range(intents):
        try:
            r = s.get(API + ruta, params=params, timeout=90)
            if r.status_code == 200 and r.text.lstrip().startswith(("{", "[")):
                return r.json()
            print(f"  {ruta}: HTTP {r.status_code}", file=sys.stderr)
        except (requests.RequestException, ValueError) as e:
            print(f"  {ruta}: {type(e).__name__}", file=sys.stderr)
        time.sleep(4 * (i + 1))
    return None


def llegeix(s, provider, sensor, ini, fi):
    """Totes les mesures entre ini i fi; l'API en dona com a màxim 200 per petició, es pagina enrere."""
    out, fins, ant = {}, fi, None
    for _ in range(2000):
        j = get(s, f"/data/{provider}/{sensor}", {"limit": 200, "from": ini.strftime(FMT), "to": fins.strftime(FMT)})
        ob = (j or {}).get("observations", [])
        if not ob:
            break
        for o in ob:
            try:
                out[int(o["time"])] = str(o.get("value")).strip()
            except (KeyError, ValueError):
                pass
        antic = min(int(o["time"]) for o in ob)
        if len(ob) < 200 or antic <= ini.timestamp() * 1000 or (ant is not None and antic >= ant):
            break
        ant, fins = antic, datetime.fromtimestamp(antic / 1000 - 1, timezone.utc)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dest", required=True)
    ap.add_argument("--dies", type=float, default=100)
    ap.add_argument("--fils", type=int, default=8)
    ap.add_argument("--max", type=int, default=0, help="només els N primers sensors (proves)")
    a = ap.parse_args()
    dest = Path(a.dest); dest.mkdir(parents=True, exist_ok=True)
    s = sessio()
    cat = get(s, "/catalog")
    if not cat:
        print("No s'ha pogut llegir el catàleg"); return 1
    (dest / "cataleg.json").write_text(json.dumps(cat, ensure_ascii=False))
    sensors = []
    for p in cat.get("providers", []):
        for se in p.get("sensors", []):
            txt = " ".join(str(se.get(k, "")) for k in ("type", "sensorType", "unit", "description", "componentType", "componentDesc"))
            if INTERES.search(txt + " " + p.get("provider", "")):
                sensors.append((p.get("provider"), se))
    tipus = {}
    for p, se in sensors:
        tipus[(p, se.get("type"))] = tipus.get((p, se.get("type")), 0) + 1
    print(f"Catàleg: {sum(len(p.get('sensors', [])) for p in cat.get('providers', []))} sensors, {len(sensors)} hidrològics")
    for k, n in sorted(tipus.items(), key=lambda x: -x[1]):
        print(f"  {k[0]} / {k[1]}: {n}")
    if a.max:
        sensors = sensors[:a.max]
    fi = datetime.now(timezone.utc); ini = fi - timedelta(days=a.dies)
    files, t0 = [], time.time()

    def baixa(p, se):
        m = llegeix(sessio(), p, se["sensor"], ini, fi)
        if m:
            d = dest / p; d.mkdir(exist_ok=True)
            with gzip.open(d / f"{se['sensor']}.csv.gz", "wt") as f:
                f.write("temps_utc,valor\n")
                for ms in sorted(m):
                    f.write(f"{datetime.fromtimestamp(ms / 1000, timezone.utc):%Y-%m-%dT%H:%M:%S},{m[ms]}\n")
        loc = str(se.get("location") or "").split()
        return [p, se["sensor"], se.get("type"), se.get("unit"), se.get("component"), se.get("componentDesc") or se.get("description"),
                loc[0] if len(loc) > 1 else "", loc[1] if len(loc) > 1 else "", len(m),
                datetime.fromtimestamp(min(m) / 1000, timezone.utc).isoformat() if m else "",
                datetime.fromtimestamp(max(m) / 1000, timezone.utc).isoformat() if m else ""]

    with ThreadPoolExecutor(a.fils) as ex:
        fut = [ex.submit(baixa, p, se) for p, se in sensors]
        for k, f in enumerate(as_completed(fut)):
            try:
                files.append(f.result())
            except Exception as e:
                print(f"  error: {type(e).__name__} {e}", file=sys.stderr)
            if k % 50 == 0:
                print(f"{k + 1}/{len(fut)} sensors ({time.time() - t0:.0f} s)", flush=True)
    with open(dest / "sensors.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["provider", "sensor", "tipus", "unitat", "component", "descripcio", "lat", "lon", "mesures", "primera", "ultima"])
        w.writerows(sorted(files))
    print(f"Fet: {sum(1 for r in files if r[8])} sensors amb dades, {sum(r[8] for r in files)} mesures ({time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
