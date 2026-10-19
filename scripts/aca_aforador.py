#!/usr/bin/env python3
"""Recull les mesures de l'aforador i el pluviòmetre de l'ACA d'una conca (API Sentilo pública, sense clau).

Desa les mesures a Supabase (taula inundacions.observacions, upsert: els duplicats no es tornen a inserir).
Cal SUPABASE_URL i SUPABASE_SECRET_KEY a l'entorn (secrets de GitHub Actions).

Ús:  python3 scripts/aca_aforador.py --config config/ridaura.json [--hores 6]
     (--hores 72 per recuperar històric la primera vegada; l'ACA en guarda uns 3 mesos)"""
import argparse, json, os, sys, time
from datetime import datetime, timedelta, timezone
from pathlib import Path
import requests

API = "https://aplicacions.aca.gencat.cat/sdim2/apirest"
FMT = "%d/%m/%YT%H:%M:%S"
S = requests.Session()
S.headers.update({"Accept": "application/json", "User-Agent": "inundacions/0.1 (github.com/correvents/inundacions)"})


def get(ruta, params, intents=3):
    for i in range(intents):
        try:
            r = S.get(API + ruta, params=params, timeout=60)
            if r.status_code == 200 and r.text.lstrip().startswith(("{", "[")):
                return r.json()
            print(f"  ACA {ruta}: HTTP {r.status_code}", file=sys.stderr)
        except (requests.RequestException, ValueError) as e:
            print(f"  ACA {ruta}: {type(e).__name__}", file=sys.stderr)
        time.sleep(3 * (i + 1))
    return None


def llegeix(sensor, provider, ini, fi):
    """Totes les mesures d'un sensor entre ini i fi. L'API en dona com a màxim 200 per petició: es pagina enrere."""
    out, fins, ant = {}, fi, None
    for _ in range(100):
        j = get(f"/data/{provider}/{sensor}", {"limit": 200, "from": ini.strftime(FMT), "to": fins.strftime(FMT)})
        ob = (j or {}).get("observations", [])
        if not ob:
            break
        for o in ob:
            try:
                out[int(o["time"])] = float(str(o.get("value")).strip())
            except (KeyError, ValueError):
                pass
        antic = min(int(o["time"]) for o in ob)
        if len(ob) < 200 or antic <= ini.timestamp() * 1000 or (ant is not None and antic >= ant):
            break
        ant, fins = antic, datetime.fromtimestamp(antic / 1000 - 1, timezone.utc)
    return out


def desa(conca, var, mesures):
    files = [{"conca": conca, "var": var, "temps": datetime.fromtimestamp(ms / 1000, timezone.utc).isoformat(), "valor": v}
             for ms, v in mesures.items()]
    clau = os.environ["SUPABASE_SECRET_KEY"]
    for i in range(0, len(files), 500):
        r = requests.post(os.environ["SUPABASE_URL"].rstrip("/") + "/rest/v1/observacions", json=files[i:i + 500], timeout=60,
                          headers={"apikey": clau, "Authorization": f"Bearer {clau}", "Content-Profile": "inundacions",
                                   "Prefer": "resolution=ignore-duplicates,return=minimal"})
        if not r.ok:
            sys.exit(f"Supabase: HTTP {r.status_code} {r.text[:300]}")
    return len(files)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--hores", type=float, default=6, help="finestra que es demana a l'API (es solapa amb l'anterior; els duplicats es descarten)")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    fi = datetime.now(timezone.utc); ini = fi - timedelta(hours=a.hores)
    ok = 0
    for sens in cfg["sensors"]:
        m = llegeix(sens["id"], sens["provider"], ini, fi)
        if m:
            ok += 1
            desa(cfg["conca"], sens["var"], m)
        print(f"{sens['var']}: {len(m)} mesures")
    if not ok:
        print("ACA: cap sensor ha retornat dades"); return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
