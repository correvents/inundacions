#!/usr/bin/env python3
"""Baixa les dades que fa servir humitat.py (absorció inicial segons la pluja dels dies anteriors). Corre a GitHub Actions
(des del contenidor no s'hi arriba).

  1. Evapotranspiració de referència diària (FAO-56 Penman-Monteith, reanàlisi ERA5 via Open-Meteo) al punt de l'aforador,
     del 2009 fins fa una setmana → data/<conca>/humitat/et0.csv. Només se'n fa servir la climatologia (mitjana de cada dia
     de l'any): en directe no cal cap dada nova d'evapotranspiració.
  2. Pluja de 30 min de les estacions XEMA de cfg["humitat"]["estacions"] (sèrie llarga; s'afegeix al final del CSV que ja hi ha)
     → data/<conca>/xema/<codi>.csv.
  3. Cabal de l'aforador de cfg["humitat"]["sensor_aca"] dels arxius de l'ACA al bucket privat de Supabase
     (aca/<dia>/<provider>.tar, aca_arxiu.yml) → data/<conca>/humitat/aforador.csv.

Ús:  python3 scripts/humitat_dades.py --config config/granollers.json [--sense-xema] [--sense-aca]"""
import argparse, csv, io, json, os, sys, tarfile, gzip, time
from datetime import date, timedelta
from pathlib import Path
import requests

OM = "https://archive-api.open-meteo.com/v1/archive"
SOCRATA = "https://analisi.transparenciacatalunya.cat/resource/nzvn-apee.json"


def get(url, params=None, headers=None, intents=5, timeout=180):
    for i in range(intents):
        try:
            r = requests.get(url, params=params, headers=headers, timeout=timeout); r.raise_for_status(); return r
        except Exception as e:
            print("reintent", url, e); time.sleep(10 * (i + 1))
    raise SystemExit(f"No respon {url}")


def et0(cfg, d):
    lat, lon = cfg["aforador"]["lat"], cfg["aforador"]["lon"]
    fi = date.today() - timedelta(days=7)
    r = get(OM, {"latitude": lat, "longitude": lon, "start_date": "2009-01-01", "end_date": fi.isoformat(),
                 "daily": "et0_fao_evapotranspiration", "timezone": "UTC"}).json()
    dies, v = r["daily"]["time"], r["daily"]["et0_fao_evapotranspiration"]
    with open(d / "et0.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["dia", "et0_mm"])
        for t, x in zip(dies, v):
            if x is not None: w.writerow([t, round(x, 2)])
    print(f"ET0: {len(dies)} dies, {dies[0]} → {dies[-1]}, mitjana {sum(x for x in v if x) / len(v) * 365:.0f} mm/any")


def xema(cfg, conca):
    dx = Path("data") / conca / "xema"; dx.mkdir(parents=True, exist_ok=True)
    for codi in cfg["humitat"]["estacions"]:
        f = dx / f"{codi}.csv"
        files = list(csv.DictReader(open(f))) if f.exists() else []
        if files and files[0]["temps"] > "2009-01-02":       # sèrie curta (només un episodi): es torna a baixar sencera
            files = []
        des = files[-1]["temps"] if files else "2009-01-01T00:00"
        nous, off = [], 0
        while True:
            r = get(SOCRATA, {"codi_estacio": codi, "codi_variable": "35", "$select": "data_lectura,valor_lectura",
                              "$where": f"data_lectura > '{des}:00.000'", "$order": "data_lectura", "$limit": 50000, "$offset": off}).json()
            nous += r; off += len(r)
            if len(r) < 50000: break
        with open(f, "a" if files else "w", newline="") as fh:
            w = csv.writer(fh)
            if not files: w.writerow(["temps", "mm"])
            for x in nous: w.writerow([x["data_lectura"][:16], x["valor_lectura"]])
        print(f"XEMA {codi}: {len(files)} + {len(nous)} lectures (des de {des})")


def aca(cfg, d):
    url, key = os.environ["SUPABASE_URL"], os.environ["SUPABASE_SECRET_KEY"]
    h = {"Authorization": f"Bearer {key}", "apikey": key}
    prov, sensor = cfg["humitat"]["sensor_aca"]
    objs = requests.post(f"{url}/storage/v1/object/list/arxiu-inundacions", headers=h,
                         json={"prefix": "aca/", "limit": 1000}, timeout=60).json()
    dies = sorted(o["name"] for o in objs)
    mesures = {}
    for dia in dies:
        r = requests.get(f"{url}/storage/v1/object/arxiu-inundacions/aca/{dia}/{prov}.tar", headers=h, timeout=300)
        if r.status_code != 200:
            continue
        with tarfile.open(fileobj=io.BytesIO(r.content)) as t:
            try:
                m = t.extractfile(f"{prov}/{sensor}.csv.gz")
            except KeyError:
                continue
            for x in csv.DictReader(io.TextIOWrapper(gzip.open(m), "utf-8")):
                try: mesures[x["temps_utc"][:16]] = float(x["valor"])
                except ValueError: pass
        print(f"ACA {dia}: {len(mesures)} mesures acumulades")
    with open(d / "aforador.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["temps", "cabal_m3s"])
        for t in sorted(mesures): w.writerow([t, mesures[t]])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--sense-xema", action="store_true"); ap.add_argument("--sense-aca", action="store_true")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    d = Path("data") / cfg["conca"] / "humitat"; d.mkdir(parents=True, exist_ok=True)
    et0(cfg, d)
    if not a.sense_xema: xema(cfg, cfg["conca"])
    if not a.sense_aca and cfg["humitat"].get("sensor_aca"): aca(cfg, d)
    return 0


if __name__ == "__main__":
    sys.exit(main())
