#!/usr/bin/env python3
"""Prova temporal: què es pot baixar de pluja de Catalunya per al 2026 (XEMA, OPERA, Open-Meteo). No escriu res a producció."""
import re, json, requests
from collections import defaultdict
S = requests.Session()
def get(url, **kw):
    try:
        r = S.get(url, timeout=90, **kw); return r
    except Exception as e:
        print("  ERROR", url[:90], e); return None

print("\n### XEMA (Socrata nzvn-apee, variable 35 = pluja 30 min), 2026 per mes")
r = get("https://analisi.transparenciacatalunya.cat/resource/nzvn-apee.json", params={
    "$select": "date_trunc_ym(data_lectura) as mes, count(*) as n, count(distinct codi_estacio) as est, max(valor_lectura) as max30",
    "$where": "codi_variable='35' AND data_lectura >= '2026-01-01T00:00:00'", "$group": "mes", "$order": "mes"})
if r is not None: print(" HTTP", r.status_code); [print(" ", x.get("mes","")[:7], "lectures", x.get("n"), "estacions", x.get("est"), "max30", x.get("max30")) for x in (r.json() if r.ok else [])] or print(r.text[:300])
r = get("https://analisi.transparenciacatalunya.cat/resource/nzvn-apee.json", params={"$select": "min(data_lectura) as primer, max(data_lectura) as darrer", "$where": "codi_variable='35'"})
if r is not None: print(" rang total var 35:", r.status_code, r.text[:200])

print("\n### OPERA arxiu (openradar-archive), dia 15 de cada mes de 2026")
for base in ["https://s3.waw3-1.cloudferro.com/openradar-archive"]:
    for m in range(1, 12):
        pre = f"2026/{m:02d}/15/OPERA/COMP/"
        r = get(base, params={"list-type": 2, "prefix": pre, "max-keys": 1000})
        if r is None: continue
        ks = re.findall(r"<Key>([^<]*)</Key>", r.text)
        tipus = sorted({re.sub(r"\d{8}T\d{4}", "*", k.split("/")[-1]) for k in ks})
        print(f"  2026-{m:02d}-15 HTTP {r.status_code} fitxers {len(ks)} tipus {tipus[:6]}")
for d in ["2026/11/10", "2026/11/11", "2026/11/12", "2026/11/17", "2026/11/18"]:
    r = get("https://s3.waw3-1.cloudferro.com/openradar-archive", params={"list-type": 2, "prefix": d + "/OPERA/COMP/", "max-keys": 1000})
    if r is not None:
        ks = re.findall(r"<Key>([^<]*)</Key>", r.text); print(f"  {d} fitxers {len(ks)} darrer {ks[-1] if ks else '-'}")
print("\n### OPERA 24 h i domini de reserva")
r = get("https://s3.waw3-1.cloudferro.com/openradar-24h", params={"list-type": 2, "max-keys": 5})
if r is not None: print(" 24h HTTP", r.status_code, re.findall(r"<Key>([^<]*)</Key>", r.text)[:5], r.text[:200] if not r.ok else "")
for u in ["https://radar-backup.meteogate.eu/", "https://api.meteogate.eu/eu-eumetnet-weather-radar"]:
    r = get(u); print(" ", u, r.status_code if r is not None else "-", (r.text[:300].replace("\n", " ") if r is not None else ""))

pts = {"Platja d'Aro": (41.82, 3.06), "Granollers": (41.61, 2.29), "Cabrera de Mar": (41.53, 2.40), "Lleida": (41.62, 0.62), "Tortosa": (40.81, 0.52), "Puigcerdà": (42.43, 1.93)}
lat = ",".join(str(p[0]) for p in pts.values()); lon = ",".join(str(p[1]) for p in pts.values())
def mensual(js):
    out = []
    for j in (js if isinstance(js, list) else [js]):
        d = j.get("daily", {}); acc = defaultdict(float); n = defaultdict(int)
        for t, v in zip(d.get("time", []), d.get("precipitation_sum", [])):
            if v is not None: acc[t[:7]] += v; n[t[:7]] += 1
        out.append((j.get("latitude"), j.get("longitude"), {k: round(acc[k]) for k in sorted(acc)}, sum(n.values())))
    return out
for nom, url, extra in [
    ("ERA5-Land (archive)", "https://archive-api.open-meteo.com/v1/archive", {"models": "era5_land"}),
    ("AROME France HD previsions arxivades", "https://historical-forecast-api.open-meteo.com/v1/forecast", {"models": "meteofrance_arome_france_hd"}),
    ("ICON-D2 previsions arxivades", "https://historical-forecast-api.open-meteo.com/v1/forecast", {"models": "icon_d2"}),
]:
    print(f"\n### {nom}: pluja mensual 2026 (mm) a 6 punts")
    r = get(url, params=dict(latitude=lat, longitude=lon, start_date="2026-01-01", end_date="2026-11-17", daily="precipitation_sum", timezone="UTC", **extra))
    if r is None: continue
    print(" HTTP", r.status_code)
    if not r.ok: print(" ", r.text[:300]); continue
    for (nm, (la, lo)), (gla, glo, mm, n) in zip(pts.items(), mensual(r.json())):
        print(f"  {nm} graella {gla},{glo} dies amb dada {n}: {mm}")
print("\n### AROME HD 15 min (previsió ara)")
for model in ["meteofrance_arome_france_hd", "meteofrance_seamless"]:
    r = get("https://api.open-meteo.com/v1/forecast", params=dict(latitude=lat, longitude=lon, minutely_15="precipitation", models=model, forecast_hours=6))
    if r is not None and r.ok:
        for nm, j in zip(pts, r.json()):
            v = j.get("minutely_15", {}).get("precipitation", [])
            print(f"  {model} {nm}: {len(v)} valors, nuls {sum(x is None for x in v)}")
    elif r is not None: print(" ", model, r.status_code, r.text[:200])
print("\n### AROME HD 15 min arxivat (4/10/2026, Platja d'Aro)")
r = get("https://historical-forecast-api.open-meteo.com/v1/forecast", params=dict(latitude=41.82, longitude=3.06, start_date="2026-10-03", end_date="2026-10-04", minutely_15="precipitation", hourly="precipitation", models="meteofrance_arome_france_hd"))
if r is not None:
    print(" HTTP", r.status_code)
    if r.ok:
        j = r.json(); m = j.get("minutely_15", {}).get("precipitation", []); h = j.get("hourly", {}).get("precipitation", [])
        print("  15min valors", len(m), "no nuls", sum(x is not None for x in m), "| horari total", round(sum(x or 0 for x in h)), "mm, max", max([x or 0 for x in h] or [0]), "mm/h")
print("\n### Altres (només si responen)")
for u in ["https://api.meteo.cat/xema/v1/estacions/metadades", "https://opendata.aemet.es/opendata/api/red/radar/nacional",
          "https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGHHE.07/2026/", "https://portail-api.meteofrance.fr/"]:
    r = get(u); print(" ", u, r.status_code if r is not None else "-")
