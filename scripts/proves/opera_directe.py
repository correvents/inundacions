#!/usr/bin/env python3
"""Prova 2: compostos OPERA en directe al bucket de reserva (ECMWF) i a l'API MeteoGate."""
import re, json, requests, io
from datetime import datetime, timezone, timedelta
ara = datetime.now(timezone.utc)
BK = "https://object-store.os-api.cci1.ecmwf.int/openradar-24h-backup"
print("Ara:", ara)
for h in (0, 1, 3):
    t = ara - timedelta(hours=h)
    r = requests.get(BK, params={"list-type": 2, "prefix": t.strftime("%Y/%m/%d/OPERA/COMP/OPERA@%Y%m%dT%H")}, timeout=60)
    ks = re.findall(r"<Key>([^<]+)</Key>", r.text)
    print(f"backup {t:%Y-%m-%dT%H}: HTTP {r.status_code}, {len(ks)} claus", ks[-4:] or r.text[:300])
r = requests.get(BK, params={"list-type": 2, "prefix": ara.strftime("%Y/%m/%d/"), "delimiter": "/"}, timeout=60)
print("prefixos avui:", re.findall(r"<Prefix>([^<]+)</Prefix>", r.text)[:30])
r = requests.get(BK, params={"list-type": 2, "max-keys": 10}, timeout=60); print("arrel:", re.findall(r"<Key>([^<]+)</Key>", r.text))
# Primer i últim dia disponibles al backup
for d in range(0, 12):
    t = ara - timedelta(days=d)
    n = len(re.findall(r"<Key>", requests.get(BK, params={"list-type": 2, "prefix": t.strftime("%Y/%m/%d/OPERA/COMP/")}, timeout=60).text))
    print(f"backup COMP {t:%Y-%m-%d}: {n} claus (màx 1000 per pàgina)")
# Prova de lectura per rangs d'un DBZH recent
ks = [k for h in (0, 1) for k in re.findall(r"<Key>([^<]+)</Key>", requests.get(BK, params={"list-type": 2, "prefix": (ara - timedelta(hours=h)).strftime("%Y/%m/%d/OPERA/COMP/OPERA@%Y%m%dT%H")}, timeout=60).text) if "DBZH" in k]
if ks:
    u = BK + "/" + ks[-1].replace("@", "%40")
    hd = requests.head(u, timeout=60); print("HEAD", ks[-1], hd.status_code, hd.headers.get("Content-Length"), hd.headers.get("Accept-Ranges"))
    print("Range ->", requests.get(u, headers={"Range": "bytes=0-99"}, timeout=60).status_code)
A = "https://api.meteogate.eu/eu-eumetnet-weather-radar/collections/observations"
for p in ({"bbox": "0.15,40.5,3.35,42.9", "standard_names": "DBZH"}, {"bbox": "0.15,40.5,3.35,42.9"}):
    r = requests.get(A + "/items", params=p, timeout=60)
    print("\nEDR items", p, r.status_code, r.text[:1500])
for plat in ("0-20010-0-OPERA", "OPERA"):
    r = requests.get(A + "/items", params={"platform": plat}, timeout=60); print("\nEDR platform", plat, r.status_code, r.text[:1200])
loc = requests.get(A + "/locations", timeout=60).json()
op = [f for f in loc["features"] if "OPERA" in json.dumps(f["properties"]).upper() or "COMP" in json.dumps(f["properties"]).upper()]
print("\nLocations amb OPERA/COMP:", json.dumps(op[:3])[:1500])
