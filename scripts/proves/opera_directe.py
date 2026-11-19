#!/usr/bin/env python3
"""Prova: on és ara el radar OPERA en directe (servidor de reserva / MeteoGate EDR)?"""
import re, requests
from datetime import datetime, timezone, timedelta
ara = datetime.now(timezone.utc)
def mostra(u, **kw):
    try:
        r = requests.get(u, timeout=40, **kw)
        print(f"\n== {r.status_code} {r.url} [{r.headers.get('Content-Type')}] {len(r.content)} B\n{r.text[:1500]}")
        return r
    except Exception as e:
        print(f"\n== ERR {u}: {e}")
doc = mostra("https://eumetnet.github.io/openradardata-documentation/")
for p in ("1-ORD-API-overview/", "1-ORD-API-overview.html"):
    r = mostra("https://eumetnet.github.io/openradardata-documentation/" + p)
    if r is not None and r.ok:
        t = re.sub(r"<[^>]+>", " ", r.text); t = re.sub(r"\s+", " ", t)
        print("\n---- TEXT DOC ----\n", t[:8000]); break
for u in re.findall(r"https?://[^\s\"'<>]*(?:meteogate|cloudferro|backup)[^\s\"'<>]*", (r.text if r is not None else "")):
    print("URL a la doc:", u)
B = "https://radar-backup.meteogate.eu"
mostra(B + "/")
for b in ("openradar-24h", "openradar-archive", ""):
    for h in (0, 1):
        t = ara - timedelta(hours=h)
        mostra(f"{B}/{b}", params={"list-type": 2, "prefix": t.strftime("%Y/%m/%d/OPERA/COMP/OPERA@%Y%m%dT%H")})
A = "https://api.meteogate.eu/eu-eumetnet-weather-radar"
mostra(A); mostra(A + "/collections"); mostra(A + "/collections?f=json")
for base in (A, B):
    mostra(base + "/collections/observations/locations")
    mostra(base + "/collections/observations/items", params={"limit": 3})
