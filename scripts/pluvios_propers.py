"""Pluviòmetres (ACA i XEMA) a prop d'una conca, i la seva pluja horària recent, per corregir el radar en directe.

La llista es desa a /tmp per no tornar a demanar els catàlegs a cada pas de 10 min dins de la mateixa execució."""
import json, math, re, time
from datetime import datetime, timedelta, timezone
from pathlib import Path
import numpy as np, requests
import aca_arxiu as aca

SOC = "https://analisi.transparenciacatalunya.cat/resource/"


def _km(la, lo, la0, lo0):
    return math.hypot((la - la0) * 110.6, (lo - lo0) * 111.3 * math.cos(math.radians(la0)))


def centre(cfg):
    """[lat, lon] del centre de la conca: el de la configuració o, si no n'hi ha, el centroide de la primera conca de
    cfg["geometria"] (abans es feia servir Platja d'Aro per defecte, i una conca nova sense "centre" agafava les estacions d'allà)."""
    if "centre" in cfg:
        return cfg["centre"]
    from shapely.geometry import shape
    fc = json.loads(Path(cfg["geometria"]).read_text())
    c = shape(next(f for f in fc["features"] if f["geometry"]["type"].endswith("Polygon"))["geometry"]).centroid
    return [round(c.y, 4), round(c.x, 4)]


def propers(cfg, radi_km=40):
    cache = Path(f"/tmp/pluvios_{cfg['conca']}_{int(radi_km)}.json")
    if cache.exists() and time.time() - cache.stat().st_mtime < 6 * 3600:
        return json.loads(cache.read_text())
    la0, lo0 = centre(cfg)
    out = []
    try:
        cat = aca.get(aca.sessio(), "/catalog") or {}
        for p in cat.get("providers", []):
            for se in p.get("sensors", []):
                txt = " ".join(str(se.get(k, "")) for k in ("type", "sensorType", "unit", "description", "componentType", "componentDesc"))
                loc = str(se.get("location") or "").split()
                if re.search(r"plu|precip|rain", txt, re.I) and str(se.get("unit")) == "mm/h" and len(loc) > 1:
                    la, lo = float(loc[0]), float(loc[1])
                    if _km(la, lo, la0, lo0) <= radi_km:
                        out.append({"id": "ACA-" + se["sensor"], "font": "ACA", "provider": p["provider"], "sensor": se["sensor"], "lat": la, "lon": lo})
    except Exception as e:
        print("Pluviòmetres ACA: catàleg no disponible", type(e).__name__)
    try:
        r = requests.get(SOC + "yqwd-vj5e.json", params={"$limit": 5000}, timeout=60).json()
        for e in r:
            if e.get("data_fi"): continue
            la, lo = float(e["latitud"]), float(e["longitud"])
            if _km(la, lo, la0, lo0) <= radi_km:
                out.append({"id": "XEMA-" + e["codi_estacio"], "font": "XEMA", "codi": e["codi_estacio"], "lat": la, "lon": lo})
    except Exception as e:
        print("Pluviòmetres XEMA: llista no disponible", type(e).__name__)
    if out:
        cache.write_text(json.dumps(out))
    return out


def horaria(pluvios, ini, H):
    """{id: array de H valors (mm a l'hora que acaba a ini+(h+1)), NaN = sense dada}. ACA: mitjana de mm/h; XEMA: suma de 30 min."""
    fi = ini + timedelta(hours=H)
    out = {}
    s = aca.sessio()
    for p in [x for x in pluvios if x["font"] == "ACA"]:
        try:
            m = aca.llegeix(s, p["provider"], p["sensor"], ini, fi)
        except Exception:
            continue
        sm, n = np.zeros(H), np.zeros(H)
        for ms, v in m.items():
            try: v = float(v)
            except ValueError: continue
            h = math.ceil((ms / 1000 - ini.timestamp()) / 3600) - 1
            if 0 <= h < H and 0 <= v < 300: sm[h] += v; n[h] += 1
        tip = np.median(n[n > 0]) if (n > 0).any() else 0
        out[p["id"]] = np.where(n >= max(1, 0.75 * tip), sm / np.maximum(n, 1), np.nan)
    codis = [p["codi"] for p in pluvios if p["font"] == "XEMA"]
    if codis:
        try:
            r = requests.get(SOC + "nzvn-apee.json", timeout=90, params={"codi_variable": "35", "$limit": 50000,
                "$select": "codi_estacio,data_lectura,valor_lectura",
                "$where": f"data_lectura >= '{ini:%Y-%m-%dT%H:%M:%S}' AND codi_estacio in({','.join(repr(c) for c in codis)})"}).json()
            acc = {}
            for x in r:
                t = datetime.fromisoformat(x["data_lectura"][:19]).replace(tzinfo=timezone.utc) + timedelta(minutes=30)
                h = math.ceil((t - ini).total_seconds() / 3600) - 1
                if 0 <= h < H:
                    a = acc.setdefault(x["codi_estacio"], np.zeros((H, 2))); a[h] += (float(x["valor_lectura"]), 1)
            for c, a in acc.items():
                out["XEMA-" + c] = np.where(a[:, 1] == 2, a[:, 0], np.nan)
        except Exception as e:
            print("XEMA: sense dades recents", type(e).__name__)
    return out
