#!/usr/bin/env python3
"""Pluja prevista sobre una conca per a les hores següents (opció --pluja-prevista de previsio.py).

Dues fonts públiques i gratuïtes:
  1. Nowcast del radar OPERA (0-2 h). Amb els últims compostos DBZH (cada 5 min) d'un quadrat de ~240 km al voltant de
     la conca s'estima el moviment de la pluja (correlació de fase entre parells d'imatges de 15 min, mediana dels
     vectors) i s'extrapola l'última imatge (persistència lagrangiana). Per a cada minut de previsió es mostreja la
     pluja als punts de la conca desplaçats cap enrere. La intensitat és la mateixa relació Z-R que radar_conca.py.
  2. AROME France HD de Météo-France (1,3 km, pluja cada 15 min) via Open-Meteo, fins a 6 h. En directe, l'última
     execució publicada. Per a episodis passats (--ara a previsio.py), l'execució que hi hauria hagut disponible a
     aquella hora (Single Runs API d'Open-Meteo, execució de fa LATENCIA_H hores o més).
Combinació ('tot'): fins a 60 min el nowcast; de 60 a 120 min, pes lineal del nowcast a AROME; després, AROME.
Si una font falla, es fa servir l'altra; si fallen les dues, la pluja prevista és 0 (com l'escenari "no plou més").

Ús (prova):  python3 scripts/pluja_prevista.py --config config/ridaura.json [--ara 2026-10-04T06:00] [--font tot]"""
import argparse, json, re, sys, tempfile, urllib.parse
from datetime import datetime, timedelta, timezone
from pathlib import Path
import numpy as np, h5py, requests
from pyproj import CRS, Transformer
from shapely.geometry import shape, Point
from shapely.prepared import prep

DT = 5                                   # minuts, com el model de cabal
S3_24H = "https://object-store.os-api.cci1.ecmwf.int/openradar-24h-backup"
S3_ARXIU = "https://s3.waw3-1.cloudferro.com/openradar-archive"
FI_ARXIU = datetime(2026, 11, 11, 20, 0, tzinfo=timezone.utc)   # l'arxiu s'atura aquí; la resta és a la reserva (24 h)
MITJA_KM = 120                           # mig costat del quadrat de radar al voltant de la conca
NOWCAST_MIN, FUSIO_MIN = 60, 120         # nowcast sol fins a 60 min; fusió lineal fins a 120 min
OM_ACTUAL = "https://api.open-meteo.com/v1/forecast"
OM_RUNS = "https://single-runs-api.open-meteo.com/v1/forecast"
MODEL_NWP = "meteofrance_arome_france_hd"
LATENCIA_H = 3                           # hores entre l'inici d'una execució d'AROME i que és disponible (prudent)


# ---------- punts de la conca ----------
def punts_conca(cfg, pas=0.0025, subconca="aforador"):
    fc = json.loads(Path(cfg["geometria"]).read_text())
    g = next((shape(f["geometry"]) for f in fc["features"] if f["properties"].get("subconca") == subconca), None)
    g = g or shape(fc["features"][0]["geometry"])
    w, s, e, n = g.bounds
    pg = prep(g)
    lon, lat = np.meshgrid(np.arange(w, e, pas) + pas / 2, np.arange(s, n, pas) + pas / 2)
    dins = np.array([pg.contains(Point(x, y)) for x, y in zip(lon.ravel(), lat.ravel())])
    return lon.ravel()[dins], lat.ravel()[dins]


# ---------- nowcast del radar ----------
def _llista(s3, ini, fi):
    out, h = {}, ini.replace(minute=0)
    while h <= fi:
        xml = requests.get(s3, params={"list-type": 2, "prefix": h.strftime("%Y/%m/%d/OPERA/COMP/OPERA@%Y%m%dT%H")}, timeout=30).text
        for k in re.findall(r"<Key>([^<]+)</Key>", xml):
            m = re.search(r"OPERA@(\d{8}T\d{4})@0@DBZH(_QIND)?\.h5$", k)
            if m:
                t = datetime.strptime(m.group(1), "%Y%m%dT%H%M").replace(tzinfo=timezone.utc)
                if ini <= t <= fi:
                    out[t] = s3 + "/" + urllib.parse.quote(k, safe="/@")
        h += timedelta(hours=1)
    return out


def _retall(h5, lon0, lat0):
    """Quadrat de dBZ al voltant de (lon0, lat0) i la funció que passa lon/lat a fila/columna fraccionària del quadrat."""
    with h5py.File(h5, "r") as h:
        grup, wt = None, None
        for ds_nom in [k for k in h if k.startswith("dataset")]:
            for cami in (f"{ds_nom}/data1/what", f"{ds_nom}/what"):
                if cami in h:
                    q = h[cami].attrs.get("quantity")
                    if q is not None and (q.decode() if isinstance(q, bytes) else str(q)) == "DBZH":
                        grup, wt = f"{ds_nom}/data1", h[cami].attrs
                        break
            if grup:
                break
        if grup is None:
            grup, wt = "dataset1/data1", h["dataset1/data1/what"].attrs
        wh = h["where"].attrs
        gain, off, nodata, undet = (float(wt[k]) for k in ("gain", "offset", "nodata", "undetect"))
        pdef = wh["projdef"].decode() if isinstance(wh["projdef"], bytes) else str(wh["projdef"])
        tr = Transformer.from_crs("EPSG:4326", CRS.from_proj4(pdef), always_xy=True)
        xul, yul = tr.transform(float(wh["UL_lon"]), float(wh["UL_lat"]))
        xs, ys = float(wh["xscale"]), float(wh["yscale"])
        x0, y0 = tr.transform(lon0, lat0)
        c0, r0 = int((x0 - xul) / xs), int((yul - y0) / ys)
        m = int(MITJA_KM * 1000 / xs)
        ds = h[grup + "/data"]
        rr, cc = max(r0 - m, 0), max(c0 - m, 0)
        v = ds[rr: r0 + m, cc: c0 + m].astype(np.float64)
    z = v * gain + off
    z[v == undet] = -np.inf
    z[v == nodata] = np.nan
    def idx(lon, lat):
        X, Y = tr.transform(lon, lat)
        return (yul - Y) / ys - rr - 0.5, (X - xul) / xs - cc - 0.5
    return z, idx, xs


def _desplacament(a, b):
    """Desplaçament (files, columnes) de la imatge a a la b per correlació de fase."""
    fa, fb = np.fft.fft2(a), np.fft.fft2(b)
    r = fb * np.conj(fa)
    c = np.fft.ifft2(r / np.maximum(np.abs(r), 1e-9)).real
    i, j = np.unravel_index(np.argmax(c), c.shape)
    n, m = c.shape
    di, dj = (i if i <= n // 2 else i - n), (j if j <= m // 2 else j - m)
    return float(di), float(dj), float(c.max())


def nowcast(cfg, ara, minuts, zr=None, cache=None):
    """mm per pas de DT minuts sobre la conca, de ara a ara+minuts (o None si no hi ha radar). També torna info."""
    from radar_conca import mmh
    lon, lat = punts_conca(cfg)
    lon0, lat0 = float(lon.mean()), float(lat.mean())
    s3 = S3_ARXIU if ara <= FI_ARXIU else S3_24H
    cache = {} if cache is None else cache
    disp = _llista(s3, ara - timedelta(minutes=30), ara)
    if s3 == S3_ARXIU and ara > FI_ARXIU - timedelta(minutes=30):
        disp.update({t: u for t, u in _llista(S3_24H, ara - timedelta(minutes=30), ara).items() if t not in disp})
    if not disp:
        return None, {"error": "sense compostos"}
    ts = sorted(disp)[-4:]
    ts = [ts[0], ts[-1]] if len(ts) > 1 else ts    # l'últim i el de ~15 min abans (només calen aquests dos)
    camps = []
    with tempfile.TemporaryDirectory() as tmp:
        for t in ts:
            if t not in cache:
                f = Path(tmp) / "c.h5"
                r = requests.get(disp[t], timeout=120); r.raise_for_status(); f.write_bytes(r.content)
                cache[t] = _retall(f, lon0, lat0)
            camps.append(cache[t])
    z, idx, xs = camps[-1]
    rr = mmh(z, zr)                                # mm/h; NaN fora de cobertura
    # moviment: imatge de fa 15 min (o la més antiga) contra l'última, en log(1+R)
    vi = vj = 0.0; qual = 0.0
    if len(camps) >= 2 and all(c[0].shape == z.shape for c in camps):
        a = np.log1p(np.nan_to_num(mmh(camps[0][0], zr))); b = np.log1p(np.nan_to_num(rr))
        dt_min = (ts[-1] - ts[0]).total_seconds() / 60
        if a.max() > 0 and b.max() > 0 and dt_min > 0:
            di, dj, qual = _desplacament(a, b)
            if np.hypot(di, dj) * xs / 1000 / dt_min * 60 <= 120:       # rebutja > 120 km/h (soroll)
                vi, vj = di / dt_min, dj / dt_min                       # píxels per minut
    fi_, ci_ = idx(lon, lat)
    n = minuts // DT
    out = np.zeros(n)
    rr0 = np.nan_to_num(rr)
    from scipy.ndimage import map_coordinates
    for k in range(n):
        tmig = (k + 0.5) * DT                                              # minut central del pas
        f_ = fi_ - vi * tmig; c_ = ci_ - vj * tmig
        v = map_coordinates(rr0, [f_, c_], order=1, mode="constant", cval=0.0)
        out[k] = float(v.mean()) * DT / 60
    vel = np.hypot(vi, vj) * xs / 1000 * 60
    return out, {"compostos": [t.strftime("%H:%M") for t in ts], "km_h": round(vel, 1),
                 "direccio_graus": round(float(np.degrees(np.arctan2(vj, -vi)) % 360), 0), "qualitat": round(qual, 3)}


# ---------- AROME via Open-Meteo ----------
def _punts_mostra(cfg):
    """Punts repartits per la conca cada 0,015° (~1,5 km; la malla d'AROME HD és de 0,01°)."""
    lon, lat = punts_conca(cfg, pas=0.015)
    return (lon, lat) if len(lon) else punts_conca(cfg, pas=0.005)


def nwp(cfg, ara, minuts, directe=True):
    """mm per pas de DT minuts (mitjana de punts de la conca) d'AROME HD. directe=False: execució disponible a 'ara'."""
    lon, lat = _punts_mostra(cfg)
    p = {"latitude": ",".join(f"{v:.4f}" for v in lat), "longitude": ",".join(f"{v:.4f}" for v in lon),
         "minutely_15": "precipitation", "models": MODEL_NWP, "timezone": "GMT"}
    run = None
    if directe:
        p.update({"start_minutely_15": (ara - timedelta(minutes=15)).strftime("%Y-%m-%dT%H:%M"),
                  "end_minutely_15": (ara + timedelta(minutes=minuts + 15)).strftime("%Y-%m-%dT%H:%M")})
        r = requests.get(OM_ACTUAL, params=p, timeout=60)
    else:
        base = ara - timedelta(hours=LATENCIA_H)
        r = None
        for k in range(0, 7):                   # l'última execució abans de 'base' que existeixi
            cand = base.replace(minute=0, second=0) - timedelta(hours=k)
            q = {**p, "run": cand.strftime("%Y-%m-%dT%H:%M"), "forecast_hours": int(LATENCIA_H + k + minuts / 60 + 2)}
            r = requests.get(OM_RUNS, params=q, timeout=60)
            if r.ok and "minutely_15" in (r.json()[0] if isinstance(r.json(), list) else r.json()):
                run = cand
                break
        if run is None:
            return None, {"error": f"sense execució ({r.status_code if r is not None else '-'}: {r.text[:150] if r is not None else ''})"}
    r.raise_for_status()
    js = r.json(); js = js if isinstance(js, list) else [js]
    t = [datetime.fromisoformat(x).replace(tzinfo=timezone.utc) for x in js[0]["minutely_15"]["time"]]
    v = np.nanmean(np.array([[np.nan if y is None else y for y in j["minutely_15"]["precipitation"]] for j in js], float), axis=0)
    # Open-Meteo: el valor de les hh:mm és la pluja dels 15 min anteriors → 1/3 a cada pas de 5 min d'aquell interval
    out = np.zeros(minuts // DT)
    for ti, vi in zip(t, v):
        if not np.isfinite(vi):
            continue
        for s in range(3):
            fi = ti - timedelta(minutes=DT * s)                 # final del pas de 5 min
            k = int((fi - ara).total_seconds() // 60 // DT) - 1
            if 0 <= k < len(out):
                out[k] += vi / 3
    return out, {"model": MODEL_NWP, "execucio": run.strftime("%Y-%m-%d %H:%M") if run else "última"}


def prevista(cfg, ara, minuts, font="tot", directe=True, zr=None, cache=None):
    """Pluja prevista (mm per pas de DT min) de ara a ara+minuts i un resum de les fonts fetes servir."""
    n = minuts // DT
    info, nc, am = {}, None, None
    if font in ("nowcast", "tot"):
        try:
            nc, info["nowcast"] = nowcast(cfg, ara, min(minuts, FUSIO_MIN), zr, cache)
        except Exception as e:
            info["nowcast"] = {"error": f"{type(e).__name__} {e}"[:200]}
    if font in ("arome", "tot"):
        try:
            am, info["arome"] = nwp(cfg, ara, minuts, directe)
        except Exception as e:
            info["arome"] = {"error": f"{type(e).__name__} {e}"[:200]}
    out = np.zeros(n)
    t = (np.arange(n) + 0.5) * DT
    if nc is not None:
        nc = np.r_[nc, np.zeros(n - len(nc))]
    if nc is not None and am is not None:
        w = np.clip((FUSIO_MIN - t) / (FUSIO_MIN - NOWCAST_MIN), 0, 1)       # 1 fins a 60 min, 0 des de 120 min
        out = w * nc + (1 - w) * am
    elif nc is not None:
        out = np.where(t <= FUSIO_MIN, nc, 0.0)
    elif am is not None:
        out = am
    info["mm_1h"], info["mm_3h"], info["mm_6h"] = (round(float(out[: h * 60 // DT].sum()), 1) for h in (1, 3, 6))
    return out, info


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--ara")
    ap.add_argument("--font", default="tot", choices=["nowcast", "arome", "tot"])
    ap.add_argument("--hores", type=int, default=6)
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    ara = datetime.fromisoformat(a.ara).replace(tzinfo=timezone.utc) if a.ara else datetime.now(timezone.utc)
    ara = ara.replace(second=0, microsecond=0); ara -= timedelta(minutes=ara.minute % DT)
    out, info = prevista(cfg, ara, a.hores * 60, a.font, directe=a.ara is None)
    print(json.dumps(info, ensure_ascii=False))
    print("mm per hora:", [round(float(x), 1) for x in out.reshape(-1, 60 // DT).sum(1)])


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).parent))
    sys.exit(main())
