#!/usr/bin/env python3
"""Cabal base dels rius de Catalunya a partir dels aforadors de l'ACA (arxiu de scripts/aca_arxiu.py).

Serveix per a scripts/estrenyiments.py: el LiDAR no entra a l'aigua i, on el riu porta cabal permanent, el MET dona la
làmina i no el fons. Sota la làmina hi passa el cabal base; amb Manning i l'amplada de la làmina se'n treu la fondària.

Per a cada sensor de cabal (m³/s) de l'arxiu: mediana dels cabals mitjans diaris (dies amb ≥ 80 % de dades, ≥ 30 dies),
conca drenada al punt de l'aforador, amb les conques d'aigües amunt (MDT 5 m de la Release mdt5-cat, cel·la de riu ≥ 2 km² més propera a menys de 500 m) i cabal
específic (l/s/km²). La finestra de l'arxiu (uns 100 dies) no coincideix amb els vols LiDAR: és una aproximació.

Sortida: data/catalunya/cabal_base/aforadors.json

Ús:  python3 scripts/cabal_base.py --arxiu /tmp/aca     (carpeta amb sensors.csv i AFORAMENT-EST/<sensor>.csv.gz)"""
import argparse, csv, gzip, json, sys, time
from pathlib import Path
import numpy as np
import requests
from pyproj import Transformer
from shapely.geometry import shape, Point

sys.path.insert(0, str(Path(__file__).parent))
CONQUES = Path("data/catalunya/conques/conques_100km2.geojson")
SORT = Path("data/catalunya/cabal_base/aforadors.json")
MDT_URL = "https://github.com/correvents/inundacions/releases/download/mdt5-cat/mdt5_c{:03d}.tif"


KM2_MIN, RADI = 2.0, 500.0


def riu_proper(D, x, y):
    """Conca (km²) de la cel·la de riu (≥ KM2_MIN) més propera a (x, y) dins RADI m: les coordenades dels aforadors
    poden caure a la vora i no a la llera."""
    tr = D["tr"]; cel = abs(tr.a * tr.e); c, r = ~tr * (x, y); r, c = int(r), int(c); n = int(RADI / abs(tr.a))
    nf, nc = D["acc"].shape
    r0, r1, c0, c1 = max(r - n, 0), min(r + n + 1, nf), max(c - n, 0), min(c + n + 1, nc)
    if r0 >= r1 or c0 >= c1:
        return None
    km2 = D["acc"][r0:r1, c0:c1] * cel / 1e6
    rr, cc = np.mgrid[r0:r1, c0:c1]; d = np.hypot(rr - r, cc - c) * abs(tr.a)
    d = np.where((km2 >= KM2_MIN) & (d <= RADI), d, np.inf)
    if not np.isfinite(d.min()):
        return None
    return float(km2.flat[np.argmin(d)])


def cabals(arxiu):
    """{sensor: (descripció, lon, lat, q50 m³/s, dies)} dels sensors de cabal."""
    out = {}
    for r in csv.DictReader(open(Path(arxiu) / "sensors.csv")):
        u = (r.get("unitat") or "").lower().replace("³", "3")
        if "m3/s" not in u and "m3s" not in u:
            continue
        f = Path(arxiu) / r["provider"] / f"{r['sensor']}.csv.gz"
        if not f.exists():
            continue
        import pandas as pd
        try:
            s = pd.read_csv(f, parse_dates=["temps_utc"]).set_index("temps_utc")["valor"]
        except Exception as e:
            print(f"  {r['sensor']}: {e}"); continue
        s = pd.to_numeric(s, errors="coerce").dropna()
        s = s[s >= 0]
        if s.empty:
            continue
        pas = s.index.to_series().diff().median()
        per_dia = int(round(pd.Timedelta("1D") / pas)) if pd.notna(pas) and pas.total_seconds() > 0 else 0
        dia = s.resample("1D").agg(["mean", "count"])
        dia = dia[dia["count"] >= 0.8 * max(per_dia, 1)]
        if len(dia) < 30:
            continue
        try:
            lat, lon = float(r["lat"]), float(r["lon"])
        except ValueError:
            continue
        out[r["sensor"]] = (r.get("descripcio", ""), lon, lat, float(dia["mean"].median()), len(dia))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arxiu", required=True)
    ap.add_argument("--tmp", default="/tmp/cabal_base")
    a = ap.parse_args()
    import estrenyiments as E
    q = cabals(a.arxiu)
    print(f"{len(q)} aforadors amb ≥ 30 dies de cabal")
    C = [(f["properties"]["id"], shape(f["geometry"])) for f in json.loads(CONQUES.read_text())["features"]]
    per_conca = {}
    for sid, (desc, lon, lat, q50, dies) in q.items():
        cid = next((c for c, g in C if g.contains(Point(lon, lat))), None)
        if cid is not None:
            per_conca.setdefault(cid, []).append(sid)
    tmp = Path(a.tmp); tmp.mkdir(parents=True, exist_ok=True)
    res = []
    for cid, sids in sorted(per_conca.items()):
        f = tmp / f"mdt5_c{cid:03d}.tif"
        if not f.exists():
            for i in range(4):
                try:
                    r = requests.get(MDT_URL.format(cid), timeout=300); r.raise_for_status(); f.write_bytes(r.content); break
                except requests.RequestException as e:
                    print(f"  mdt5 c{cid:03d}: {e}"); time.sleep(10 * (i + 1))
            else:
                continue
        import rasterio
        t = Transformer.from_crs(4326, rasterio.open(f).crs, always_xy=True).transform
        D = E.drenatge(str(f), None, [(*t(lo, la), k) for lo, la, k in E.entrades(cid)])
        for sid in sids:
            desc, lon, lat, q50, dies = q[sid]
            km2 = riu_proper(D, *t(lon, lat))
            if km2 is None:
                print(f"  c{cid:03d} {sid} {desc[:40]}: cap riu ≥ {KM2_MIN} km² a {RADI} m"); continue
            dr = {"km2": km2}
            res.append({"sensor": sid, "nom": desc, "lon": round(lon, 5), "lat": round(lat, 5), "conca": cid,
                        "km2": round(dr["km2"], 1), "q50_m3s": round(q50, 3), "l_s_km2": round(1000 * q50 / dr["km2"], 2),
                        "dies": dies})
            print(f"  c{cid:03d} {sid} {desc[:40]:40s} {dr['km2']:8.1f} km²  Q50 {q50:7.3f} m³/s  {1000 * q50 / dr['km2']:6.2f} l/s/km²")
        f.unlink(missing_ok=True)
    SORT.parent.mkdir(parents=True, exist_ok=True)
    SORT.write_text(json.dumps({"font": "ACA, aforadors (API Sentilo, arxiu de scripts/aca_arxiu.py)",
                                "metode": "mediana dels cabals mitjans diaris; conca drenada amb l'MDT 5 m",
                                "aforadors": res}, ensure_ascii=False, indent=0))
    print(f"{SORT}: {len(res)} aforadors")
    return 0 if res else 1


if __name__ == "__main__":
    sys.exit(main())
