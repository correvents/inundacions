#!/usr/bin/env python3
"""Compara la simulació 2D de cada episodi amb els carrers que la premsa diu que es van inundar.

Per a cada carrer amb nom de la plana (model/carrers.geojson) i cada episodi simulat (sim2d/<id>/calat_max_cm.u2):
  calat_cm   calat màxim típic al carrer (percentil 75 de les cel·les a menys de 15 m de l'eix; el percentil evita que
             un pont sobre el riu compti com a carrer inundat)
  hora_utc   primer moment en què una quarta part del carrer supera 15 cm
Escriu data/<conca>/model/sim2d_carrers.csv i treu, per episodi, on queden els carrers citats a la premsa dins la
llista de tots els carrers ordenats per calat (rang percentil: 100 = el més inundat).

Ús:  python3 scripts/valida_sim2d.py --config config/ridaura.json"""
import argparse, json, sys, csv
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np, shapely
from shapely.geometry import shape
from pyproj import Transformer

PREMSA = {"2026-10-04": "2026-10-04", "2026-05-18": "2026-05-18", "2026-01-20": "2026-01-20", "2024-09-23": "2024-09-23",
          "2020-01-21": "2020-01-21", "2019-10-22": "2019-10-22", "2018-11-18": "2018-11-18"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    d = Path("data") / cfg["conca"]
    carrers = json.loads((d / "model" / "carrers.geojson").read_text())["features"]
    premsa = json.loads((d / "episodis" / "premsa.geojson").read_text())["features"]
    files, resum = [], {}
    for idx in sorted((d / "sim2d").glob("*/index.json")):
        ix = json.loads(idx.read_text()); g = ix.get("graella_utm")
        if not g: continue
        ep = idx.parent.name
        nf, nc, pas = g["nf"], g["nc"], ix["pas_m"]
        h = np.fromfile(idx.parent / "calat_max_cm.u2", "<u2").reshape(nf, nc).astype(float)
        tm = np.fromfile(idx.parent / "minut_15cm.u2", "<u2").reshape(nf, nc).astype(float); tm[tm == 65535] = np.inf
        t0 = datetime.fromisoformat(ix["inici_utc"])
        X, Y = np.meshgrid(g["x0"] + (np.arange(nc) + .5) * pas, g["y0"] - (np.arange(nf) + .5) * pas)
        to_utm = Transformer.from_crs(4326, g["crs"], always_xy=True)
        utm = lambda geom: shapely.transform(geom, lambda c: np.column_stack(to_utm.transform(c[:, 0], c[:, 1])))
        val = {}
        for f in carrers:
            gg = utm(shape(f["geometry"])).buffer(15)
            x0, y0, x1, y1 = gg.bounds
            c0, c1 = int((x0 - g["x0"]) // pas), int((x1 - g["x0"]) // pas) + 1
            r0, r1 = int((g["y0"] - y1) // pas), int((g["y0"] - y0) // pas) + 1
            sub = (slice(max(r0, 0), r1), slice(max(c0, 0), c1))
            m = shapely.contains_xy(gg, X[sub], Y[sub])
            if not m.any(): continue
            hh, tt = h[sub][m], tm[sub][m]
            cal = float(np.percentile(hh, 75)); tq = float(np.percentile(tt, 25))
            nom = f["properties"]["nom"]; val[nom] = cal
            files.append({"episodi": ep, "carrer": nom, "calat_cm": round(cal),
                          "hora_utc": (t0 + timedelta(minutes=tq)).strftime("%Y-%m-%dT%H:%M") if np.isfinite(tq) else ""})
        noms = sorted(val, key=val.get); rang = {n: 100 * (k + 1) / len(noms) for k, n in enumerate(noms)}
        citats = sorted({p["properties"]["nom"] for p in premsa if p["properties"]["episodi"] == PREMSA.get(ep)} & set(val))
        resum[ep] = {"carrers": len(val), "amb_30cm": sum(v >= 30 for v in val.values()),
                     "premsa": {n: (round(val[n]), round(rang[n])) for n in citats}}
        print(ep, resum[ep])
    with open(d / "model" / "sim2d_carrers.csv", "w", newline="") as f:
        w = csv.DictWriter(f, ["episodi", "carrer", "calat_cm", "hora_utc"]); w.writeheader(); w.writerows(files)
    (d / "model" / "sim2d_premsa.json").write_text(json.dumps(resum, ensure_ascii=False, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
