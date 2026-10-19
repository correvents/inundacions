#!/usr/bin/env python3
"""Característiques fixes de cada carrer de la plana de Platja d'Aro per al model d'inundació (el que no depèn de la pluja).

Per a cada carrer amb nom de l'OpenStreetMap dins la plana, sobre la graella de 20 m (graella_web.py):
  conca_km2      àrea màxima que drena cap al carrer (aigua que li arriba de dalt)
  pendent_pct    pendent mitjà del terreny
  clotada_m      fondària màxima de clotada (relleu omplert − terreny): on s'embassa l'aigua sense sortida
  hand_m         alçada mínima sobre el curs d'aigua més proper aigües avall (HAND): com de baix queda respecte a riu i rieres
  dist_ridaura_m distància al Ridaura
  aca_T10/T100/T500  part del carrer dins la zona inundable de l'ACA
  urba           part urbanitzada
  punt_baix      pla (< 3,5 %), en clotada (≥ 0,4 m) i a ≤ 1 m sobre el curs d'aigua: on s'inunda amb els xàfecs. Regla ajustada
                 amb la premsa 2018-2026: els 6 carrers que s'inunden per xàfec hi són tots, i només marca el 8 % dels carrers
Escriu data/<conca>/model/carrers.csv i carrers.geojson.

Ús:  python3 scripts/punts_model.py --config config/ridaura.json"""
import argparse, json, sys, csv
from pathlib import Path
import numpy as np, shapely
from shapely.geometry import shape, mapping
from shapely.ops import unary_union

CAIXA = (3.02, 41.79, 3.085, 41.84)
DR, DC = [-1, -1, 0, 1, 1, 1, 0, -1], [0, 1, 1, 1, 0, -1, -1, -1]
RIU_KM2 = 0.5                                   # a partir d'aquí la cel·la és un curs d'aigua (riera)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    d = Path("data") / cfg["conca"]
    G = json.loads((d / "web" / "graella.json").read_text()); bin_ = (d / "web" / "graella.bin").read_bytes()
    nf, nc, pas = G["nf"], G["nc"], G["pas_m"]
    cap = lambda nom: next(np.frombuffer(bin_, np.dtype(c["tipus"]), c["bytes"] // np.dtype(c["tipus"]).itemsize, c["offset"]).reshape(nf, nc)
                           for c in G["capes"] if c["nom"] == nom)
    z, zf, dr, urba = cap("z_cm").astype(float), cap("zomplert_cm").astype(float), cap("dir"), cap("urba")
    terra = z != 65535; z /= 100; zf /= 100
    # receptors i drenatge acumulat
    N = nf * nc; rcv = np.full(N, -1)
    rr, cc = np.divmod(np.arange(N), nc); dd = dr.ravel()
    ok = (dd <= 7) & terra.ravel()
    r2 = rr + np.where(ok, np.take(DR, np.minimum(dd, 7)), 0); c2 = cc + np.where(ok, np.take(DC, np.minimum(dd, 7)), 0)
    dins = ok & (r2 >= 0) & (r2 < nf) & (c2 >= 0) & (c2 < nc)
    rcv[dins] = (r2 * nc + c2)[dins]; rcv[rcv >= 0] = np.where(terra.ravel()[rcv[rcv >= 0]], rcv[rcv >= 0], -1)
    entra = np.bincount(rcv[rcv >= 0], minlength=N); ordre = list(np.flatnonzero(entra == 0)); k = 0
    entra = entra.copy()
    while k < len(ordre):
        j = rcv[ordre[k]]; k += 1
        if j >= 0:
            entra[j] -= 1
            if entra[j] == 0: ordre.append(j)
    acc = np.ones(N)
    for i in ordre:
        if rcv[i] >= 0: acc[rcv[i]] += acc[i]
    km2 = acc * pas * pas / 1e6
    # HAND: alçada sobre la primera cel·la de curs d'aigua aigües avall (recorregut invers de l'ordre)
    riu = km2 >= RIU_KM2; zr = z.ravel(); base = np.where(riu, zr, np.nan)
    for i in reversed(ordre):
        if not riu[i] and rcv[i] >= 0: base[i] = base[rcv[i]]
    hand = (zr - base).reshape(nf, nc)
    gy, gx = np.gradient(np.where(terra, z, np.nan), pas); pendent = 100 * np.hypot(gx, gy)
    clot = np.where(terra, zf - z, 0)
    # coordenades del centre de cada cel·la
    (s, w), (n, e) = G["limits"]
    LA, LO = np.meshgrid(n - (np.arange(nf) + .5) * (n - s) / nf, w + (np.arange(nc) + .5) * (e - w) / nc, indexing="ij")
    caixa = (LO > CAIXA[0]) & (LO < CAIXA[2]) & (LA > CAIXA[1]) & (LA < CAIXA[3])
    zones = {k: unary_union([shape(f["geometry"]) for f in json.loads((d / "zones" / f"aca_{k}.geojson").read_text())["features"]]).buffer(0)
             for k in ("T10", "T100", "T500")}
    zm = {k: shapely.contains_xy(g, LO, LA) for k, g in zones.items()}
    ridaura = (km2 >= 30).reshape(nf, nc)
    from scipy.ndimage import distance_transform_edt
    dist_r = distance_transform_edt(~ridaura) * pas
    osm = json.loads((d / "zones" / "osm.geojson").read_text())["features"]
    per = {}
    for f in osm:
        nom = f["properties"].get("name")
        if not nom or not f["geometry"]["type"].endswith("LineString"): continue
        g = shape(f["geometry"]); c = g.centroid
        if CAIXA[0] < c.x < CAIXA[2] and CAIXA[1] < c.y < CAIXA[3]: per.setdefault(nom, []).append(g)
    files, feats = [], []
    K = km2.reshape(nf, nc)
    for nom, gs in sorted(per.items()):
        g = unary_union(gs); m = shapely.contains_xy(g.buffer(0.00012), LO, LA) & terra & caixa
        if m.sum() < 2: continue
        r = {"nom": nom, "cel": int(m.sum()), "conca_km2": round(float(K[m].max()), 3),
             "pendent_pct": round(float(np.nanmean(pendent[m])), 2), "clotada_m": round(float(clot[m].max()), 2),
             "hand_m": round(float(np.nanmin(hand[m])), 2), "dist_ridaura_m": round(float(dist_r[m].min())),
             **{f"aca_{k}": round(float(zm[k][m].mean()), 2) for k in zm}, "urba": round(float((urba[m] > 0).mean()), 2)}
        r["punt_baix"] = bool(r["pendent_pct"] < 3.5 and r["clotada_m"] >= 0.4 and r["hand_m"] <= 1.0)
        files.append(r); feats.append({"type": "Feature", "properties": r, "geometry": mapping(g.simplify(0.00002))})
    with open(d / "model" / "carrers.csv", "w", newline="") as f:
        w_ = csv.DictWriter(f, files[0].keys()); w_.writeheader(); w_.writerows(files)
    (d / "model" / "carrers.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": feats}, ensure_ascii=False, separators=(",", ":")))
    print(len(files), "carrers")
    return 0


if __name__ == "__main__":
    sys.exit(main())
