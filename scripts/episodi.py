#!/usr/bin/env python3
"""Reconstrucció pas a pas d'un episodi per a l'animació del visor: cabal, període de retorn i què s'inunda a cada moment.

Les zones inundables de l'ACA de Platja d'Aro venen de dos estudis amb motors diferents:
  - codi 009, el Ridaura: es mou amb el cabal observat a l'aforador de Santa Cristina. Aigües avall el període de retorn
    és el mateix però arriba amb retard, proporcional a la distància fins a la desembocadura (--retard minuts).
  - codi 640, la plana de Platja d'Aro (riera de Fenals, ~3,5 km² molt urbanitzats): no té aforador. El cabal s'estima amb
    el mètode racional: pluja neta (llindar P0 del mètode SCS) mitjana en el temps de concentració × superfície.
Cada element s'inunda quan el període de retorn del seu motor arriba al de la primera zona que el toca, i deixa d'estar-ho
quan baixa per sota.

Ús:  python3 scripts/episodi.py --config config/ridaura.json --csv data/ridaura/episodi_2026-10-02.csv \\
         --des 2026-10-03T18:00 --fins 2026-10-04T14:00
Escriu data/<conca>/episodis/<data de --fins>.json."""
import argparse, csv, json, sys
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
from shapely.geometry import shape, mapping
from shapely.ops import unary_union

DT = 5
LON_AFORADOR, LON_MAR = 2.979, 3.067
# Riera de Fenals a Platja d'Aro (ACA, tram FENALS020, estudi HIR640)
FENALS = {"area_km2": 3.53, "QMCO": 15.6, "Q10": 40.52, "Q50": 71.35, "Q100": 86.74, "Q500": 126.52}


def periode_retorn(q, qt):
    """Interpolació en log T entre els cabals de l'ACA (com al visor). 0 = per sota de la crescuda ordinària."""
    p = [(2.33, qt["QMCO"]), (10, qt["Q10"]), (50, qt["Q50"]), (100, qt["Q100"]), (500, qt["Q500"])]
    if not np.isfinite(q) or q < p[0][1]:
        return 0.0
    for (t0, q0), (t1, q1) in zip(p, p[1:]):
        if q <= q1:
            return float(np.exp(np.log(t0) + (q - q0) / (q1 - q0) * (np.log(t1) - np.log(t0))))
    return 501.0


def cabal_racional(pluja, area_km2, p0, tc_min):
    """Cabal (m³/s): pluja neta SCS amb llindar p0, mitjana mòbil en el temps de concentració."""
    p = np.cumsum(pluja)
    neta = np.diff(np.where(p > p0, (p - p0) ** 2 / (p + 4 * p0), 0.0), prepend=0.0)
    n = max(int(round(tc_min / DT)), 1)
    i_mmh = np.convolve(neta, np.ones(n) / n, mode="full")[: len(neta)] * 60 / DT
    return i_mmh * area_km2 / 3.6


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--csv", required=True, help="episodi exportat amb model_cabal.py --exporta")
    ap.add_argument("--des", required=True)
    ap.add_argument("--fins", required=True)
    ap.add_argument("--retard", type=float, default=45, help="minuts que triga la punta de l'aforador a la desembocadura")
    ap.add_argument("--p0-fenals", type=float, default=26, help="llindar d'escorrentia (mm) de la conca urbana de Fenals")
    ap.add_argument("--tc-fenals", type=float, default=45, help="temps de concentració de Fenals (min)")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    d = Path("data") / cfg["conca"]

    files = {r["temps"]: r for r in csv.DictReader(open(a.csv))}
    ini, fi = datetime.fromisoformat(a.des), datetime.fromisoformat(a.fins)
    temps = [ini + timedelta(minutes=DT * i) for i in range(int((fi - ini).total_seconds() // (DT * 60)) + 1)]
    claus = [t.strftime("%Y-%m-%dT%H:%M") for t in temps]
    num = lambda k, c: float(files[k][c]) if k in files and files[k][c] not in ("", None) else np.nan
    pluvio = np.nan_to_num(np.array([num(k, "pluja_pluviometre_mm") for k in claus]))
    radar = np.nan_to_num(np.array([num(k, "pluja_radar_mm") for k in claus]))
    q_af = np.array([num(k, "cabal_m3s") for k in claus])
    q_af = np.where(np.isnan(q_af), np.interp(np.arange(len(q_af)), np.flatnonzero(~np.isnan(q_af)), q_af[~np.isnan(q_af)]), q_af)

    qt_af, qt_mar = cfg["cabals_T"]["aforador"], cfg["cabals_T"]["desembocadura"]
    t_af = np.array([periode_retorn(q, qt_af) for q in q_af])
    # la pluja de tota l'avantconca s'afegeix aigües avall: a la desembocadura es manté el període de retorn, no el cabal
    pas_retard = int(round(a.retard / DT))
    t_mar = np.concatenate([np.zeros(pas_retard), t_af[: len(t_af) - pas_retard]])
    q_mar = np.array([np.interp(np.log(max(t, 2.33)), np.log([2.33, 10, 50, 100, 500]),
                                [qt_mar[k] for k in ("QMCO", "Q10", "Q50", "Q100", "Q500")]) if t > 0 else np.nan for t in t_mar])
    q_fen = cabal_racional(pluvio, FENALS["area_km2"], a.p0_fenals, a.tc_fenals)
    t_fen = np.array([periode_retorn(q, FENALS) for q in q_fen])

    # zones acumulades de cada motor (T10 ⊂ T50 ⊂ T100 ⊂ T500), per dibuixar-les a l'animació
    zones, motors = [], {"ridaura": "009", "fenals": "640"}
    for motor, codi in motors.items():
        acum = None
        for T in (10, 50, 100, 500):
            g = unary_union([shape(f["geometry"]).buffer(0) for f in json.loads((d / "zones" / f"aca_T{T}.geojson").read_text())["features"]
                             if f["properties"].get("CODI") == codi])
            acum = g if acum is None else acum.union(g)
            zones.append({"type": "Feature", "properties": {"motor": motor, "T": T},
                          "geometry": mapping(acum.simplify(0.00003, preserve_topology=True))})
    (d / "zones_per_motor.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": zones[::-1]}, separators=(",", ":")))
    # motor de cada element: zona 640 (Fenals) o 009 (Ridaura, amb retard segons la longitud)
    z640 = shape([z for z in zones if z["properties"] == {"motor": "fenals", "T": 500}][0]["geometry"]).buffer(0.0002)
    af = json.loads((d / "afectats.geojson").read_text())["features"]
    elements, retards = [], []
    for f in af:
        g = shape(f["geometry"]); T = f["properties"]["T"]
        if z640.intersects(g):
            serie, motor = t_fen, "fenals"
        else:
            frac = min(max((g.centroid.x - LON_AFORADOR) / (LON_MAR - LON_AFORADOR), 0), 1)
            r = int(round(frac * a.retard / DT))
            serie, motor = np.concatenate([np.zeros(r), t_af[: len(t_af) - r]]), "ridaura"
            retards.append(r)
        dins = np.flatnonzero(serie >= T)
        elements.append({"motor": motor, "ini": int(dins[0]) if len(dins) else None, "fi": int(dins[-1]) if len(dins) else None})

    r1 = lambda v: [None if not np.isfinite(x) else round(float(x), 1) for x in v]
    out = {"titol": f"Episodi del {fi:%d/%m/%Y}", "temps_utc": claus, "pas_min": DT,
           "pluja_pluviometre_mm": r1(pluvio), "pluja_radar_conca_mm": [round(float(x), 2) for x in radar],
           "q_aforador_obs": r1(q_af), "T_aforador": r1(t_af), "q_desembocadura_est": r1(q_mar), "T_desembocadura": r1(t_mar),
           "q_fenals_est": r1(q_fen), "T_fenals": r1(t_fen), "elements": elements,
           # la zona del Ridaura es dibuixa amb el retard mitjà dels elements que hi ha
           "T_ridaura_zona": r1(np.concatenate([np.zeros(int(np.median(retards))), t_af[: len(t_af) - int(np.median(retards))]])),
           "hipotesis": {"retard_aforador_mar_min": a.retard, "p0_fenals_mm": a.p0_fenals, "tc_fenals_min": a.tc_fenals,
                         "cabals_T_fenals": FENALS}}
    dest = d / "episodis"; dest.mkdir(parents=True, exist_ok=True)
    (dest / f"{fi:%Y-%m-%d}.json").write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")))
    for nom, s, q in (("Ridaura aforador", t_af, q_af), ("Ridaura mar", t_mar, q_mar), ("Fenals", t_fen, q_fen)):
        i = int(np.nanargmax(q)); print(f"{nom}: pic {np.nanmax(q):.1f} m³/s a les {claus[i]} UTC, T màx {s.max():.0f}")
    n = sum(e["ini"] is not None for e in elements)
    print(f"elements inundats: {n} de {len(elements)} "
          f"(Ridaura {sum(e['ini'] is not None for e in elements if e['motor'] == 'ridaura')}, Fenals {sum(e['ini'] is not None for e in elements if e['motor'] == 'fenals')})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
