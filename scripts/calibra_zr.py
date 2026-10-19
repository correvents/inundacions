#!/usr/bin/env python3
"""Ajusta la relació reflectivitat-pluja (Z = a·R^b) del radar OPERA a la plana de Platja d'Aro amb les estacions XEMA.

El radar es baixa convertit amb Marshall-Palmer (Z = 200 R^1.6). Aquí es desfà la conversió (R → Z) i es busca una
relació per a cada tipus de pluja, separats per la reflectivitat de cada píxel i fotograma:
  estratiforme (llevantades, pluja d'hivern; Z < LLINDAR dBZ): Z = a_s R^1.6
  convectiva   (tempestes, calamarsa; Z ≥ LLINDAR dBZ):           Z = a_c R^1.4
Es minimitza l'error (log) entre la pluja total de cada episodi a cada estació i la del radar. Escriu config/<conca>_zr.json.

Ús:  python3 scripts/calibra_zr.py --config config/ridaura.json"""
import argparse, json, math, sys, importlib.util
from pathlib import Path
import numpy as np, pandas as pd

LLINDAR_DBZ = 38.0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    d = Path("data") / json.loads(Path(a.config).read_text())["conca"]
    spec = importlib.util.spec_from_file_location("premsa", "scripts/premsa.py"); P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)
    estacions = json.loads((d / "xema" / "estacions.json").read_text())
    mostres = []                                       # (episodi, estació, mm estació, Z de cada fotograma, minuts)
    for ep in P.EPISODIS:
        f = d / "episodis" / ep["radar"]
        if not f.exists(): continue
        R = json.loads(f.read_text()); LO, LA = np.meshgrid(R["lons"], R["lats"]); mm = np.array(R["mm"])
        t = pd.to_datetime(R["temps"]); dt = R["interval_min"]
        for e in estacions:
            fe = d / "xema" / f"{e['codi_estacio']}.csv"
            if not fe.exists(): continue
            lat, lon = float(e["latitud"]), float(e["longitud"])
            km = np.hypot((LA - lat) * 111, (LO - lon) * 111 * math.cos(math.radians(lat))).ravel()
            if not (km <= 3).any(): continue
            s = pd.read_csv(fe, parse_dates=["temps"]).set_index("temps").mm
            s.index = s.index + pd.Timedelta(minutes=30)
            sel = s[(s.index > t[0] - pd.Timedelta(minutes=dt)) & (s.index <= t[-1])]
            if len(sel) < 0.8 * (len(t) * dt / 30) or sel.sum() < 10: continue
            rate = mm[:, km <= 3] * 60 / dt                        # mm/h de cada píxel i fotograma (Marshall-Palmer)
            Z = 200 * np.maximum(rate, 1e-6) ** 1.6
            mostres.append((ep["id"], e["codi_estacio"], float(sel.sum()), Z, rate > 0, dt))
    def total(Z, plou, dt, a_s, a_c):
        dbz = 10 * np.log10(Z)
        r = np.where(dbz >= LLINDAR_DBZ, (Z / a_c) ** (1 / 1.4), (Z / a_s) ** (1 / 1.6))
        r = np.minimum(np.where(plou, r, 0), 125)
        return float(r.mean(axis=1).sum() * dt / 60)
    def cost(a_s, a_c):
        return sum(math.log((total(Z, p, dt, a_s, a_c) + 1) / (obs + 1)) ** 2 for _, _, obs, Z, p, dt in mostres)
    millor = min(((cost(a_s, a_c), a_s, a_c) for a_s in (30, 50, 75, 100, 130, 160, 200, 250, 300)
                  for a_c in (200, 300, 400, 500, 700, 1000, 1500, 2000, 3000)))
    _, a_s, a_c = millor
    print(f"{len(mostres)} episodi-estació. Millor: Z = {a_s} R^1.6 (estratiforme), Z = {a_c} R^1.4 (convectiva, ≥ {LLINDAR_DBZ} dBZ)")
    files = []
    for ep, est, obs, Z, p, dt in mostres:
        mp, nou = total(Z, p, dt, 200, 200 * 1.0), total(Z, p, dt, a_s, a_c)
        mp = float((np.where(p, (Z / 200) ** (1 / 1.6), 0)).mean(axis=1).sum() * dt / 60)
        files.append((ep, est, round(obs), round(mp), round(nou))); print(f"  {ep} {est}: estació {obs:.0f}, Marshall-Palmer {mp:.0f}, nova {nou:.0f} mm")
    err = lambda k: float(np.median([abs(math.log((f[k] + 1) / (f[2] + 1))) for f in files]))
    print(f"error típic (factor): Marshall-Palmer ×{math.exp(err(3)):.2f}, nova ×{math.exp(err(4)):.2f}")
    (Path("config") / f"{d.name}_zr.json").write_text(json.dumps({
        "descripcio": "Relació Z-R del radar OPERA ajustada amb les estacions XEMA de Castell d'Aro i Cassà en els episodis d'inundació i pluja forta (calibra_zr.py)",
        "llindar_dbz": LLINDAR_DBZ, "estratiforme": {"a": a_s, "b": 1.6}, "convectiva": {"a": a_c, "b": 1.4},
        "error_tipic_factor": {"marshall_palmer": round(math.exp(err(3)), 2), "ajustada": round(math.exp(err(4)), 2)},
        "mostres": [dict(zip(("episodi", "estacio", "mm_estacio", "mm_marshall_palmer", "mm_ajustada"), f)) for f in files]}, ensure_ascii=False, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
