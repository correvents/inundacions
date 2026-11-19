#!/usr/bin/env python3
"""Radar corregit amb els pluviòmetres propers, en directe (es crida a cada pas de vigila.yml, després de radar_conca.py).

Per a cada hora de les últimes --hores: factor F = Σ pluviòmetres / Σ radar al punt de cada pluviòmetre (font radar_punt,
que desa radar_conca.py), entre 0,25 i 4, només si hi ha ≥ 3 pluviòmetres i plou prou (Σ ≥ 1 mm); si no, 1.
L'hora en curs, que encara no té pluviòmetres complets, fa servir el factor de l'última hora vàlida (com simula_2d.py).
Desa a pluja_conca la font radar_corregit: les mateixes files que radar (cada subconca, cada 5 min) multiplicades per F.
Es reescriu a cada pas (upsert), de manera que quan arriba la XEMA (amb ~1 h de retard) el factor s'afina.
És el mateix biaix mitjà que correccio_cat.py valida a l'històric (octubre 2026: biaix del radar -26% → +1%).

Ús:  python3 scripts/radar_corregit.py --config config/ridaura.json [--hores 3]"""
import argparse, json, sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from sb import insereix, llegeix
import pluvios_propers as pp


def tot(params):
    out, off = [], 0
    while True:
        r = llegeix("pluja_conca", {**params, "order": "temps", "offset": off, "limit": 1000})
        out += r; off += 1000
        if len(r) < 1000: return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--hores", type=int, default=3)
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text()); conca = cfg["conca"]
    ara = datetime.now(timezone.utc).replace(minute=0, second=0, microsecond=0) + timedelta(hours=1)   # fi de l'hora en curs
    ini = ara - timedelta(hours=a.hores + 1)
    H = a.hores + 1
    rang = {"and": f"(temps.gt.{ini.isoformat()},temps.lte.{ara.isoformat()})"}
    radar = tot({"select": "subconca,temps,minuts,mm,mmh_max,cobertura", "conca": f"eq.{conca}", "font": "eq.radar", **rang})
    if not radar:
        print("Radar corregit: no hi ha radar a la finestra"); return 0
    punts = tot({"select": "subconca,temps,mm", "conca": f"eq.{conca}", "font": "eq.radar_punt", **rang})
    hora = lambda t: int(np.ceil((datetime.fromisoformat(t).astimezone(timezone.utc) - ini).total_seconds() / 3600)) - 1
    R = {}
    for r in punts:
        h = hora(r["temps"])
        if 0 <= h < H and r["mm"] is not None:
            R.setdefault(r["subconca"][2:], np.zeros(H))[h] += r["mm"]
    G = pp.horaria(pp.propers(cfg), ini, H)
    F = np.ones(H); valid = np.zeros(H, bool)
    for h in range(H - 1):                         # l'última hora és la que està en curs
        g = [G[k][h] for k in G if k in R and np.isfinite(G[k][h])]
        r = [R[k][h] for k in G if k in R and np.isfinite(G[k][h])]
        if len(g) >= 3 and (sum(g) >= 1 or sum(r) >= 1):
            F[h] = float(np.clip(sum(g) / max(sum(r), 0.2), 0.25, 4)); valid[h] = True
    for h in range(1, H):                           # hores sense prou dades i l'hora en curs: l'últim factor vàlid
        if not valid[h] and (h == H - 1 or not any(np.isfinite(G[k][h]) for k in G)):
            F[h] = F[h - 1] if (valid[h - 1] or h == H - 1) else 1.0
    files = []
    for r in radar:
        h = hora(r["temps"])
        if not 0 <= h < H: continue
        f = F[h]
        files.append({**r, "conca": conca, "font": "radar_corregit",
                      "mm": None if r["mm"] is None else round(r["mm"] * f, 3),
                      "mmh_max": None if r["mmh_max"] is None else round(r["mmh_max"] * f, 2)})
    insereix("pluja_conca", files, actualitza=True)
    print(f"Radar corregit: {len(files)} files; {len(G)} pluviòmetres amb dades; factor per hora " +
          ", ".join(f"{(ini + timedelta(hours=h + 1)):%H}h {F[h]:.2f}{'' if valid[h] else '*'}" for h in range(H)) + " (* = heretat)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
