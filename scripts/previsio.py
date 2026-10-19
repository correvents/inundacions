#!/usr/bin/env python3
"""Previsió de cabal en temps real (cada 10 min): Ridaura a l'aforador i a la desembocadura, i riera de Fenals.

1. Llegeix de Supabase els últims 7 dies de pluja (pluviòmetre de l'ACA a l'aforador; els forats s'omplen amb el radar
   OPERA sobre la conca) i de cabal observat.
2. Busca l'inici de l'episodi de pluja: després de l'última finestra de 48 h amb menys d'1 mm. La pluja acumulada des
   d'aleshores és la que omple el llindar Ia del model (el 4/10 la conca no va respondre fins a ~150 mm).
3. Simula amb el model calibrat (model_cabal.py, pluviòmetre) fins a 3 h endavant en dos escenaris:
     - aforador: no plou més,
     - aforador_continua: plou 1 h més amb la intensitat de l'última hora.
   La diferència entre el cabal observat i el simulat ara es manté a la previsió i s'esvaeix en 2 h.
4. Desembocadura: mateix període de retorn que a l'aforador, 45 min més tard (com a episodi.py).
   Fenals: mètode racional amb la mateixa pluja (episodi.py).
5. Desa a previsions_cabal (des de 6 h enrere fins a 3 h endavant), un avís quan canvia el nivell, i esborra
   les previsions de més de 2 dies.

Ús:  python3 scripts/previsio.py --config config/ridaura.json"""
import argparse, json, sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from sb import insereix, llegeix, esborra
from model_cabal import simula, DT
from episodi import periode_retorn, cabal_racional, FENALS

DIES, ENDAVANT_H, ENRERE_H, RETARD_MAR = 7, 3, 6, 45
NIVELLS = [(500, 5, "T500"), (100, 4, "T100"), (50, 3, "T50"), (10, 2, "T10"), (2.33, 1, "crescuda ordinària")]


def tot(taula, params):
    out, off = [], 0
    while True:
        r = llegeix(taula, {**params, "order": "temps", "offset": off, "limit": 1000})
        out += r; off += 1000
        if len(r) < 1000:
            return out


def serie(files, temps, camp, f=lambda v: v):
    s = np.full(len(temps), np.nan); ix = {t: i for i, t in enumerate(temps)}
    for r in files:
        t = datetime.fromisoformat(r["temps"]).astimezone(timezone.utc).replace(second=0, microsecond=0)
        if t in ix and r[camp] is not None:
            s[ix[t]] = f(r[camp])
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--ara", help="fa la previsió com si fos aquesta hora UTC (YYYY-MM-DDTHH:MM), per provar amb episodis passats")
    ap.add_argument("--sense-desar", action="store_true")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    mod = json.loads(Path(f"config/{cfg['conca']}_model_pluviometre.json").read_text())
    conca = cfg["conca"]

    ara = (datetime.fromisoformat(a.ara).replace(tzinfo=timezone.utc) if a.ara else datetime.now(timezone.utc)).replace(second=0, microsecond=0)
    ara -= timedelta(minutes=ara.minute % DT)
    ini = ara - timedelta(days=DIES)
    temps = [ini + timedelta(minutes=DT * i) for i in range(int((ara - ini).total_seconds() // (DT * 60)) + 1)]
    rang = {"and": f"(temps.gt.{(ini - timedelta(minutes=1)).isoformat()},temps.lte.{ara.isoformat()})"}
    obs = {v: tot("observacions", {"select": "temps,valor", "conca": f"eq.{conca}", "var": f"eq.{v}", **rang})
           for v in ("pluja_mmh", "cabal_m3s")}
    rad = tot("pluja_conca", {"select": "temps,mm", "conca": f"eq.{conca}", "subconca": "eq.aforador", "font": "eq.radar", **rang})
    pluvio = serie(obs["pluja_mmh"], temps, "valor", lambda v: v * DT / 60)
    radar = serie(rad, temps, "mm")
    q_obs = serie(obs["cabal_m3s"], temps, "valor")
    pluja = np.where(np.isnan(pluvio), np.nan_to_num(radar), pluvio)

    # l'última dada de l'aforador marca el "ara" del model (l'ACA publica amb uns minuts de retard)
    ok = np.flatnonzero(~np.isnan(q_obs))
    if not len(ok):
        sys.exit("Sense cabal observat als últims dies")
    i_ara = int(ok[-1])
    pluja = pluja[: i_ara + 1]
    sec = np.convolve(pluja, np.ones(48 * 60 // DT), mode="full")[: len(pluja)]   # pluja de les 48 h anteriors a cada pas
    secs = np.flatnonzero(sec < 1.0)
    i0 = int(secs[-1]) + 1 if len(secs) else 0
    i0 = min(i0, i_ara)
    q_base = float(np.nanmin(q_obs[max(i0 - 12, 0): i0 + 1])) if np.isfinite(q_obs[max(i0 - 12, 0): i0 + 1]).any() else float(q_obs[i_ara])
    n_fut = ENDAVANT_H * 60 // DT
    darrera_h = pluja[max(i_ara - 11, 0): i_ara + 1].mean()
    escenaris = {"aforador": np.zeros(n_fut), "aforador_continua": np.r_[np.full(12, darrera_h), np.zeros(n_fut - 12)]}
    p = dict(ia=mod["Ia_mm"], c=mod["C"], k=mod["K"], m=mod["m"], n=mod["n"])
    acum = float(pluja[i0:].sum())
    t_run = temps[i_ara]
    t_ep = [temps[i0] + timedelta(minutes=DT * k) for k in range(i_ara - i0 + 1 + n_fut)]
    desa = [t for t in t_ep if t >= t_run - timedelta(hours=ENRERE_H)]
    k_desa = len(t_ep) - len(desa)
    files, resum = [], {}
    for nom, fut in escenaris.items():
        q = simula(np.r_[pluja[i0:], fut], mod["area_km2"], q_base, **p)
        k_ara = i_ara - i0
        offset = q_obs[i_ara] - q[k_ara]                       # assimilació: el model parteix del cabal mesurat ara
        dt_min = (np.arange(len(q)) - k_ara) * DT
        q = np.where(dt_min >= 0, np.maximum(q + offset * np.exp(-np.maximum(dt_min, 0) / 120), 0), q)
        resum[nom] = q
        files += [{"conca": conca, "punt": nom, "run": t_run.isoformat(), "temps": t.isoformat(), "cabal_m3s": round(float(v), 2)}
                  for t, v in zip(desa, q[k_desa:])]
    # desembocadura: mateix període de retorn, 45 min més tard; Fenals: mètode racional
    qt_af, qt_mar = cfg["cabals_T"]["aforador"], cfg["cabals_T"]["desembocadura"]
    q_af = resum["aforador"]; r = RETARD_MAR // DT
    t_af = np.array([periode_retorn(v, qt_af) for v in q_af])
    t_mar = np.r_[np.zeros(r), t_af[: len(t_af) - r]]
    q_af_ret = np.r_[np.full(r, q_af[0]), q_af[: len(q_af) - r]]
    q_mar = np.where(t_mar > 0, np.interp(np.log(np.maximum(t_mar, 2.33)), np.log([2.33, 10, 50, 100, 500]),
                                          [qt_mar[k] for k in ("QMCO", "Q10", "Q50", "Q100", "Q500")]),
                     q_af_ret * qt_mar["QMCO"] / qt_af["QMCO"])
    q_fen = cabal_racional(np.r_[pluja[i0:], escenaris["aforador"]], FENALS["area_km2"], 26, 45)
    for nom, q in (("desembocadura", q_mar), ("fenals", q_fen)):
        files += [{"conca": conca, "punt": nom, "run": t_run.isoformat(), "temps": t.isoformat(), "cabal_m3s": round(float(v), 2)}
                  for t, v in zip(desa, q[k_desa:])]
    ins, esb = ((lambda *x, **k: None),) * 2 if a.sense_desar else (insereix, esborra)
    ins("previsions_cabal", files, actualitza=True)
    esb("previsions_cabal", {"conca": f"eq.{conca}", "run": f"lt.{(t_run - timedelta(days=2)).isoformat()}"})

    # avís: el pitjor període de retorn previst a les 3 h vinents (sense més pluja / si continua)
    fut = slice(i_ara - i0, None)
    pitjor = max(max(periode_retorn(v, qt_af) for v in resum[n][fut]) for n in escenaris)
    pitjor = max(pitjor, max(periode_retorn(v, FENALS) for v in q_fen[fut]))
    nivell, motiu = next(((n, f"Previsió fins a {m} en les 3 h vinents") for T, n, m in NIVELLS if pitjor >= T), (0, "Sense crescuda prevista"))
    darrer = llegeix("avisos", {"select": "nivell", "conca": f"eq.{conca}", "order": "temps.desc", "limit": 1})
    if not darrer or darrer[0]["nivell"] != nivell:
        ins("avisos", [{"conca": conca, "temps": t_run.isoformat(), "nivell": nivell,
                             "motiu": f"{motiu}. Pluja de l'episodi: {acum:.0f} mm (llindar de resposta {mod['Ia_mm']:.0f} mm)."}])
    print(f"Run {t_run:%Y-%m-%d %H:%M} UTC · episodi des de {temps[i0]:%Y-%m-%d %H:%M} · pluja {acum:.1f} mm · "
          f"Q ara {q_obs[i_ara]:.1f} · màx 3 h {resum['aforador'][fut].max():.1f} (si continua {resum['aforador_continua'][fut].max():.1f}) m³/s · "
          f"Fenals {q_fen[fut].max():.1f} · nivell {nivell}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
