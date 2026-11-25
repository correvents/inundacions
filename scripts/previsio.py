#!/usr/bin/env python3
"""Previsió de cabal en temps real (cada 10 min): Ridaura a l'aforador i a la desembocadura, i riera de Fenals.

1. Llegeix de Supabase els últims 7 dies de pluja (pluviòmetre de l'ACA a l'aforador; els forats s'omplen amb el radar
   OPERA sobre la conca) i de cabal observat.
2. Busca l'inici de l'episodi de pluja: després de l'última finestra de 48 h amb menys d'1 mm. La pluja acumulada des
   d'aleshores és la que omple el llindar Ia del model. Ia depèn de la pluja dels 60 dies anteriors (humitat.py: dipòsit
   d'humitat del sòl calibrat amb el 4/10, quan després de l'estiu sec la conca no va respondre fins a ~160 mm). Es calcula
   un cop per inici d'episodi i es guarda a /tmp mentre dura l'execució de vigila.yml.
3. Simula amb el model calibrat (humitat.py si hi ha config/<conca>_humitat.json; si no, model_cabal.py, pluviòmetre) fins a 3 h endavant en dos escenaris:
     - aforador: no plou més,
     - aforador_continua: plou 1 h més amb la intensitat de l'última hora.
   La diferència entre el cabal observat i el simulat ara es manté a la previsió i s'esvaeix en 2 h.
4. Desembocadura: mateix període de retorn que a l'aforador, 45 min més tard (com a episodi.py).
   Fenals: mètode racional amb la mateixa pluja (episodi.py).
5. Desa a previsions_cabal (des de 6 h enrere fins a 3 h endavant), un avís quan canvia el nivell, i esborra
   les previsions de més de 2 dies.

Opció --pluja-prevista (per defecte 'no', i aleshores tot és igual que abans): afegeix l'escenari aforador_previst,
amb la pluja prevista de les 6 h següents (pluja_prevista.py: nowcast del radar OPERA i AROME HD via Open-Meteo), i
allarga totes les previsions a 6 h. L'avís té en compte també aquest escenari i ho diu al motiu.
Amb --avis-previst possible, en canvi, l'avís principal només compta la pluja ja mesurada, i el que afegeix la pluja
prevista es desa a part com a "avís possible" (avisos, conca "<conca>:possible"): és menys segur (el 6/10/2026 AROME va
preveure 46 mm en 3 h i en van caure 5).

Ús:  python3 scripts/previsio.py --config config/ridaura.json [--pluja-prevista tot]"""
import argparse, json, sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from sb import insereix, llegeix, esborra
from model_cabal import simula, DT
from episodi import periode_retorn, cabal_racional, FENALS
import pluja_prevista as pp
import humitat

DIES, ENDAVANT_H, ENRERE_H, RETARD_MAR = 7, 3, 6, 45
pp_H = 6                     # hores de previsió amb --pluja-prevista
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


def absorcio_inicial(conca, t0, mod):
    """Ia (mm) en començar l'episodi a t0, amb el pluviòmetre de l'ACA dels 60 dies anteriors (sense dada = sec).
    Sense config/<conca>_humitat.json, l'Ia fix del model."""
    par = humitat.parametres(conca)
    if not par:
        return mod["Ia_mm"], "fixa"
    cau = Path(f"/tmp/ia_{conca}.json")
    memo = json.loads(cau.read_text()) if cau.exists() else {}
    if t0.isoformat() in memo:
        return memo[t0.isoformat()], "60 dies"
    import pandas as pd
    ini = t0 - timedelta(days=humitat.DIES)
    files = tot("observacions", {"select": "temps,valor", "conca": f"eq.{conca}", "var": "eq.pluja_mmh",
                                 "and": f"(temps.gte.{ini.isoformat()},temps.lt.{t0.isoformat()})"})
    s = pd.Series([r["valor"] * DT / 60 for r in files], index=pd.to_datetime([r["temps"] for r in files], utc=True).tz_localize(None), dtype=float)
    ph = s.resample("1h").sum() if len(s) else pd.Series(dtype=float)
    ia = round(humitat.ia(t0.replace(tzinfo=None), ph, par), 1)
    memo[t0.isoformat()] = ia; cau.write_text(json.dumps(memo))
    return ia, "60 dies"


def preveu(cfg, ara=None, prevista="no", desar=True, cache=None, verb=True, avis_previst="junt"):
    """Una previsió. ara=None: ara mateix (directe). Torna un resum (nivell, cabals màxims i períodes de retorn per escenari)."""
    mod = json.loads(Path(f"config/{cfg['conca']}_model_pluviometre.json").read_text())
    conca = cfg["conca"]
    directe = ara is None
    ara = (ara or datetime.now(timezone.utc)).replace(second=0, microsecond=0)
    ara -= timedelta(minutes=ara.minute % DT)
    ini = ara - timedelta(days=DIES)
    temps = [ini + timedelta(minutes=DT * i) for i in range(int((ara - ini).total_seconds() // (DT * 60)) + 1)]
    rang = {"and": f"(temps.gt.{(ini - timedelta(minutes=1)).isoformat()},temps.lte.{ara.isoformat()})"}
    obs = {v: tot("observacions", {"select": "temps,valor", "conca": f"eq.{conca}", "var": f"eq.{v}", **rang})
           for v in ("pluja_mmh", "cabal_m3s")}
    rad = tot("pluja_conca", {"select": "temps,mm,font", "conca": f"eq.{conca}", "subconca": "eq.aforador", "font": "in.(radar,radar_corregit)", **rang})
    corr = {r["temps"] for r in rad if r["font"] == "radar_corregit"}       # el radar corregit amb pluviòmetres té preferència
    rad = [r for r in rad if r["font"] == "radar_corregit" or r["temps"] not in corr]
    pluvio = serie(obs["pluja_mmh"], temps, "valor", lambda v: v * DT / 60)
    radar = serie(rad, temps, "mm")
    q_obs = serie(obs["cabal_m3s"], temps, "valor")
    pluja = np.where(np.isnan(pluvio), np.nan_to_num(radar), pluvio)

    # l'última dada de l'aforador marca el "ara" del model (l'ACA publica amb uns minuts de retard)
    ok = np.flatnonzero(~np.isnan(q_obs))
    if not len(ok):
        sys.exit("Sense cabal observat als últims dies")
    i_ara = int(ok[-1])
    # la pluja ja mesurada després de l'última dada de l'aforador (el radar arriba abans) també entra al model, fins a 1 h
    amb_pluja = np.flatnonzero(~np.isnan(pluvio) | ~np.isnan(radar))
    i_fi = min(max(int(amb_pluja[-1]) if len(amb_pluja) else i_ara, i_ara), i_ara + 60 // DT)
    pluja = pluja[: i_fi + 1]
    sec = np.convolve(pluja, np.ones(48 * 60 // DT), mode="full")[: len(pluja)]   # pluja de les 48 h anteriors a cada pas
    secs = np.flatnonzero(sec < 1.0)
    i0 = int(secs[-1]) + 1 if len(secs) else 0
    i0 = min(i0, i_ara)
    q_base = float(np.nanmin(q_obs[max(i0 - 12, 0): i0 + 1])) if np.isfinite(q_obs[max(i0 - 12, 0): i0 + 1]).any() else float(q_obs[i_ara])
    endavant_h = ENDAVANT_H if prevista == "no" else pp_H
    n_fut = endavant_h * 60 // DT
    darrera_h = pluja[max(i_fi - 11, 0): i_fi + 1].mean()
    escenaris = {"aforador": np.zeros(n_fut), "aforador_continua": np.r_[np.full(12, darrera_h), np.zeros(n_fut - 12)]}
    info_prev = None
    if prevista != "no":
        # la pluja prevista comença just després de l'última pluja mesurada (i_fi)
        t_fi = temps[i_fi]
        from radar_conca import relacio_zr
        prev, info_prev = pp.prevista(cfg, t_fi, n_fut * DT, prevista, directe=directe,
                                      zr=relacio_zr(conca) if directe else None, cache=cache)
        escenaris["aforador_previst"] = prev
    ia_ep, via_ia = absorcio_inicial(conca, temps[i0], mod)
    hum = humitat.parametres(conca)
    if hum and "model" in hum:        # model calibrat junt amb el dipòsit d'humitat: Ia variable + pèrdua contínua phi
        mh = hum["model"]
        p = dict(ia=ia_ep, c=mh["C"], k=mh["K"], m=mh["m"], n=mh["n"], phi=mh["phi_mm_h"])
    else:
        p = dict(ia=ia_ep, c=mod["C"], k=mod["K"], m=mod["m"], n=mod["n"])
    acum = float(pluja[i0:].sum())
    t_run = temps[i_ara]
    t_ep = [temps[i0] + timedelta(minutes=DT * k) for k in range(i_fi - i0 + 1 + n_fut)]
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
    ins, esb = (insereix, esborra) if desar else ((lambda *x, **k: None),) * 2
    ins("previsions_cabal", files, actualitza=True)
    esb("previsions_cabal", {"conca": f"eq.{conca}", "run": f"lt.{(t_run - timedelta(days=2)).isoformat()}"})

    # avís: el pitjor període de retorn previst a les 3 h vinents (sense més pluja / si continua)
    fut = slice(i_ara - i0, None)
    T_esc = {n: max(periode_retorn(v, qt_af) for v in resum[n][fut]) for n in escenaris}
    # Fenals: amb la pluja observada i, si n'hi ha, amb la prevista
    q_fen_prev = (cabal_racional(np.r_[pluja[i0:], escenaris["aforador_previst"]], FENALS["area_km2"], 26, 45)
                  if "aforador_previst" in escenaris else q_fen)
    T_fen = max(periode_retorn(v, FENALS) for v in q_fen[fut])
    T_fen_prev = max(periode_retorn(v, FENALS) for v in q_fen_prev[fut])
    pitjor_obs = max(max(T_esc["aforador"], T_esc["aforador_continua"]), T_fen)
    pitjor = max(max(T_esc.values()), T_fen, T_fen_prev)
    nivell, motiu = next(((n, f"Previsió fins a {m} en les {endavant_h} h vinents") for T, n, m in NIVELLS if pitjor >= T), (0, "Sense crescuda prevista"))
    nivell_obs = next((n for T, n, m in NIVELLS if pitjor_obs >= T), 0)
    nivell_possible, motiu_possible = 0, "Sense avís possible per la pluja prevista"
    if nivell > nivell_obs and avis_previst == "possible":
        # l'avís principal es queda amb la pluja mesurada; el que afegeix la pluja prevista va a part
        nivell_possible = nivell
        motiu_possible = (f"Possible {NIVELLS[len(NIVELLS) - nivell][2]} en les {endavant_h} h vinents si cau la pluja prevista "
                          f"({info_prev.get('mm_3h', 0):.0f} mm en 3 h, {info_prev.get('mm_6h', 0):.0f} mm en 6 h)")
        nivell = nivell_obs
        motiu = next((f"Previsió fins a {m} en les {endavant_h} h vinents" for T, n, m in NIVELLS if pitjor_obs >= T), "Sense crescuda prevista")
    elif nivell > nivell_obs:
        motiu += f" amb la pluja prevista ({info_prev.get('mm_3h', 0):.0f} mm en 3 h; sense comptar-la, nivell {nivell_obs})"
    darrer = llegeix("avisos", {"select": "nivell", "conca": f"eq.{conca}", "order": "temps.desc", "limit": 1}) if desar else []
    if desar and (not darrer or darrer[0]["nivell"] != nivell):
        ins("avisos", [{"conca": conca, "temps": t_run.isoformat(), "nivell": nivell,
                             "motiu": f"{motiu}. Pluja de l'episodi: {acum:.0f} mm (el terreny n'absorbeix {ia_ep:.0f} mm, segons la pluja dels 60 dies anteriors"
                                      + (f", i després {p['phi']:.0f} mm/h" if p.get("phi") else "") + ")."}])
    if desar and avis_previst == "possible":
        darrer_p = llegeix("avisos", {"select": "nivell", "conca": f"eq.{conca}:possible", "order": "temps.desc", "limit": 1})
        if (darrer_p[0]["nivell"] if darrer_p else 0) != nivell_possible:
            ins("avisos", [{"conca": f"{conca}:possible", "temps": t_run.isoformat(), "nivell": nivell_possible, "motiu": motiu_possible + "."}])
    if verb:
        print(f"Run {t_run:%Y-%m-%d %H:%M} UTC · episodi des de {temps[i0]:%Y-%m-%d %H:%M} · pluja {acum:.1f} mm · "
              f"Q ara {q_obs[i_ara]:.1f} · màx {endavant_h} h {resum['aforador'][fut].max():.1f} (si continua {resum['aforador_continua'][fut].max():.1f}"
              + (f", amb la prevista {resum['aforador_previst'][fut].max():.1f}" if "aforador_previst" in resum else "")
              + f") m³/s · Fenals {q_fen[fut].max():.1f} · Ia {ia_ep:.0f} mm ({via_ia}) · nivell {nivell}"
              + (f" (possible {nivell_possible})" if nivell_possible else ""))
        if info_prev is not None:
            print("Pluja prevista:", json.dumps(info_prev, ensure_ascii=False))
    return {"run": t_run, "ultima_pluja": temps[i_fi], "nivell": nivell, "nivell_sense_prevista": nivell_obs, "nivell_possible": nivell_possible, "pluja_episodi_mm": acum, "ia_mm": ia_ep,
            "q_ara": float(q_obs[i_ara]), "q_max": {n: float(resum[n][fut].max()) for n in resum},
            "T": {**T_esc, "fenals": T_fen, "fenals_previst": T_fen_prev}, "prevista": info_prev,
            "prevista_mm": escenaris["aforador_previst"] if "aforador_previst" in escenaris else None}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--ara", help="fa la previsió com si fos aquesta hora UTC (YYYY-MM-DDTHH:MM), per provar amb episodis passats")
    ap.add_argument("--sense-desar", action="store_true")
    ap.add_argument("--pluja-prevista", default="no", choices=["no", "nowcast", "arome", "tot"],
                    help="afegeix l'escenari amb la pluja prevista de les 6 h següents (per defecte, no)")
    ap.add_argument("--avis-previst", default="junt", choices=["junt", "possible"],
                    help="junt: la pluja prevista pot pujar l'avís principal; possible: va a un avís a part (menys segur)")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    ara = datetime.fromisoformat(a.ara).replace(tzinfo=timezone.utc) if a.ara else None
    preveu(cfg, ara, a.pluja_prevista, desar=not a.sense_desar, avis_previst=a.avis_previst)
    return 0


if __name__ == "__main__":
    sys.exit(main())
