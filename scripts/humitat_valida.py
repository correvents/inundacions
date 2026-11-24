#!/usr/bin/env python3
"""Validació de l'absorció inicial variable (humitat.py) contra l'absorció fixa, amb el cabal observat.

Ridaura (aforador ACA de Santa Cristina, model de cabal de la previsió, config/ridaura_model_pluviometre.json):
  a) 4/10/2026, l'episodi de calibratge: NSE.
  b) 17-21/1/2026, episodi independent: pic observat 28,6 m³/s.
  c) Sèrie contínua 11/8 – avui (pluviòmetre i cabal de l'ACA cada 5 min), episodi a episodi com la previsió: NSE i
     falses alarmes.
  d) Els episodis de pluja forta de la XEMA 2009-2026: quants passen de la crescuda ordinària (QMCO 37 m³/s); la premsa
     no en recull cap desbordament des del 2005 fora del 4/10/2026.
Congost (aforador ACA de La Garriga): el mateix model ràpid ajustat al 3-5/10/2026 amb l'Ia de cada variant, i
  sèrie contínua juliol-novembre 2026.
Escriu data/<conca>/humitat/validacio.json.

Ús:  python3 scripts/humitat.py valida"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from model_cabal import simula, nse, DT
from humitat import pluja_xema, pluja_conca, ia, parametres

SEC_H = 48


def a5(h):
    """Sèrie horària → passos de 5 min (repartida uniformement dins l'hora)."""
    idx = pd.date_range(h.index[0], h.index[-1] + pd.Timedelta(minutes=55), freq=f"{DT}min")
    return h.reindex(idx, method="ffill") / (60 // DT)


def inici_episodi(p5, t):
    """Com la previsió: just després de l'última finestra de 48 h amb menys d'1 mm abans de t."""
    s = p5[:t]
    sec = s.rolling(SEC_H * 60 // DT, min_periods=1).sum()
    secs = sec[sec < 1.0]
    return (secs.index[-1] + pd.Timedelta(minutes=DT)) if len(secs) else s.index[0]


def episodis(p5, llindar_mm=10):
    """Talla la sèrie en episodis separats per 48 h seques (< 1 mm) amb ≥ llindar_mm."""
    sec = p5.rolling(SEC_H * 60 // DT, min_periods=1).sum()
    out, ini = [], None
    for t, plou in zip(p5.index, p5.values > 0):
        if plou and ini is None:
            ini = t
        if ini is not None and sec[t] < 1.0 and t > ini:
            fi = t
            if p5[ini:fi].sum() >= llindar_mm: out.append((ini, fi))
            ini = None
    if ini is not None and p5[ini:].sum() >= llindar_mm: out.append((ini, p5.index[-1]))
    return out


def simula_ep(p5, t0, fi, ia_mm, mod, q_base=0.0):
    p = p5[t0:fi].values
    return pd.Series(simula(p, mod["area_km2"], q_base, ia_mm, mod["C"], mod["K"], mod["m"], mod["n"], mod.get("phi_mm_h", 0.0)),
                     index=p5[t0:fi].index)


def continu(p5, q, mods, ias, area):
    """Simulació episodi a episodi sobre tot el període amb cabal observat. ias: {variant: funció t0 → Ia}."""
    sims = {v: pd.Series(np.nan, index=q.index) for v in ias}
    det = []
    for t0, fi in episodis(p5[q.index[0]:q.index[-1]]):
        fi2 = min(fi + pd.Timedelta(hours=24), q.index[-1])
        qb = q[t0 - pd.Timedelta(hours=1):t0].min()
        qb = 0.0 if np.isnan(qb) else float(qb)
        fila = {"inici": f"{t0:%Y-%m-%d %H:%M}", "pluja_mm": round(float(p5[t0:fi].sum()), 1), "q_obs_max": round(float(q[t0:fi2].max()), 1)}
        for v, f in ias.items():
            x = f(t0)
            s = simula_ep(p5, t0, fi2, x, mods[v], qb).reindex(q[t0:fi2].index)
            sims[v][s.index] = s.values
            fila[f"ia_{v}"] = round(x); fila[f"q_{v}_max"] = round(float(s.max()), 1)
        det.append(fila)
    res = {}
    for v, s in sims.items():
        ok = s.notna() & q.notna()
        res[v] = {"NSE": round(float(nse(q[ok].values, s[ok].values)), 3),
                  "NSE_sense_4_10": round(float(nse(q[ok & ~q.index.to_series().between("2026-10-03", "2026-10-07")].values,
                                                    s[ok & ~q.index.to_series().between("2026-10-03", "2026-10-07")].values)), 3)}
    return res, det


def ridaura():
    cfg = json.loads(Path("config/ridaura.json").read_text()); par = parametres("ridaura")
    mod = json.loads(Path("config/ridaura_model_pluviometre.json").read_text())
    qt = cfg["cabals_T"]["aforador"]
    ob = pd.read_csv("data/ridaura/humitat/observacions.csv", parse_dates=["temps"]).set_index("temps")
    ob.index = ob.index.tz_localize(None) if ob.index.tz is not None else ob.index
    aca5 = (ob.pluja_mmh * DT / 60).asfreq(f"{DT}min")
    q = ob.cabal_m3s.asfreq(f"{DT}min")
    xh = pluja_xema("ridaura", par["estacions_xema"])
    ph = pluja_conca("ridaura", par["estacions_xema"])          # com en directe: pluviòmetre ACA on n'hi ha
    ac = aca5.dropna(); p5 = a5(xh); p5 = p5.reindex(p5.index.union(ac.index)).fillna(0.0); p5[ac.index] = ac
    fix = mod["Ia_mm"]
    ias = {"fixa": lambda t0: fix, "diposit": lambda t0: ia(t0, ph, par)}
    mods = {"fixa": mod, "diposit": par["model"]}
    out = {"Ia_fixa_mm": fix, "smax_mm": par["smax_mm"], "model_diposit": par["model"]}
    # a) 4/10
    t0, fi = pd.Timestamp("2026-10-03T12:00"), pd.Timestamp("2026-10-05T00:00")
    a = {}
    for v, f in ias.items():
        x = f(t0); s = simula_ep(p5, t0, fi, x, mods[v]).reindex(q[t0:fi].index)
        a[v] = {"ia_mm": round(x, 1), "NSE": round(float(nse(q[t0:fi].values, s.values)), 3), "pic_sim": round(float(s.max()), 1)}
    out["4_10_2026"] = {"pic_obs": round(float(q[t0:fi].max()), 1), **a}
    # b) gener 2026 (només XEMA)
    px = a5(xh)
    t0 = inici_episodi(px, pd.Timestamp("2026-01-19T00:00")); fi = pd.Timestamp("2026-01-22T00:00")
    out["gener_2026"] = {"inici": f"{t0:%Y-%m-%d %H:%M}", "pluja_xema_mm": round(float(px[t0:fi].sum()), 1), "pic_obs": 28.6,
                         **{v: {"ia_mm": round(f(t0)), "pic_sim": round(float(simula_ep(px, t0, fi, f(t0), mods[v]).max()), 1)} for v, f in ias.items()}}
    # c) continu
    out["continu_ago_nov_2026"], out["episodis_2026"] = continu(p5, q, mods, ias, mod["area_km2"])
    # d) catàleg XEMA 2009-2026
    cat = pd.read_csv("data/ridaura/xema/episodis_pluja.csv"); cat = cat[cat.estacio == "DO"]
    files, n = [], {v: 0 for v in ias}
    for r in cat.itertuples():
        t0 = inici_episodi(px, pd.Timestamp(r.inici_utc)); fi = pd.Timestamp(r.fi_utc) + pd.Timedelta(hours=24)
        fila = {"inici": r.inici_utc[:10], "pluja_mm": round(float(px[t0:fi].sum()))}
        for v, f in ias.items():
            x = f(t0); qm = float(simula_ep(px, t0, fi, x, mods[v]).max())
            fila[f"ia_{v}"] = round(x); fila[f"q_{v}"] = round(qm, 1)
            if qm >= qt["QMCO"] and r.inici_utc[:10] not in ("2026-10-03", "2026-10-04"): n[v] += 1
        files.append(fila)
    out["cataleg_xema"] = {"episodis": len(files), "sobre_QMCO_fora_4_10": n,
                           "mes_humits": sorted(files, key=lambda f: f["ia_diposit"])[:10],
                           "mes_cabal_diposit": sorted(files, key=lambda f: -f["q_diposit"])[:10]}
    return out


def congost():
    par = parametres("granollers")
    ob = pd.read_csv("data/granollers/humitat/aforador.csv", parse_dates=["temps"]).set_index("temps").cabal_m3s
    q = ob[~ob.index.duplicated()].asfreq(f"{DT}min")
    xh = pluja_xema("granollers", par["estacions_xema"])
    p5 = a5(xh)
    area = 100.0      # nominal: C és lliure (sense sostre d'1), i només compta el producte C·àrea
    t0c, fic = pd.Timestamp(par["calibratge"]["t0"]), pd.Timestamp("2026-10-06T00:00")
    ias = {"fixa": lambda t0: par["calibratge"]["ia_mm"], "diposit": lambda t0: ia(t0, xh, par)}
    mods = {}
    for v, f in ias.items():                     # model ràpid ajustat al 3-5/10 amb l'Ia de cada variant (C, K, m lliures)
        x = f(t0c); pp = p5[t0c:fic].values; oo = q.reindex(p5[t0c:fic].index).values
        qb = float(np.nanmin(oo[:12])) if np.isfinite(oo[:12]).any() else 0.0
        cost = lambda y: 1e9 if not (0 < y[0] and y[1] > 0 and 0.15 < y[2] < 1.5) else -nse(oo, simula(pp, area, qb, x, y[0], y[1], y[2], 2))
        r = min((minimize(cost, [c, k, m], method="Nelder-Mead", options={"maxiter": 800}) for c in (0.3, 0.7) for k in (1e3, 1e5) for m in (0.4, 1.0)), key=lambda r: r.fun)
        mods[v] = {"area_km2": area, "C": float(r.x[0]), "K": float(r.x[1]), "m": float(r.x[2]), "n": 2, "NSE_calibratge": round(-r.fun, 3)}
    res, det = continu(p5, q, mods, ias, area)
    return {"Ia_fixa_mm": par["calibratge"]["ia_mm"], "smax_mm": par["smax_mm"], "area_km2": area,
            "models": {v: {k: (round(x, 3) if isinstance(x, float) else x) for k, x in m.items()} for v, m in mods.items()},
            "continu_jul_nov_2026": res, "episodis": det}


def valida():
    for nom, f in (("ridaura", ridaura), ("granollers", congost)):
        try:
            r = f()
        except FileNotFoundError as e:
            print(nom, "sense dades:", e); continue
        Path(f"data/{nom}/humitat/validacio.json").write_text(json.dumps(r, ensure_ascii=False, indent=1, default=str) + "\n")
        print(f"== {nom}"); print(json.dumps(r, ensure_ascii=False, indent=1, default=str)[:6000])
    return 0
