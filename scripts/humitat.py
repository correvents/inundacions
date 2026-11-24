#!/usr/bin/env python3
"""Absorció inicial del terreny (Ia) segons la pluja dels 60 dies anteriors: dipòsit d'humitat del sòl.

El terreny es tracta com un dipòsit amb un dèficit D (mm d'aigua que encara pot empassar-se). Cada hora:
    D ← D + ET0 · (1 − D/Smax)          l'evapotranspiració l'asseca, cada cop menys com més sec és (Thornthwaite-Mather)
    D ← max(D − pluja, 0)                la pluja l'omple
i l'absorció inicial d'un episodi és el dèficit en començar-lo: Ia = D(t0).

  - ET0: evapotranspiració de referència FAO-56 de la reanàlisi ERA5 (Open-Meteo, humitat_dades.py), en climatologia:
    la mitjana de cada dia de l'any 2009-2025 al punt de l'aforador. En directe no cal cap dada nova.
  - Pluja: la dels 60 dies anteriors (pluviòmetre ACA / radar en directe; estacions XEMA per als episodis històrics).
  - Estat 60 dies abans: el dèficit típic (mediana) d'aquell dia de l'any, del dipòsit fet córrer amb tota la sèrie XEMA
    des del 2009 (precalculat).
  - Smax és l'únic paràmetre i es calibra amb el cabal observat: el que fa que el dipòsit doni, en començar l'episodi
    de calibratge, l'Ia que el model de cabal hi va ajustar (Ridaura 4/10/2026: 159 mm; Congost 3-4/10/2026: 65 mm).
    Res més no es toca a mà.

Ús:  python3 scripts/humitat.py calibra --config config/ridaura.json     → config/<conca>_humitat.json
     python3 scripts/humitat.py ia --config config/ridaura.json --t0 2026-10-03T12:00
     python3 scripts/humitat.py valida                                    → informe de validació (stdout i JSON)"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

DIES = 60


# ---------- dades ----------
def pluja_xema(conca, estacions):
    """Pluja horària (mm) mitjana de les estacions XEMA disponibles a cada hora."""
    ss = []
    for c in estacions:
        f = Path("data") / conca / "xema" / f"{c}.csv"
        if f.exists():
            s = pd.read_csv(f, parse_dates=["temps"]).set_index("temps").mm.astype(float)
            ss.append(s[~s.index.duplicated()].resample("1h").sum(min_count=1))
    return pd.concat(ss, axis=1).mean(axis=1).fillna(0.0)


def pluja_conca(conca, estacions):
    """XEMA, substituïda pel pluviòmetre de l'ACA a l'aforador on n'hi ha (observacions.csv), que és el que llegeix la previsió."""
    p = pluja_xema(conca, estacions)
    f = Path("data") / conca / "humitat" / "observacions.csv"
    if f.exists():
        o = pd.read_csv(f, parse_dates=["temps"]).set_index("temps").pluja_mmh.dropna() / 12
        o.index = o.index.tz_localize(None) if o.index.tz is not None else o.index
        oh = o.resample("1h").sum(min_count=1).dropna()
        p = p.reindex(p.index.union(oh.index)).fillna(0.0); p[oh.index] = oh
    return p


def et0_clima(conca):
    """ET0 mitjana (mm/dia) per a cada dia de l'any (1..366), suavitzada amb una finestra circular de 31 dies."""
    e = pd.read_csv(Path("data") / conca / "humitat" / "et0.csv", parse_dates=["dia"])
    e = e[e.dia < "2026-01-01"]
    m = e.groupby(e.dia.dt.dayofyear).et0_mm.mean().reindex(range(1, 367)).interpolate().bfill().ffill().values
    k = 15
    ext = np.r_[m[-k:], m, m[:k]]
    return np.convolve(ext, np.ones(2 * k + 1) / (2 * k + 1), mode="valid")


# ---------- dipòsit ----------
def corre(pluja_h, et0_h, smax, d0):
    """Dèficit (mm) al final de cada hora."""
    d, out = float(d0), np.empty(len(pluja_h))
    for i, (p, e) in enumerate(zip(pluja_h, et0_h)):
        d += e * (1 - d / smax)
        d = max(d - p, 0.0)
        out[i] = d
    return out


def et0_hores(index, et0c):
    return et0c[np.minimum(index.dayofyear.values, 366) - 1] / 24.0


def clima_deficit(pluja_h, et0c, smax):
    """Mediana del dèficit de cada dia de l'any amb el dipòsit fet córrer amb tota la sèrie (el 1r any s'escalfa)."""
    d = corre(pluja_h.values, et0_hores(pluja_h.index, et0c), smax, smax)
    s = pd.Series(d, index=pluja_h.index)
    s = s[s.index >= s.index[0] + pd.Timedelta(days=365)]
    m = s.groupby(s.index.dayofyear).median().reindex(range(1, 367)).interpolate().bfill().ffill().values
    return m


def ia(t0, pluja_h, par):
    """Absorció inicial (mm) en començar un episode a t0 (UTC), amb la pluja horària dels 60 dies anteriors.
    pluja_h: sèrie horària (pd.Series indexada per hora); les hores sense dada compten com a seques."""
    t0 = pd.Timestamp(t0).tz_localize(None) if pd.Timestamp(t0).tzinfo else pd.Timestamp(t0)
    t0 = t0.floor("1h")                    # la sèrie és horària: l'hora en curs encara no compta
    ini = t0 - pd.Timedelta(days=DIES)
    idx = pd.date_range(ini, t0 - pd.Timedelta(hours=1), freq="1h")
    p = pluja_h.reindex(idx).fillna(0.0).values
    et0c, dc = np.array(par["et0_clima_mm_dia"]), np.array(par["deficit_clima_mm"])
    d0 = dc[min(ini.dayofyear, 366) - 1]
    return float(corre(p, et0_hores(idx, et0c), par["smax_mm"], d0)[-1]) if len(p) else float(d0)


def parametres(conca):
    f = Path("config") / f"{conca}_humitat.json"
    return json.loads(f.read_text()) if f.exists() else None


# ---------- calibratge ----------
def calibra(cfg):
    conca, h = cfg["conca"], cfg["humitat"]
    pluja_h = pluja_conca(conca, h["estacions"])
    et0c = et0_clima(conca)
    t0, ia_obj = pd.Timestamp(h["calibratge"]["t0"]), h["calibratge"]["ia_mm"]

    def fa(smax):
        par = {"smax_mm": smax, "et0_clima_mm_dia": et0c, "deficit_clima_mm": clima_deficit(pluja_h, et0c, smax)}
        return ia(t0, pluja_h, par), par
    lo, hi = ia_obj, 5 * ia_obj
    for _ in range(30):                                   # bisecció: Ia(t0) creix amb Smax
        mid = (lo + hi) / 2
        if fa(mid)[0] < ia_obj: lo = mid
        else: hi = mid
    v, par = fa((lo + hi) / 2)
    out = {"descripcio": "Dipòsit d'humitat del sòl (scripts/humitat.py): Ia = dèficit en començar l'episodi.",
           "smax_mm": round((lo + hi) / 2, 1), "estacions_xema": h["estacions"],
           "calibratge": {**h["calibratge"], "ia_dipòsit_mm": round(v, 1)},
           "et0_clima_mm_dia": [round(x, 3) for x in et0c], "deficit_clima_mm": [round(x, 1) for x in par["deficit_clima_mm"]]}
    Path(f"config/{conca}_humitat.json").write_text(json.dumps(out, ensure_ascii=False) + "\n")
    print(f"{conca}: Smax {out['smax_mm']} mm → Ia({t0:%Y-%m-%d %H:%M}) = {v:.1f} mm (objectiu {ia_obj}); "
          f"ET0 {et0c.sum():.0f} mm/any; dèficit típic 1 gen {par['deficit_clima_mm'][0]:.0f}, 1 jul {par['deficit_clima_mm'][181]:.0f}, "
          f"1 oct {par['deficit_clima_mm'][273]:.0f} mm")
    return out


# ---------- calibratge conjunt amb el cabal (conques amb model de previsió) ----------
def calibra_cabal(cfg):
    """Ridaura: Smax, la pèrdua contínua phi i la propagació (K, m) del model de cabal de la previsió, junts.

    Amb un sol episodi (4/10, sòl sec) l'absorció inicial i la pèrdua contínua es confonen: tot el que el terreny s'empassa
    es pot posar a Ia (159 mm, phi 0) o repartir-ho. Un episodi amb el sòl moll les separa: el 17-21/1/2026 (pic 28,6 m³/s
    a l'aforador). Per a cada phi (0-20 mm/h) s'ajusten Ia(4/10), K i m al hidrograma del 4/10 (NSE, pluviòmetre ACA), Smax
    surt d'Ia(4/10), i es tria la phi que dona el pic de gener (pluja XEMA, Ia del dipòsit). C = 1.
    Escriu config/<conca>_humitat.json amb el model sencer (la previsió el fa servir en lloc de <conca>_model_pluviometre.json)."""
    from scipy.optimize import minimize
    from model_cabal import simula, nse, DT
    conca, h = cfg["conca"], cfg["humitat"]; cal = h["calibratge"]
    area = cfg["cabals_T"]["aforador"]["area_km2"]
    ph = pluja_conca(conca, h["estacions"]); et0c = et0_clima(conca)
    o = pd.read_csv(Path("data") / conca / "humitat" / "observacions.csv", parse_dates=["temps"]).set_index("temps")
    o.index = o.index.tz_localize(None) if o.index.tz is not None else o.index
    t0, fi = pd.Timestamp(cal["t0"]), pd.Timestamp(cal["fi"])
    w = o[t0:fi]; p1 = (w.pluja_mmh.fillna(0) * DT / 60).values; q1 = w.cabal_m3s.values
    ind = cal["independent"]
    xh = pluja_xema(conca, h["estacions"])
    idx = pd.date_range(ind["t0"], ind["fi"], freq=f"{DT}min")
    p2 = (xh.reindex(idx.floor("1h")).fillna(0).values) / (60 // DT)
    files = []
    for phi in np.arange(0, 20.1, 1.0):
        cost = lambda y: 1e9 if not (y[0] >= 0 and y[1] > 0 and 0.15 < y[2] < 1.5) else -nse(q1, simula(p1, area, 0, y[0], 1, y[1], y[2], 2, phi))
        r = min((minimize(cost, [a, k, m], method="Nelder-Mead", options={"maxiter": 800})
                 for a in (40, 100, 150) for k in (200, 600) for m in (0.8, 1.3)), key=lambda r: r.fun)
        ia1, k, m = (float(v) for v in r.x)
        lo, hi = max(ia1, 1.0), 6 * max(ia1, 30.0)
        for _ in range(25):                                  # Smax tal que el dipòsit doni ia1 el 4/10
            mid = (lo + hi) / 2
            par = {"smax_mm": mid, "et0_clima_mm_dia": et0c, "deficit_clima_mm": clima_deficit(ph, et0c, mid)}
            if ia(t0, ph, par) < ia1: lo = mid
            else: hi = mid
        ia2 = ia(pd.Timestamp(ind["t0"]), ph, par)
        q2 = float(simula(p2, area, 0, ia2, 1, k, m, 2, phi).max())
        files.append({"phi": float(phi), "ia_4_10": ia1, "K": k, "m": m, "NSE_4_10": -r.fun, "smax": par["smax_mm"],
                      "ia_independent": ia2, "pic_independent": q2, "par": par})
        print(f"phi {phi:4.1f}: Ia(4/10) {ia1:5.1f} K {k:6.1f} m {m:.3f} NSE {-r.fun:.3f} → Smax {mid:5.1f}, "
              f"Ia(gener) {ia2:5.1f}, pic gener {q2:6.1f} (obs {ind['pic_m3s']})")
    b = min(files, key=lambda f: abs(np.log(max(f["pic_independent"], 0.01) / ind["pic_m3s"])) - 0.5 * f["NSE_4_10"])
    out = {"descripcio": "Dipòsit d'humitat del sòl (scripts/humitat.py): Ia = dèficit en començar l'episodi; després pèrdua "
                         "contínua phi; propagació amb 2 embassaments no lineals (model_cabal.simula, C = 1).",
           "smax_mm": round(b["smax"], 1), "estacions_xema": h["estacions"],
           "model": {"phi_mm_h": b["phi"], "K": round(b["K"], 1), "m": round(b["m"], 3), "n": 2, "C": 1.0, "area_km2": area},
           "calibratge": {**cal, "ia_4_10_mm": round(b["ia_4_10"], 1), "NSE_4_10": round(b["NSE_4_10"], 3),
                          "ia_independent_mm": round(b["ia_independent"], 1), "pic_independent_sim": round(b["pic_independent"], 1),
                          "proves_phi": [{k: round(v, 3) for k, v in f.items() if k != "par"} for f in files]},
           "et0_clima_mm_dia": [round(x, 3) for x in et0c], "deficit_clima_mm": [round(x, 1) for x in b["par"]["deficit_clima_mm"]]}
    Path(f"config/{conca}_humitat.json").write_text(json.dumps(out, ensure_ascii=False) + "\n")
    print(f"tria: phi {b['phi']} mm/h, Smax {out['smax_mm']} mm")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("accio", choices=["calibra", "ia", "valida"])
    ap.add_argument("--config"); ap.add_argument("--t0")
    a = ap.parse_args()
    if a.accio == "valida":
        from humitat_valida import valida
        return valida()
    cfg = json.loads(Path(a.config).read_text())
    if a.accio == "calibra":
        calibra_cabal(cfg) if cfg["humitat"]["calibratge"].get("independent") else calibra(cfg)
    else:
        par = parametres(cfg["conca"])
        print(f"Ia = {ia(a.t0, pluja_conca(cfg['conca'], cfg['humitat']['estacions']), par):.1f} mm")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).parent))
    sys.exit(main())
