#!/usr/bin/env python3
"""Model pluja → cabal de la Ridaura a l'aforador (capa 2 del pla), calibrat amb un episodi observat.

Model agregat a pas de 5 minuts, pensat per a crescudes ràpides de conques petites:
  - pèrdues amb el llindar d'escorrentia P0 (mètode SCS, el que fa servir l'ACA):
        E acumulada = (P - P0)² / (P + 4·P0)  si P > P0
  - la pluja neta es propaga amb una cascada de Nash (n dipòsits lineals de constant k),
  - cabal de base constant (el de l'inici de l'episodi).
Paràmetres calibrats: P0 (mm), k (h) i n (1–5). El factor del radar f es fixa (per defecte 1: el 3-4/10/2026
el radar va donar 195 mm a la conca i el pluviòmetre de l'ACA 198 mm); amb un sol episodi P0 i f no es poden separar.

Ús:  python3 scripts/model_cabal.py --config config/ridaura.json --des 2026-10-03T00:00 --fins 2026-10-05T00:00
     (llegeix Supabase; amb --csv llegeix un CSV fet amb --exporta)
Desa els paràmetres a config/<conca>_model_<font>.json i la gràfica a data/<conca>/calibratge_<font>.png."""
import argparse, csv, json, sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

DT = 5          # minuts


def simula(pluja, area_km2, q_base, p0, k_h, f=1.0, n=3):
    """Cabal (m³/s) a cada pas a partir de la pluja mitjana de la conca (mm per pas)."""
    p = np.cumsum(pluja * f)
    e = np.where(p > p0, (p - p0) ** 2 / (p + 4 * p0), 0.0)
    neta = np.diff(e, prepend=0.0)                                   # mm per pas
    entrada = neta * area_km2 * 1000 / (DT * 60)                     # m³/s
    a = DT / 60 / max(k_h, 1e-3)                                     # fracció que surt de cada dipòsit per pas
    s = np.zeros(n); q = np.empty_like(entrada)
    for i, x in enumerate(entrada):
        flux = x
        for j in range(n):
            s[j] += flux; flux = s[j] * a / (1 + a); s[j] -= flux   # dipòsit lineal (implícit, estable)
        q[i] = flux
    return q + q_base


def nse(obs, sim):
    ok = ~np.isnan(obs)
    return 1 - np.sum((obs[ok] - sim[ok]) ** 2) / np.sum((obs[ok] - obs[ok].mean()) ** 2)


def calibra(pluja, obs, area_km2, f=1.0):
    q_base = max(float(np.nanmin(obs[: 12 * 6])), 0.0) if np.isfinite(obs[: 12 * 6]).any() else 0.0
    def cost(x, n):
        p0, k = x
        if p0 < 1 or k < 0.05:
            return 1e9
        sim = simula(pluja, area_km2, q_base, p0, k, f, n)
        ok = ~np.isnan(obs)
        # error quadràtic + penalització del pic (el que importa per a l'avís)
        return np.mean((obs[ok] - sim[ok]) ** 2) + 0.5 * (np.nanmax(obs) - sim.max()) ** 2 / 10
    millor, n_millor = None, None
    for n in range(1, 6):
        for p0 in (10, 40, 80):
            for k in (0.2, 0.6, 1.5):
                r = minimize(cost, [p0, k], args=(n,), method="Nelder-Mead", options={"maxiter": 2000, "xatol": 1e-3, "fatol": 1e-3})
                if millor is None or r.fun < millor.fun:
                    millor, n_millor = r, n
    p0, k = (float(v) for v in millor.x)
    return {"P0_mm": round(p0, 2), "k_h": round(k, 3), "n": n_millor, "factor_radar": f, "q_base_m3s": round(q_base, 3)}


def llegeix_supabase(cfg, ini, fi, font="radar"):
    from sb import llegeix
    def tot(taula, params):
        out, off = [], 0
        while True:
            r = llegeix(taula, {**params, "order": "temps", "offset": off, "limit": 1000})
            out += r; off += 1000
            if len(r) < 1000:
                return out
    rang = {"and": f"(temps.gte.{ini.isoformat()},temps.lte.{fi.isoformat()})"}
    if font == "radar":
        pl = tot("pluja_conca", {"select": "temps,mm", "conca": f"eq.{cfg['conca']}", "subconca": "eq.aforador", "font": "eq.radar", **rang})
    else:   # pluviòmetre de l'ACA a l'aforador (mm/h cada 5 min): un sol punt, sense distribució espacial
        pl = [{"temps": r["temps"], "mm": r["valor"] * DT / 60} for r in
              tot("observacions", {"select": "temps,valor", "conca": f"eq.{cfg['conca']}", "var": "eq.pluja_mmh", **rang})]
    qo = tot("observacions", {"select": "temps,valor", "conca": f"eq.{cfg['conca']}", "var": "eq.cabal_m3s", **rang})
    return ({datetime.fromisoformat(r["temps"]): r["mm"] for r in pl if r["mm"] is not None},
            {datetime.fromisoformat(r["temps"]): r["valor"] for r in qo})


def llegeix_csv(path, font):
    pl, qo = {}, {}
    for r in csv.DictReader(open(path)):
        t = datetime.fromisoformat(r["temps"]).replace(tzinfo=timezone.utc)
        if r.get(f"pluja_{font}_mm"): pl[t] = float(r[f"pluja_{font}_mm"])
        if r.get("cabal_m3s"): qo[t] = float(r["cabal_m3s"])
    return pl, qo


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--des", required=True)
    ap.add_argument("--fins", required=True)
    ap.add_argument("--csv")
    ap.add_argument("--factor-radar", type=float, default=1.0)
    ap.add_argument("--font", choices=["radar", "pluviometre"], default="radar")
    ap.add_argument("--exporta", help="només desa un CSV amb temps, pluja radar, pluja pluviòmetre i cabal, i acaba")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    ini = datetime.fromisoformat(a.des).replace(tzinfo=timezone.utc)
    fi = datetime.fromisoformat(a.fins).replace(tzinfo=timezone.utc)
    temps = [ini + timedelta(minutes=DT * i) for i in range(int((fi - ini).total_seconds() // (DT * 60)) + 1)]
    if a.exporta:
        rad, qo = llegeix_supabase(cfg, ini, fi, "radar"); plu, _ = llegeix_supabase(cfg, ini, fi, "pluviometre")
        with open(a.exporta, "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["temps", "pluja_radar_mm", "pluja_pluviometre_mm", "cabal_m3s"])
            for t in temps:
                w.writerow([t.strftime("%Y-%m-%dT%H:%M"), rad.get(t, ""), plu.get(t, ""), qo.get(t, "")])
        return 0
    pl, qo = llegeix_csv(a.csv, a.font) if a.csv else llegeix_supabase(cfg, ini, fi, a.font)
    pluja = np.array([pl.get(t, 0.0) for t in temps])
    obs = np.array([qo.get(t, np.nan) for t in temps])
    print(f"{len(pl)} passos de radar ({pluja.sum():.1f} mm), {np.isfinite(obs).sum()} de cabal (pic {np.nanmax(obs):.1f} m³/s)")
    area = cfg["cabals_T"]["aforador"]["area_km2"]
    par = calibra(pluja, obs, area, a.factor_radar)
    sim = simula(pluja, area, par["q_base_m3s"], par["P0_mm"], par["k_h"], par["factor_radar"], par["n"])
    ip_o, ip_s = int(np.nanargmax(obs)), int(np.argmax(sim))
    par.update({"font_pluja": a.font, "punt": "aforador", "area_km2": area, "episodi": [a.des, a.fins], "NSE": round(float(nse(obs, sim)), 3),
                "pic_obs_m3s": round(float(obs[ip_o]), 1), "pic_sim_m3s": round(float(sim[ip_s]), 1),
                "hora_pic_obs": temps[ip_o].isoformat(), "hora_pic_sim": temps[ip_s].isoformat(),
                "pluja_radar_mm": round(float(pluja.sum()), 1)})
    print(json.dumps(par, indent=2, ensure_ascii=False))
    Path(f"config/{cfg['conca']}_model_{a.font}.json").write_text(json.dumps(par, indent=2, ensure_ascii=False) + "\n")
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 4.5))
        ax.plot(temps, obs, "k", lw=1.5, label="Cabal observat (ACA)")
        ax.plot(temps, sim, "tab:red", lw=1.5, label=f"Model (NSE {par['NSE']:.2f})")
        for nom, c in (("Q10", "#c0392b"), ("Q50", "#e67e22"), ("Q100", "#b7950b"), ("Q500", "#2e86c1")):
            ax.axhline(cfg["cabals_T"]["aforador"][nom], color=c, ls=":", lw=1); ax.text(temps[0], cfg["cabals_T"]["aforador"][nom], nom, color=c, fontsize=8, va="bottom")
        ax.set_ylabel("m³/s"); ax.legend(loc="upper right")
        ax2 = ax.twinx(); ax2.bar(temps, pluja * par["factor_radar"] * 12, width=DT / 1440, color="tab:blue", alpha=0.4)
        ax2.set_ylim(ax2.get_ylim()[1] * 3, 0); ax2.set_ylabel("pluja (mm/h)")
        ax.set_title(f"Ridaura a Santa Cristina d'Aro · {a.font} · calibratge {a.des[:10]} – {a.fins[:10]} (UTC)")
        fig.tight_layout(); fig.savefig(f"data/{cfg['conca']}/calibratge_{a.font}.png", dpi=120)
    except ImportError:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
