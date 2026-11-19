#!/usr/bin/env python3
"""Corregeix la pluja horària del radar de tot Catalunya (data/pluja_cat/AAAA-MM.bin.gz) amb els pluviòmetres (AAAA-MM_estacions.json).

Per a cada hora:
  1. Biaix mitjà del camp: F = Σ pluviòmetres / Σ radar als mateixos punts (entre 0,2 i 5; 1 si plou massa poc per saber-ho).
  2. Residus: la diferència pluviòmetre − radar×F a cada estació s'interpola amb pesos per distància (radi 30 km) i
     s'esvaeix lluny de les estacions, de manera que on no n'hi ha queda el radar×F.
Control de qualitat: es descarta una estació l'hora que marca > 150 mm/h, o tot el mes si el seu total és < 1/4 o > 4 vegades
el de les estacions veïnes (< 25 km) amb prou pluja (pluviòmetre embussat o sensor mal configurat).
Validació deixant cada estació fora (la correcció es calcula sense ella i es compara amb el que va mesurar).
Escriu AAAA-MM_corr.bin.gz (mateix format) i AAAA-MM_validacio.json, i marca el mes com a corregit a index.json.

Ús:  python3 scripts/correccio_cat.py --mes 2026-10"""
import argparse, gzip, json
from pathlib import Path
import numpy as np

RADI, ESC, NUGGET = 30.0, 4.0, 12.0     # km: radi d'influència, distància mínima efectiva, distància a què el residu ja pesa la meitat


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mes", required=True)
    ap.add_argument("--desfase", type=int, default=0, help="hores a desplaçar els pluviòmetres (per provar l'alineació)")
    a = ap.parse_args()
    d = Path("data/pluja_cat"); idx = json.loads((d / "index.json").read_text())
    F_, C_, P = idx["files"], idx["columnes"], idx["pas"]
    q = np.frombuffer(gzip.decompress((d / f"{a.mes}.bin.gz").read_bytes()), np.uint8).reshape(-1, F_, C_)
    H = q.shape[0]
    R = np.where(q < 255, (q / 20.0) ** 2, np.nan)
    # ecos fixos (muntanya, parcs eòlics): píxels amb un total mensual > 3 vegades la mediana del voltant (±5 px) i > 100 mm
    tot_r = np.nansum(R, 0); veins = np.zeros_like(tot_r)
    for i in range(F_):
        for j in range(C_):
            veins[i, j] = np.median(tot_r[max(0, i - 5):i + 6, max(0, j - 5):j + 6])
    eco = (tot_r > 3 * veins) & (tot_r > 100)
    R[:, eco] = np.nan
    print(f"Ecos fixos tapats: {eco.sum()} píxels de {eco.size}")
    est = json.loads((d / f"{a.mes}_estacions.json").read_text())
    lat0 = 41.7; kx = 111.32 * np.cos(np.radians(lat0)); ky = 110.57
    E = [e for e in est if idx["sud"] < e["lat"] < idx["nord"] and idx["oest"] < e["lon"] < idx["est"]]
    G = np.array([[np.nan if v is None else v for v in e["mm"]] for e in E], float)
    if a.desfase: G = np.roll(G, a.desfase, axis=1)
    G[G > 150] = np.nan
    fil = np.array([int((idx["nord"] - e["lat"]) / P) for e in E]); col = np.array([int((e["lon"] - idx["oest"]) / P) for e in E])
    X = np.array([e["lon"] * kx for e in E]); Y = np.array([e["lat"] * ky for e in E])
    D = np.hypot(X[:, None] - X[None], Y[:, None] - Y[None])
    # control de qualitat mensual amb els veïns
    tot = np.nansum(G, 1); cobert = np.isfinite(G).mean(1)
    dolent = np.zeros(len(E), bool)
    for i in range(len(E)):
        v = (D[i] < 25) & (D[i] > 0) & (cobert > 0.7)
        if v.sum() >= 2 and cobert[i] > 0.5:
            ref = np.median(tot[v] / np.maximum(cobert[v], 0.01)) * cobert[i]
            if ref > 20 and not (ref / 4 < tot[i] < ref * 4): dolent[i] = True
    G[dolent] = np.nan
    print(f"{len(E)} pluviòmetres dins la reixa; {dolent.sum()} descartats pel control de qualitat:",
          [E[i]["nom"] for i in np.where(dolent)[0]][:15])
    Rg = R[:, fil, col].T                                  # radar al píxel de cada estació [est, H]
    # reixa en km
    LON = idx["oest"] + (np.arange(C_) + 0.5) * P; LAT = idx["nord"] - (np.arange(F_) + 0.5) * P
    GX, GY = np.meshgrid(LON * kx, LAT * ky)
    Dg = np.hypot(GX.ravel()[None] - X[:, None], GY.ravel()[None] - Y[:, None])     # [est, píxels]
    Wg = np.where(Dg < RADI, 1 / (np.maximum(Dg, ESC) ** 2), 0.0)
    W = np.where(D < RADI, 1 / (np.maximum(D, ESC) ** 2), 0.0); np.fill_diagonal(W, 0)
    k0 = 1 / NUGGET ** 2
    out = np.full_like(R, np.nan); Fh = np.ones(H)
    loo_raw, loo_mfb, loo_cor, obs, quin = [], [], [], [], []
    for h in range(H):
        g, r = G[:, h], Rg[:, h]
        ok = np.isfinite(g) & np.isfinite(r)
        sg, sr = g[ok].sum(), r[ok].sum()
        F = float(np.clip(sg / sr, 0.25, 4)) if ok.sum() >= 5 and sr >= 2 and sg >= 2 else (1.0 if sr < 2 else float(np.clip(sg / max(sr, 0.5), 0.25, 4)))
        Fh[h] = F
        e = np.where(ok, g - r * F, 0.0)
        sostre = 1.5 * g[ok].max() if ok.sum() >= 5 else np.inf     # cap píxel corregit per sobre d'1,5 × el pluviòmetre més plujós de l'hora (si no ho deia ja el radar)
        w = Wg * ok[:, None]
        corr = (e[:, None] * w).sum(0) / (w.sum(0) + k0)
        out[h] = np.where(np.isnan(R[h]), np.nan, np.minimum(np.maximum(R[h] * F + corr.reshape(F_, C_), 0), np.maximum(R[h], sostre)))
        # deixant cada estació fora
        for i in np.where(ok)[0]:
            o = ok.copy(); o[i] = False
            sg_i, sr_i = g[o].sum(), r[o].sum()
            Fi = float(np.clip(sg_i / sr_i, 0.25, 4)) if o.sum() >= 5 and sr_i >= 2 and sg_i >= 2 else 1.0
            wi = W[i] * o
            ci = (np.where(o, g - r * Fi, 0) * wi).sum() / (wi.sum() + k0)
            loo_raw.append(r[i]); loo_mfb.append(r[i] * Fi); loo_cor.append(min(max(r[i] * Fi + ci, 0), max(r[i], 1.5 * g[o].max() if o.sum() >= 5 else np.inf))); obs.append(g[i]); quin.append((i, h // 24))
    obs, raw, mfb, cor = map(np.array, (obs, loo_raw, loo_mfb, loo_cor))
    plou = (obs > 0.1) | (raw > 0.1)
    def m(p):
        x, y = p[plou], obs[plou]
        return {"biaix_total": round(float(x.sum() / y.sum()), 2), "rmse_mm_h": round(float(np.sqrt(((x - y) ** 2).mean())), 2),
                "correlacio": round(float(np.corrcoef(x, y)[0, 1]), 2)}
    val = {"hores_estacio": int(plou.sum()), "radar_sol": m(raw), "biaix_mitja": m(mfb), "corregit": m(cor)}
    # també per dies (pluja diària de cada estació), que és el que compta per a l'acumulat d'un episodi
    from collections import defaultdict
    dia = defaultdict(lambda: np.zeros(4))
    for k, (i, dd) in enumerate(quin): dia[(i, dd)] += (obs[k], raw[k], mfb[k], cor[k])
    D4 = np.array([v for v in dia.values() if v[0] >= 1 or v[1] >= 1])
    def md(c):
        x, y = D4[:, c], D4[:, 0]
        return {"biaix_total": round(float(x.sum() / y.sum()), 2), "rmse_mm_dia": round(float(np.sqrt(((x - y) ** 2).mean())), 1),
                "correlacio": round(float(np.corrcoef(x, y)[0, 1]), 2)}
    val["diari"] = {"dies_estacio": len(D4), "radar_sol": md(1), "biaix_mitja": md(2), "corregit": md(3)}
    # acumulat mensual per estació (validació creuada)
    print("Validació deixant cada estació fora (hores amb pluja):", json.dumps(val))
    enc = np.where(np.isnan(out), 255, np.clip(np.round(np.sqrt(np.nan_to_num(out)) * 20), 0, 254)).astype(np.uint8)
    (d / f"{a.mes}_corr.bin.gz").write_bytes(gzip.compress(enc.tobytes(), 9))
    # estacions per a la web: només les útils, amb el seu control de qualitat
    web = [{"id": e["id"], "font": e["font"], "nom": e["nom"], "lat": round(e["lat"], 4), "lon": round(e["lon"], 4), "descartada": bool(dolent[i]),
            "mm": [None if not np.isfinite(v) else round(float(v), 1) for v in np.array([np.nan if x is None else x for x in e["mm"]], float)]}
           for i, e in enumerate(E)]
    (d / f"{a.mes}_estacions.json").write_text(json.dumps(web, ensure_ascii=False, separators=(",", ":")))
    (d / f"{a.mes}_validacio.json").write_text(json.dumps({**val, "factor_hora": [round(x, 2) for x in Fh],
        "estacions": len(E), "descartades": int(dolent.sum())}, ensure_ascii=False))
    idx["mesos"][a.mes].update({"corregit": True, "estacions": len(E) - int(dolent.sum()), "validacio": val})
    (d / "index.json").write_text(json.dumps(idx, ensure_ascii=False, indent=1))
    t = np.nansum(out, 0); i = np.unravel_index(np.nanargmax(t), t.shape)
    print(f"Màxim mensual corregit {t[i]:.0f} mm a {LAT[i[0]]:.2f}, {LON[i[1]]:.2f}")


if __name__ == "__main__":
    main()
