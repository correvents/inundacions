#!/usr/bin/env python3
"""Pluja forta vs actuacions dels Bombers per inundació d'aigua, dia a dia, en una conca (validació sense premsa).

Per a cada episode de pluja (estacions XEMA del domini; un episodi s'acaba després de 12 h seques) compta les actuacions
«Inundació d'aigua» dels municipis de la conca des del dia que comença fins l'endemà que acaba (hora local). Com que
aquesta alarma també inclou fuites i canonades, cada municipi té un soroll de fons (actuacions/dia els dies secs): un
episodi «amb inundacions» és el que en té clarament més del que tocaria pel soroll (Poisson, p < 0,01).
Llavors busca quina pluja (màxim en 30 min, en 1 h, total) separa millor els episodis amb i sense actuacions, igual que
la regla de Platja d'Aro (≥15 mm en 30 min o ≥70 mm), i en dona els encerts i les falses alarmes.

Entrades: data/<conca>/xema/<codi>.csv i estacions.json (pluja_xema.py), data/<conca>/model/bombers_aigua_dia.csv
(bombers.py). Sortida: data/<conca>/model/bombers_vs_pluja.csv (un episodi per fila, sense adreces) i el resum per pantalla.

Ús:  python3 scripts/valida_bombers.py --config config/granollers.json [--radi-km 8]"""
import argparse, json, math
from pathlib import Path
import numpy as np, pandas as pd


def poisson_sf(k, mu):
    """P(X ≥ k) per a una Poisson de mitjana mu."""
    if k <= 0:
        return 1.0
    return float(1 - sum(math.exp(-mu) * mu ** i / math.factorial(i) for i in range(int(k))))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--radi-km", type=float, default=8, help="estacions a menys d'aquesta distància del centre de la conca")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text()); d = Path("data") / cfg["conca"]
    est = [e for e in json.loads((d / "xema" / "estacions.json").read_text()) if e["km"] <= a.radi_km and (d / "xema" / f"{e['codi_estacio']}.csv").exists()]
    series = {e["codi_estacio"]: pd.read_csv(d / "xema" / f"{e['codi_estacio']}.csv", parse_dates=["temps"]).set_index("temps").mm for e in est}
    pl = pd.DataFrame(series).sort_index()
    pl = pl[pl.notna().sum(axis=1) > 0]
    print("estacions:", {c: (str(s.index[0].date()), len(s)) for c, s in series.items()})
    # el màxim entre estacions: el xàfec que cau a una part de la ciutat ja pot inundar
    b = pd.read_csv(d / "model" / "bombers_aigua_dia.csv", parse_dates=["dia"]).groupby("dia").actuacions.sum()
    fi_b = b.index.max()
    pl = pl[(pl.index >= max(pd.Timestamp("2009-01-01"), b.index.min())) & (pl.index <= fi_b)]
    mitj = pl.mean(axis=1)
    # episodis sobre la mitjana de les estacions
    plou = mitj[mitj > 0]
    eps, cur = [], None
    for t, v in plou.items():
        if cur and t - cur["fi"] > pd.Timedelta(hours=12):
            eps.append(cur); cur = None
        cur = cur or {"inici": t, "fi": t}
        cur["fi"] = t
    if cur: eps.append(cur)
    # soroll de fons: actuacions els dies que no plou (ni el dia abans)
    dies_pluja = set()
    for e in eps:
        for k in range(-1, (e["fi"] - e["inici"]).days + 3):
            dies_pluja.add((e["inici"] + pd.Timedelta(days=k)).normalize())
    tots = pd.date_range(pl.index.min().normalize(), fi_b, freq="1D")
    secs = [t for t in tots if t not in dies_pluja]
    fons = b.reindex(secs).fillna(0).mean()
    print(f"soroll de fons: {fons:.2f} actuacions/dia sec ({len(secs)} dies)")
    files = []
    for e in eps:
        s = pl[(pl.index >= e["inici"]) & (pl.index <= e["fi"])]
        tot = float(s.sum().max())
        if tot < 10:
            continue
        m30 = float(s.max().max()); m1h = float(s.rolling(2, min_periods=1).sum().max().max())
        loc0 = e["inici"].tz_localize("UTC").tz_convert("Europe/Madrid").normalize().tz_localize(None)
        loc1 = e["fi"].tz_localize("UTC").tz_convert("Europe/Madrid").normalize().tz_localize(None) + pd.Timedelta(days=1)
        dies = pd.date_range(loc0, loc1, freq="1D")
        n = int(b.reindex(dies).fillna(0).sum())
        mu = fons * len(dies)
        p = poisson_sf(n, mu)
        files.append({"inici_utc": e["inici"].isoformat(), "fi_utc": e["fi"].isoformat(), "mm_max_estacio": round(tot, 1),
                      "max30min": round(m30, 1), "max1h": round(m1h, 1), "actuacions": n, "esperades_soroll": round(mu, 2),
                      "inundacio": bool(p < 0.01 and n >= 2)})
    df = pd.DataFrame(files)
    df.to_csv(d / "model" / "bombers_vs_pluja.csv", index=False)
    si, no = df[df.inundacio], df[~df.inundacio]
    print(f"{len(df)} episodis de ≥10 mm; {len(si)} amb inundacions (més actuacions que el soroll, p < 0,01)")
    print(si.sort_values("actuacions", ascending=False).head(25).to_string(index=False))
    # quin llindar separa millor (índex de Youden) per a cada variable
    for var in ("max30min", "max1h", "mm_max_estacio"):
        best = max(((((si[var] >= u).mean() - (no[var] >= u).mean()), u) for u in sorted(df[var].unique())), default=(0, 0))
        u = best[1]
        print(f"{var}: millor llindar {u} → encerta {(si[var] >= u).sum()}/{len(si)}, avisa en {(no[var] >= u).sum()}/{len(no)} sense inundació")
    # regla de dues condicions com a Platja d'Aro: xàfec curt O molta pluja acumulada
    millor = (-1, None)
    for u1 in sorted(df.max30min.unique()):
        for u2 in sorted(df.mm_max_estacio.unique()):
            av = (df.max30min >= u1) | (df.mm_max_estacio >= u2)
            j = av[df.inundacio].mean() - av[~df.inundacio].mean()
            if j > millor[0]:
                millor = (j, (u1, u2))
    u1, u2 = millor[1]; av = (df.max30min >= u1) | (df.mm_max_estacio >= u2)
    print(f"regla: ≥{u1} mm en 30 min o ≥{u2} mm l'episodi → encerta {av[df.inundacio].sum()}/{len(si)}, "
          f"avisa en {av[~df.inundacio].sum()}/{len(no)} sense inundació")
    return 0


if __name__ == "__main__":
    main()
