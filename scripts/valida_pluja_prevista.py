#!/usr/bin/env python3
"""Validació de la pluja prevista (pluja_prevista.py) i de quanta antelació guanya l'avís (previsio.py --pluja-prevista).

A. Pluja: per a cada hora del període (--des/--fins, pas --pas min) es compara la pluja prevista sobre la conca a 1, 3 i
   6 h amb la mesurada (pluviòmetre de l'ACA a l'aforador; si falta, radar). Només s'avaluen les hores amb pluja propera
   (mesurada en ±6 h ≥ 0,5 mm o prevista per AROME ≥ 0,5 mm): la resta són hores seques sense res a avaluar. Taula de
   contingència per llindars (encerts, falses alarmes, no previstos) i biaix.
B. Antelació: per a cada episodi (--episodi INICI,FI) es fa la previsió cada 10 min, sense i amb pluja prevista, i es
   mira quan avisa per primer cop de crescuda ordinària i de T10, comparat amb quan ho mesura l'aforador.
   També amb la pluja que va caure de debò (previsió perfecta), que és el sostre del que es pot guanyar.
C. Falses alarmes de l'avís: totes les hores de la part A amb la previsió completa. Una alarma és falsa si l'avís puja per
   la pluja prevista i ni l'aforador (Ridaura) ni la previsió amb la pluja real (Fenals, sense aforador) hi arriben.

No desa res a Supabase. Escriu un informe Markdown (--sortida) i el detall en JSON al costat.

Ús:  python3 scripts/valida_pluja_prevista.py --config config/ridaura.json --des 2026-08-11T00:00 --fins 2026-11-11T18:00 \\
        --episodi 2026-10-04T02:00,2026-10-04T09:30 --sortida informes/pluja_prevista.md"""
import argparse, json, sys, time
from datetime import datetime, timedelta, timezone
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from sb import llegeix
import pluja_prevista as pp
import previsio

DT = 5
LLINDARS = {1: [5, 10, 20], 3: [10, 20, 40], 6: [20, 40, 70]}       # mm a 1, 3 i 6 h
FONTS = ["nowcast", "arome", "tot"]


def tot(taula, params):
    out, off = [], 0
    while True:
        r = llegeix(taula, {**params, "order": "temps", "offset": off, "limit": 1000})
        out += r; off += 1000
        if len(r) < 1000:
            return out


def utc(s):
    return datetime.fromisoformat(s).replace(tzinfo=timezone.utc) if "+" not in s[10:] else datetime.fromisoformat(s).astimezone(timezone.utc)


def series(cfg, ini, fi):
    """Pluja (mm per pas de 5 min, pluviòmetre i si falta radar) i cabal observat a l'aforador, indexats per temps."""
    rang = {"and": f"(temps.gt.{ini.isoformat()},temps.lte.{fi.isoformat()})"}
    c = cfg["conca"]
    pl = {utc(r["temps"]).replace(second=0): r["valor"] * DT / 60 for r in
          tot("observacions", {"select": "temps,valor", "conca": f"eq.{c}", "var": "eq.pluja_mmh", **rang}) if r["valor"] is not None}
    rd = {utc(r["temps"]).replace(second=0): r["mm"] for r in
          tot("pluja_conca", {"select": "temps,mm", "conca": f"eq.{c}", "subconca": "eq.aforador", "font": "eq.radar", **rang}) if r["mm"] is not None}
    q = {utc(r["temps"]).replace(second=0): r["valor"] for r in
         tot("observacions", {"select": "temps,valor", "conca": f"eq.{c}", "var": "eq.cabal_m3s", **rang}) if r["valor"] is not None}
    return pl, rd, q


def acumulat(pl, rd, t0, minuts):
    """mm mesurats entre t0 i t0+minuts (pas de 5 min, el pluviòmetre té preferència)."""
    s = 0.0
    for k in range(1, minuts // DT + 1):
        t = t0 + timedelta(minutes=DT * k)
        s += pl.get(t, rd.get(t, 0.0))
    return s


def contingencia(prev, obs, llindar):
    p, o = np.asarray(prev) >= llindar, np.asarray(obs) >= llindar
    a, b, c = int((p & o).sum()), int((p & ~o).sum()), int((~p & o).sum())
    return {"encerts": a, "falses": b, "no_previstos": c,
            "POD": round(a / (a + c), 2) if a + c else None, "FAR": round(b / (a + b), 2) if a + b else None}


def informe(cfg, res, sortida):
    L = [f"# Validació de la pluja prevista ({cfg['nom']})", "",
         f"Període {res['periode'][0]} - {res['periode'][1]} UTC, una previsió cada {res['pas_min']} min. Pluja mesurada: pluviòmetre de l'ACA a l'aforador "
         f"(si falta, radar). AROME: execució disponible a cada hora (≥ {pp.LATENCIA_H} h abans). Nowcast: radar OPERA (Marshall-Palmer).", ""]
    A = res["A"]
    L += [f"## A. Pluja prevista vs mesurada ({len(A)} hores amb pluja propera)", "",
          "| Font | Termini | Llindar | Encerts | Falses alarmes | No previstos | POD | FAR |", "|---|---|---|---|---|---|---|---|"]
    resum_biaix = []
    for font in FONTS:
        for h in (1, 3, 6):
            if font == "nowcast" and h > 1:
                continue
            ok = [f for f in A if font in f and h in f[font]]
            if not ok:
                continue
            prev = [f[font][h] for f in ok]; obs = [f["obs"][h] for f in ok]
            resum_biaix.append(f"{font} {h} h: prevista {sum(prev):.0f} mm, mesurada {sum(obs):.0f} mm (n = {len(ok)}), "
                               f"correlació {np.corrcoef(prev, obs)[0, 1]:.2f}" if np.std(prev) > 0 and np.std(obs) > 0 else f"{font} {h} h: n = {len(ok)}")
            for ll in LLINDARS[h]:
                c = contingencia(prev, obs, ll)
                L.append(f"| {font} | {h} h | ≥ {ll} mm | {c['encerts']} | {c['falses']} | {c['no_previstos']} | {c['POD']} | {c['FAR']} |")
    L += ["", "Totals (biaix) i correlació:", ""] + [f"- {x}" for x in resum_biaix] + [""]
    for ep in res["B"]:
        L += [f"## B. Antelació de l'avís, episodi {ep['episodi']}", "",
              f"Aforador: crescuda ordinària a les {ep['aforador'][1]}, T10 a les {ep['aforador'][2]} UTC.", "",
              "| Pluja prevista | Primer avís de crescuda | Antelació | Primer avís de T10 | Antelació |", "|---|---|---|---|---|"]
        def ant(x, y):
            return f"{int((utc(y) - utc(x)).total_seconds() // 60)} min" if x and y else "-"
        for font, v in ep["fonts"].items():
            nom = {"no": "cap (com ara)", "real": "la que va caure (previsió perfecta)"}.get(font, font)
            L.append(f"| {nom} | {v['crescuda'] or '-'} | {ant(v['crescuda'], ep['aforador'][1])} | {v['T10'] or '-'} | {ant(v['T10'], ep['aforador'][2])} |")
        L.append("")
    C = res["C"]
    if C:
        puja = [c for c in C if c["puja"]]
        falses = [c for c in puja if c["falsa"]]
        L += [f"## C. Falses alarmes de l'avís", "",
              f"Previsions avaluades: {len(C)} (hores amb ≥ 5 mm previstos o mesurats en 6 h). La pluja prevista fa pujar el nivell "
              f"en {len(puja)}; d'aquestes, {len(falses)} són falses (ni l'aforador ni la previsió amb la pluja real hi arriben).", ""]
        if puja:
            L += ["| Hora | Nivell amb prevista | Sense | Amb pluja real | Aforador | Falsa |", "|---|---|---|---|---|---|"]
            L += [f"| {c['t']} | {c['nivell']} | {c['sense_prevista']} | {c['amb_pluja_real']} | {c['aforador']} | {'sí' if c['falsa'] else 'no'} |" for c in puja]
            L.append("")
    out = Path(sortida); out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(L))
    out.with_suffix(".json").write_text(json.dumps(res, ensure_ascii=False, default=str))
    print("\n".join(L))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--des", required=True); ap.add_argument("--fins", required=True)
    ap.add_argument("--pas", type=int, default=60, help="minuts entre previsions a la part A i C")
    ap.add_argument("--episodi", action="append", default=[], help="INICI,FI (UTC) per a la part B")
    ap.add_argument("--sortida", default="informes/pluja_prevista.md")
    ap.add_argument("--sense-c", action="store_true", help="no facis la part C (és la més lenta)")
    ap.add_argument("--sense-a", action="store_true", help="només la part B (episodis)")
    ap.add_argument("--combina", nargs="+", help="ajunta els JSON de diverses execucions (per trossos de període) en un sol informe")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    if a.combina:
        rs = [json.loads(Path(f).read_text()) for f in a.combina]
        res = {"periode": [min(r["periode"][0] for r in rs if r["A"] or not r["B"]), max(r["periode"][1] for r in rs if r["A"] or not r["B"])],
               "pas_min": rs[0]["pas_min"], **{k: sorted(sum((r[k] for r in rs), []), key=lambda x: x.get("t", x.get("episodi", ""))) for k in "ABC"}}
        for f in res["A"]:                          # el JSON desa les claus numèriques com a text
            for k in ("obs", "arome", "nowcast", "tot"):
                if k in f:
                    f[k] = {(int(x) if x != "2" or k != "nowcast" else x): v for x, v in f[k].items()}
        informe(cfg, res, a.sortida)
        return 0
    qt = cfg["cabals_T"]["aforador"]
    ini, fi = utc(a.des), utc(a.fins)
    pl, rd, qobs = series(cfg, ini - timedelta(days=1), fi + timedelta(hours=7))
    print(f"Sèries: {len(pl)} pluviòmetre, {len(rd)} radar, {len(qobs)} cabal")
    cache = {}
    res = {"periode": [a.des, a.fins], "pas_min": a.pas, "A": [], "B": [], "C": []}

    # ---------- A. pluja ----------
    t = ini
    while t <= fi and not a.sense_a:
        obs = {h: acumulat(pl, rd, t, h * 60) for h in (1, 3, 6)}
        prop = acumulat(pl, rd, t - timedelta(hours=6), 12 * 60)
        fila = {"t": t.strftime("%Y-%m-%dT%H:%M"), "obs": obs}
        try:
            am, inf_am = pp.nwp(cfg, t, 360, directe=False)
        except Exception as e:
            am, inf_am = None, {"error": str(e)[:150]}
        am_6 = float(am.sum()) if am is not None else 0.0
        if prop >= 0.5 or am_6 >= 0.5:
            fila["arome_info"] = inf_am
            if am is not None:
                fila["arome"] = {h: float(am[: h * 12].sum()) for h in (1, 3, 6)}
            try:
                nc, inf_nc = pp.nowcast(cfg, t, pp.FUSIO_MIN, None, cache)
                fila["nowcast_info"] = inf_nc
                if nc is not None:
                    fila["nowcast"] = {h: float(nc[: h * 12].sum()) for h in (1,)} | {"2": float(nc.sum())}
            except Exception as e:
                fila["nowcast_info"] = {"error": f"{type(e).__name__} {e}"[:150]}
            if am is not None and "nowcast" in fila:
                w = np.clip((pp.FUSIO_MIN - (np.arange(72) + 0.5) * DT) / (pp.FUSIO_MIN - pp.NOWCAST_MIN), 0, 1)
                nc72 = np.r_[nc, np.zeros(72 - len(nc))]
                tt = w * nc72 + (1 - w) * am
                fila["tot"] = {h: float(tt[: h * 12].sum()) for h in (1, 3, 6)}
            res["A"].append(fila)
            print(fila["t"], "obs", {k: round(v, 1) for k, v in obs.items()}, "arome", {k: round(v, 1) for k, v in fila.get("arome", {}).items()},
                  "nowcast 1h", round(fila.get("nowcast", {}).get(1, -1), 1), flush=True)
        if len(cache) > 50:
            cache.clear()
        t += timedelta(minutes=a.pas)

    # ---------- B. antelació en episodis ----------
    def primer(llista, nivell):
        return next((r["t"] for r in llista if r["nivell"] >= nivell), None)
    for ep in a.episodi:
        e0, e1 = (utc(x) for x in ep.split(","))
        q_ep = sorted((tq, v) for tq, v in qobs.items() if e0 <= tq <= e1 + timedelta(hours=6))
        obs_n = {n: next((tq.strftime("%Y-%m-%dT%H:%M") for tq, v in q_ep if v >= qq), None) for n, qq in ((1, qt["QMCO"]), (2, qt["Q10"]))}
        sortida = {"episodi": ep, "aforador": obs_n, "fonts": {}}
        for font in ["no"] + FONTS + ["real"]:
            llista, t = [], e0
            while t <= e1:
                try:
                    if font == "real":           # previsió perfecta: la pluja que va caure després
                        orig = pp.prevista
                        pp.prevista = lambda cfg_, t0, minuts, *x, **k: (np.array([pl.get(t0 + timedelta(minutes=DT * (i + 1)), rd.get(t0 + timedelta(minutes=DT * (i + 1)), 0.0))
                                                                                  for i in range(minuts // DT)]), {"mm_3h": acumulat(pl, rd, t0, 180)})
                        try:
                            r = previsio.preveu(cfg, t, "tot", desar=False, verb=False)
                        finally:
                            pp.prevista = orig
                    else:
                        r = previsio.preveu(cfg, t, font, desar=False, cache=cache, verb=False)
                    llista.append({"t": t.strftime("%Y-%m-%dT%H:%M"), "nivell": r["nivell"], "q_max": {k: round(v, 1) for k, v in r["q_max"].items()},
                                   "prevista": r["prevista"]})
                except SystemExit as e:
                    llista.append({"t": t.strftime("%Y-%m-%dT%H:%M"), "nivell": -1, "error": str(e)})
                t += timedelta(minutes=10)
            sortida["fonts"][font] = {"crescuda": primer(llista, 1), "T10": primer(llista, 2), "detall": llista}
            print("B", ep, font, sortida["fonts"][font]["crescuda"], sortida["fonts"][font]["T10"], flush=True)
        res["B"].append(sortida)

    # ---------- C. falses alarmes de l'avís ----------
    if not a.sense_c:
        for fila in res["A"]:
            t = utc(fila["t"])
            if fila.get("arome", {}).get(6, 0) < 5 and fila.get("nowcast", {}).get("2", 0) < 5 and fila["obs"][6] < 5:
                continue                                   # sense pluja prevista ni mesurada que pugui moure el nivell
            try:
                r = previsio.preveu(cfg, t, "tot", desar=False, cache=cache, verb=False)
            except SystemExit as e:
                continue
            if r["nivell"] <= r["nivell_sense_prevista"]:
                res["C"].append({"t": fila["t"], "puja": False}); continue
            q6 = [v for tq, v in qobs.items() if t < tq <= t + timedelta(hours=6)]
            n_obs = 2 if q6 and max(q6) >= qt["Q10"] else 1 if q6 and max(q6) >= qt["QMCO"] else 0
            orig = pp.prevista
            pp.prevista = lambda cfg_, t0, minuts, *x, **k: (np.array([pl.get(t0 + timedelta(minutes=DT * (i + 1)), rd.get(t0 + timedelta(minutes=DT * (i + 1)), 0.0))
                                                                      for i in range(minuts // DT)]), {"mm_3h": acumulat(pl, rd, t0, 180)})
            try:
                rr = previsio.preveu(cfg, t, "tot", desar=False, verb=False)
            finally:
                pp.prevista = orig
            real = max(n_obs, rr["nivell"])
            res["C"].append({"t": fila["t"], "puja": True, "nivell": r["nivell"], "sense_prevista": r["nivell_sense_prevista"],
                             "amb_pluja_real": rr["nivell"], "aforador": n_obs, "falsa": r["nivell"] > real,
                             "T": {k: round(v, 1) for k, v in r["T"].items()}})
            print("C", fila["t"], res["C"][-1], flush=True)

    informe(cfg, res, a.sortida)
    return 0


if __name__ == "__main__":
    sys.exit(main())
