#!/usr/bin/env python3
"""Índex d'episodis d'una conca nova (data/<conca>/episodis/index.json) per a evolucio.html.

Per a cada episodi calcula l'absorció del terreny com premsa.py (pluja antecedent dels 60 dies anteriors, API) amb l'estació XEMA
més propera amb dades (cfg["estacio_api"] o, si és buida, la que tingui més dades als 60 dies abans del primer episodi).

Ús:  python3 scripts/index_episodis.py --config config/granollers.json \\
        --episodi '{"id":"2026-10-04","titol":"...","resum":"...","radar":"2026-10-04_radar.json"}' [--episodi ...]"""
import argparse, json
from pathlib import Path
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--episodi", action="append", required=True)
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text()); d = Path("data") / cfg["conca"]
    eps = [json.loads(e) for e in a.episodi]
    ant = d / "episodis" / "index.json"                      # conserva la sèrie mesurada (episodi) si ja hi era
    previ = {e["id"]: e for e in json.loads(ant.read_text())["episodis"]} if ant.exists() else {}
    for e in eps:
        if "episodi" not in e and previ.get(e["id"], {}).get("episodi"):
            e["episodi"] = previ[e["id"]]["episodi"]
    est = {e["codi_estacio"]: e for e in json.loads((d / "xema" / "estacions.json").read_text())}
    codi = cfg.get("estacio_api")
    t0s = {e["id"]: pd.Timestamp(json.loads((d / "episodis" / e["radar"]).read_text())["temps"][0]) for e in eps}
    def serie(c):
        f = d / "xema" / f"{c}.csv"
        return pd.read_csv(f, parse_dates=["temps"]).set_index("temps").mm if f.exists() else None
    def cobertura(c):
        s = serie(c)
        if s is None: return -1
        return sum(len(s[t - pd.Timedelta(days=60): t]) for t in t0s.values())
    if not codi:
        codi = max(est, key=cobertura)
    s = serie(codi)
    print(f"estació per a la pluja antecedent: {codi} {est[codi]['nom_estacio']} ({est[codi]['km']} km)")
    for e in eps:
        t0 = t0s[e["id"]]
        dies = s[t0 - pd.Timedelta(days=60): t0 - pd.Timedelta(minutes=1)].resample("1D", closed="right", label="right").sum()
        api = float(sum(v * 0.92 ** ((t0 - t).days) for t, v in dies.items()))
        e.update({"api_mm": round(api), "ia_mm": round(150 - 100 * min(api / 100, 1)), "perdua_mm_h": 12, "estacio_api": codi})
        print(f"{e['id']}: API {api:.0f} mm → absorbeix {e['ia_mm']} mm")
    nous = {e["id"] for e in eps}                              # els episodis d'abans que no es refan es conserven
    eps = sorted(eps + [e for i, e in previ.items() if i not in nous], key=lambda e: e["id"])
    (d / "episodis" / "index.json").write_text(json.dumps({"episodis": eps, "sense_radar": [], "pendents_de_radar": []}, ensure_ascii=False, indent=1))
    if not (d / "episodis" / "premsa.geojson").exists():
        (d / "episodis" / "premsa.geojson").write_text('{"type":"FeatureCollection","features":[]}')
    return 0


if __name__ == "__main__":
    main()
