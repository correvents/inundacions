#!/usr/bin/env python3
"""Resum curt (≤ 20 línies, Markdown) d'una simulació 2D, per llegir-lo en lloc dels fitxers grans o del log.
Llegeix data/<conca>/sim2d/<id>/index.json (i soterranis.json si hi és), escriu resum.md a la mateixa carpeta i l'imprimeix.
Ús:  python3 scripts/resum.py data/ridaura/sim2d/2026-10-04 [més carpetes...]"""
import json, sys
from pathlib import Path


def resum_sim2d(carpeta):
    d = Path(carpeta); ix = json.loads((d / "index.json").read_text())
    fot = ix.get("fotogrames") or []; q = ix.get("q_aforador_sim") or []
    l = [f"# Sim 2D {d.name} ({d.parent.parent.name})", "",
         f"- Pas {ix.get('pas_m')} m, inici {ix.get('inici_utc')} UTC, {len(fot)} fotogrames"]
    if q and max(q) > 0:
        k = max(range(len(q)), key=q.__getitem__)
        l.append(f"- Pic a l'aforador: {q[k]} m³/s a les {fot[k]['temps'][:16].replace('T', ' ') if k < len(fot) else '?'} UTC")
    elif q:
        l.append("- Sense cabal a l'aforador")
    g = ix.get("graella_utm"); f = d / "calat_max_cm.u2"
    if g and f.exists():
        import numpy as np
        h = np.fromfile(f, "<u2"); a = (ix.get("pas_m") or 0) ** 2 / 1e6
        l.append("- Superfície amb calat màxim > 15 / 30 / 60 cm: " + " / ".join(f"{(h > c).sum() * a:.2f}" for c in (15, 30, 60))
                 + f" km²; màxim {h[h < 65534].max() / 100:.1f} m")
    if ix.get("assecat_final_mullat_km2") is not None:
        l.append(f"- Encara mullat al final (> 5 cm): {ix['assecat_final_mullat_km2']} km²")
    b = ix.get("balanc_hm3")
    if b:
        l.append("- Balanç (hm³): " + ", ".join(f"{k} {v}" for k, v in b.items()))
    e = ix.get("edificis")
    if e:
        l.append(f"- Edificis: ocupen {e.get('ocupacio_pct')} %, soterranis amb ≥ 15 cm a tocar: {e.get('soterranis_risc')}")
    if ix.get("plantes_baixes"):
        l.append(f"- Plantes baixes: {json.dumps(ix['plantes_baixes'], ensure_ascii=False)[:300]}")
    v = ix.get("velocitat")
    if v:
        l.append(f"- Velocitat p99 {v.get('velocitat_p99_ms')} m/s; perill (km²) {v.get('km2_per_classe')}")
    h = ix.get("hipotesis") or {}
    claus = ["ia_natural_mm", "ia_urba_mm", "perdua_natural_mm_h", "clavegueram_escenari", "clavegueram_desti", "pluja"]
    l.append("- Hipòtesis: " + "; ".join(f"{c} {h[c]}" for c in claus if c in h))
    p = h.get("passos")
    if p:
        l.append(f"- Passos: {p if isinstance(p, str) else p.get('mode')}")
    text = "\n".join(l) + "\n"
    (d / "resum.md").write_text(text)
    return text


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for c in sys.argv[1:]:
        print(resum_sim2d(c))
