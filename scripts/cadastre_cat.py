#!/usr/bin/env python3
"""Edificis del Cadastre de tot Catalunya, municipi per municipi (INSPIRE Buildings, mateix lector que cadastre.py).

Per a cada municipi de les gerències 08, 17, 25 i 43 baixa el zip INSPIRE i en desa les parts d'edifici
(BuildingPart) en UTM 31N, amb plantes sobre i sota rasant:
  /tmp/cadastre_cat/edificis_MMMMM.json.gz   {"crs","camps":["ref","plantes_sobre","plantes_sota","anells"],"parts":[...]}
El workflow els puja a la Release `cadastre-cat`. Cada lot escriu també un índex amb el requadre de cada municipi,
i `index` els ajunta a data/catalunya/cadastre/municipis.json (per triar quins municipis toquen cada conca).

Ús:  python3 scripts/cadastre_cat.py llista --lots 16                 → matriu JSON per al workflow
     python3 scripts/cadastre_cat.py baixa --lot 3 --lots 16 [--salta fitxer_amb_noms_ja_pujats]
     python3 scripts/cadastre_cat.py index /tmp/lots/*.json"""
import argparse, gzip, json, re, sys
from pathlib import Path
import numpy as np
from lxml import etree

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cadastre as C

PROVS = [8, 17, 25, 43]
TMP = Path("/tmp/cadastre_cat")
SORT = Path("data/catalunya/cadastre")
TOT = [-1e12, -1e12, 1e12, 1e12]                  # sense filtre de caixa: tot el municipi


def tots():
    """[(prov, codi, nom, url)] de tots els municipis de Catalunya, en ordre estable."""
    out = []
    for p in PROVS:
        xml = C.xml_de(C.feed_provincial(p))
        for e in xml.iter():
            if C.ln(e) != "entry":
                continue
            titol = next((x.text for x in e if C.ln(x) == "title"), "") or ""
            enl = [x.get("href") for x in e if C.ln(x) == "link" and (x.get("href") or "").lower().endswith((".zip", ".xml"))]
            if not enl:
                continue
            m = re.search(r"(\d{5})", enl[0].rsplit("/", 1)[-1]) or re.search(r"(\d{5})", titol)
            if not m:
                continue
            out.append((p, m.group(1), titol.strip(), enl[0]))
    out.sort(key=lambda t: t[1])
    print(f"{len(out)} municipis", file=sys.stderr)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="ordre", required=True)
    a1 = sub.add_parser("llista"); a1.add_argument("--lots", type=int, default=16)
    a2 = sub.add_parser("baixa"); a2.add_argument("--lot", type=int, required=True); a2.add_argument("--lots", type=int, required=True)
    a2.add_argument("--salta", help="fitxer amb els noms d'assets que ja hi ha a la Release")
    a3 = sub.add_parser("index"); a3.add_argument("fitxers", nargs="+")
    a = ap.parse_args()

    if a.ordre == "llista":
        print(json.dumps({"lot": list(range(a.lots))}))
        return 0

    if a.ordre == "index":
        idx = {}
        for f in a.fitxers:
            idx.update(json.loads(Path(f).read_text()))
        SORT.mkdir(parents=True, exist_ok=True)
        (SORT / "municipis.json").write_text(json.dumps(dict(sorted(idx.items())), ensure_ascii=False, separators=(",", ":")))
        print(f"índex: {len(idx)} municipis, {sum(v['parts'] for v in idx.values())} parts d'edifici")
        return 0

    TMP.mkdir(parents=True, exist_ok=True)
    ja = set(Path(a.salta).read_text().split()) if a.salta and Path(a.salta).exists() else set()
    prev = json.loads((SORT / "municipis.json").read_text()) if (SORT / "municipis.json").exists() else {}
    llista = tots()[a.lot::a.lots]
    idx, fallats = {}, []
    for prov, codi, nom, url in llista:
        f = TMP / f"edificis_{codi}.json.gz"
        if f.name in ja and codi in prev:
            idx[codi] = prev[codi]; continue
        try:
            parts, _ = C.llegeix(C.get(C.zip_de(url)).content, TOT)
        except Exception as ex:
            print(f"{codi} {nom}: {type(ex).__name__} {ex}"); fallats.append(codi); continue
        if not parts:
            print(f"{codi} {nom}: cap edifici"); continue
        xy = np.concatenate([np.asarray(p[3][0]) for p in parts])
        b = [round(float(v), 1) for v in (*xy.min(axis=0), *xy.max(axis=0))]
        with gzip.open(f, "wt", compresslevel=9) as g:
            json.dump({"crs": "EPSG:25831", "municipi": nom, "camps": ["ref", "plantes_sobre", "plantes_sota", "anells"], "parts": parts},
                      g, separators=(",", ":"))
        idx[codi] = {"nom": nom, "prov": prov, "bbox": b, "parts": len(parts)}
        print(f"{codi} {nom}: {len(parts)} parts, {f.stat().st_size / 1e6:.1f} MB")
    (TMP / f"index_{a.lot}.json").write_text(json.dumps(idx, ensure_ascii=False))
    (TMP / f"fallats_{a.lot}.txt").write_text(" ".join(fallats))
    print(f"Lot {a.lot}: {len(idx)} municipis, fallats: {' '.join(fallats) or 'cap'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
