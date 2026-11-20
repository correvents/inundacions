#!/usr/bin/env python3
"""Mapa de cobertes del sòl de Catalunya (ICGC, CC BY 4.0) a la reixa de l'MDT de 5 m, conca per conca.

La descàrrega directa (datacloud.icgc.cat) no respon des de GitHub Actions, però el WMS de l'ICGC sí, i és ràster:
cada coberta té un color únic i GetFeatureInfo en torna el codi i el nom. Per tant:

1. llegenda: recorre Catalunya a baixa resolució, i per a cada color nou demana GetFeatureInfo en aquell píxel →
   config/cobertes_llegenda.json (color RGB → codi, nom). Es fa un cop; després es reutilitza.
2. conca --id N: demana el WMS (capa cobertes_AAAA, EPSG:25831) a 2,5 m a trossos de 5 km sobre la caixa de la conca,
   tradueix colors a codis i passa a la reixa exacta de l'MDT de la conca (Release mdt5-cat, EPSG:25830, 5 m) per
   moda → /tmp/cobertes/cobertes5_cNNN.tif (uint8, 0 = sense dades), que el workflow puja a la Release cobertes-cat.
   També escriu la capa web data/catalunya/cobertes/cNNN/{cobertes.png, resum.json}.
3. index: data/catalunya/cobertes/index.json per a conques.html.

Els paràmetres hidrològics de cada coberta (n de Manning i número de corba) són a config/cobertes_parametres.json,
de taules publicades (Chow 1959 / HEC-RAS; USDA TR-55 1986).

Ús:  python3 scripts/cobertes.py llegenda
     python3 scripts/cobertes.py lots --lots 30 [--conques totes|pendents|1,2,3]
     python3 scripts/cobertes.py conca --id 124
     python3 scripts/cobertes.py index"""
import argparse, io, json, re, sys, time
from pathlib import Path
import numpy as np, requests

WMS = "https://geoserveis.icgc.cat/servei/catalunya/cobertes-sol/wms"
CAPA = "cobertes_2024"
LLEGENDA = Path("config/cobertes_llegenda.json")
PARAM = Path("config/cobertes_parametres.json")
CONQUES = Path("data/catalunya/conques/conques_100km2.geojson")
SORT = Path("data/catalunya/cobertes")
TMP = Path("/tmp/cobertes")
MDT_URL = "https://github.com/correvents/inundacions/releases/download/mdt5-cat/mdt5_{}.tif"
CAT = (258000, 4485000, 530000, 4752000)        # caixa de Catalunya en EPSG:25831


def getmap(x0, y0, x1, y1, w, h, capa=CAPA):
    p = {"SERVICE": "WMS", "VERSION": "1.3.0", "REQUEST": "GetMap", "LAYERS": capa, "STYLES": "", "CRS": "EPSG:25831",
         "BBOX": f"{x0},{y0},{x1},{y1}", "WIDTH": w, "HEIGHT": h, "FORMAT": "image/png"}
    from PIL import Image
    for intent in range(6):
        try:
            r = requests.get(WMS, params=p, timeout=180)
            if r.content[:4] == b"\x89PNG":
                return np.array(Image.open(io.BytesIO(r.content)).convert("RGBA"))
        except requests.RequestException:
            pass
        time.sleep(5 * (intent + 1))
    raise RuntimeError(f"WMS {x0},{y0}: {getattr(r, 'status_code', '?')} {getattr(r, 'text', '')[:200]}")


def featureinfo(x0, y0, x1, y1, w, h, i, j, capa=CAPA):
    p = {"SERVICE": "WMS", "VERSION": "1.3.0", "REQUEST": "GetFeatureInfo", "LAYERS": capa, "QUERY_LAYERS": capa,
         "STYLES": "", "CRS": "EPSG:25831", "BBOX": f"{x0},{y0},{x1},{y1}", "WIDTH": w, "HEIGHT": h, "I": i, "J": j,
         "INFO_FORMAT": "text/plain"}
    r = requests.get(WMS, params=p, timeout=60)
    v = re.search(r"value_0 = '([^']*)'", r.text)
    c = re.search(r"class = '(.*)'\s*$", r.text, re.M)
    return (int(float(v.group(1))) if v else None), (c.group(1) if c else None)


def clau(rgb):
    return "%d,%d,%d" % tuple(int(x) for x in rgb[:3])


def llegenda():
    """Color → (codi, nom), mostrejant Catalunya a 50 m en trossos de 50 km."""
    lleg = json.loads(LLEGENDA.read_text()) if LLEGENDA.exists() else {}
    x0, y0, x1, y1 = CAT
    pas, mida = 50, 1000
    for ya in range(y0, y1, pas * mida):
        for xa in range(x0, x1, pas * mida):
            box = (xa, ya, xa + pas * mida, ya + pas * mida)
            a = getmap(*box, mida, mida)
            ple = a[..., 3] > 0
            if not ple.any():
                continue
            cols = np.unique(a[ple][:, :3], axis=0)
            for c in cols:
                k = clau(c)
                if k in lleg:
                    continue
                fi, ci = np.nonzero(np.all(a[..., :3] == c, axis=-1) & ple)
                # un píxel envoltat del mateix color (no a la vora d'un polígon)
                m = len(fi) // 2
                codi, nom = featureinfo(*box, mida, mida, int(ci[m]), int(fi[m]))
                if codi is not None:
                    lleg[k] = {"codi": codi, "nom": nom}
                    print(f"{k:>12} → {codi:3d} {nom}", flush=True)
    LLEGENDA.write_text(json.dumps(dict(sorted(lleg.items(), key=lambda kv: kv[1]["codi"])), ensure_ascii=False, indent=1))
    print(f"{len(lleg)} colors a la llegenda")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="ordre", required=True)
    sub.add_parser("llegenda")
    a = ap.parse_args()
    if a.ordre == "llegenda":
        llegenda()
    return 0


if __name__ == "__main__":
    sys.exit(main())
