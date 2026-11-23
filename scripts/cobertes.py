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
COB_URL = "https://github.com/correvents/inundacions/releases/download/cobertes-cat/cobertes5_{}.tif"
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
    """Color → (codi, nom), mostrejant Catalunya a 25 m en trossos de 25 km. Per a cada color es pregunta GetFeatureInfo
    en 3 píxels del mig de taques homogènies (5×5 del mateix color) i es pren el codi majoritari."""
    from scipy import ndimage
    from collections import Counter
    vots = {}
    noms = {}
    x0, y0, x1, y1 = CAT
    pas, mida = 25, 1000
    for ya in range(y0, y1, pas * mida):
        for xa in range(x0, x1, pas * mida):
            box = (xa, ya, xa + pas * mida, ya + pas * mida)
            a = getmap(*box, mida, mida)
            ple = a[..., 3] > 0
            if not ple.any():
                continue
            for c in np.unique(a[ple][:, :3], axis=0):
                k = clau(c)
                if sum(vots.get(k, Counter()).values()) >= 6:
                    continue
                m = ndimage.binary_erosion(np.all(a[..., :3] == c, axis=-1) & ple, structure=np.ones((5, 5)))
                fi, ci = np.nonzero(m)
                if not len(fi):
                    continue
                for q in np.random.default_rng(0).choice(len(fi), size=min(3, len(fi)), replace=False):
                    codi, nom = featureinfo(*box, mida, mida, int(ci[q]), int(fi[q]))
                    if codi is not None:
                        vots.setdefault(k, Counter())[codi] += 1; noms[codi] = nom
    lleg = {}
    for k, v in vots.items():
        codi, n = v.most_common(1)[0]
        lleg[k] = {"codi": codi, "nom": noms[codi], "vots": f"{n}/{sum(v.values())}"}
        print(f"{k:>12} → {codi:3d} {noms[codi]}  ({n}/{sum(v.values())})", flush=True)
    LLEGENDA.write_text(json.dumps(dict(sorted(lleg.items(), key=lambda kv: kv[1]["codi"])), ensure_ascii=False, indent=1))
    print(f"{len(lleg)} colors a la llegenda")

def tradueix(a, lleg):
    """RGBA (WMS) → codi de coberta (uint8, 0 = sense dades o color desconegut)."""
    rgb = a[..., 0].astype(np.int32) << 16 | a[..., 1].astype(np.int32) << 8 | a[..., 2].astype(np.int32)
    u, inv = np.unique(rgb, return_inverse=True)
    taula = np.zeros(len(u), np.uint8)
    for k, v in enumerate(u):
        e = lleg.get("%d,%d,%d" % (v >> 16, (v >> 8) & 255, v & 255))
        taula[k] = e["codi"] if e else 0
    out = taula[inv.reshape(rgb.shape)]
    out[a[..., 3] == 0] = 0
    return out


def conques():
    fc = json.loads(CONQUES.read_text())
    return {f["properties"]["id"]: f for f in fc["features"]}


def wms_reixa(x0, y0, x1, y1, pas=2.5, tros=5000):
    """Codis de coberta en EPSG:25831 a «pas» m sobre la caixa, a trossos de «tros» m. Torna (codis, transform)."""
    from rasterio.transform import from_origin
    lleg = json.loads(LLEGENDA.read_text())
    x0, y0 = np.floor(x0 / tros) * tros, np.floor(y0 / tros) * tros
    x1, y1 = np.ceil(x1 / tros) * tros, np.ceil(y1 / tros) * tros
    n = int(tros / pas)
    nf, nc = int((y1 - y0) / pas), int((x1 - x0) / pas)
    out = np.zeros((nf, nc), np.uint8)
    for ya in np.arange(y0, y1, tros):
        for xa in np.arange(x0, x1, tros):
            a = getmap(xa, ya, xa + tros, ya + tros, n, n)
            r0, c0 = int((y1 - ya - tros) / pas), int((xa - x0) / pas)
            out[r0:r0 + n, c0:c0 + n] = tradueix(a, lleg)
    return out, from_origin(x0, y1, pas, pas)


def conca(i):
    import rasterio
    from rasterio.warp import reproject, Resampling, transform_bounds
    nom = f"c{i:03d}"
    TMP.mkdir(parents=True, exist_ok=True)
    with rasterio.open(MDT_URL.format(nom)) as f:
        prof = f.profile
        tr, crs, forma = f.transform, f.crs, (f.height, f.width)
        buit = f.read(1) == f.nodata
    caixa = transform_bounds(crs, "EPSG:25831", tr.c, tr.f + forma[0] * tr.e, tr.c + forma[1] * tr.a, tr.f)
    t = time.time()
    src, str_ = wms_reixa(caixa[0] - 20, caixa[1] - 20, caixa[2] + 20, caixa[3] + 20)
    print(f"{nom}: WMS {src.shape[1]}×{src.shape[0]} a 2,5 m en {time.time() - t:.0f} s; sense codi {np.mean(src == 0) * 100:.1f} %")
    cob = np.zeros(forma, np.uint8)
    reproject(src, cob, src_transform=str_, src_crs="EPSG:25831", src_nodata=0, dst_transform=tr, dst_crs=crs,
              dst_nodata=0, resampling=Resampling.mode)
    cob[buit] = 0
    fit = TMP / f"cobertes5_{nom}.tif"
    prof.update(dtype="uint8", nodata=0, compress="deflate", predictor=1, zlevel=9, tiled=True, blockxsize=256, blockysize=256)
    with rasterio.open(fit, "w", **prof) as f:
        f.write(cob, 1)
    print(f"{nom}: {fit} {fit.stat().st_size / 1e6:.1f} MB")
    web(nom, cob, tr, crs)


def web(nom, cob, tr, crs):
    """PNG de paleta en lat/lon (10 m) amb els grups de coberta, i resum de km² per grup i paràmetres mitjans."""
    import rasterio
    from rasterio.warp import reproject, Resampling, transform_bounds
    from PIL import Image
    par = json.loads(PARAM.read_text())
    grups = par["grups"]
    gcod = np.zeros(256, np.uint8)
    for c, e in par["classes"].items():
        gcod[int(c)] = list(grups).index(e["grup"]) + 1
    g = gcod[cob]
    w, s, e, n = transform_bounds(crs, "EPSG:4326", tr.c, tr.f + cob.shape[0] * tr.e, tr.c + cob.shape[1] * tr.a, tr.f)
    pas = 10.0
    pas_y = pas / 111320
    pas_x = pas / (111320 * np.cos(np.radians((s + n) / 2)))
    ncol, nfil = int(np.ceil((e - w) / pas_x)), int(np.ceil((n - s) / pas_y))
    gl = np.zeros((nfil, ncol), np.uint8)
    reproject(g, gl, src_transform=tr, src_crs=crs, src_nodata=0, dst_transform=rasterio.Affine(pas_x, 0, w, 0, -pas_y, n),
              dst_crs="EPSG:4326", dst_nodata=0, resampling=Resampling.mode)
    d = SORT / nom
    d.mkdir(parents=True, exist_ok=True)
    im = Image.fromarray(gl, "P")
    pal = [0, 0, 0] * 256
    for k, gr in enumerate(grups.values(), 1):
        pal[3 * k:3 * k + 3] = gr["color"]
    im.putpalette(pal)
    im.save(d / "cobertes.png", optimize=True, transparency=bytes([0] + [210] * len(grups)))
    cel = abs(tr.a * tr.e) / 1e6
    km2 = {k: round(float((g == j).sum() * cel), 2) for j, k in enumerate(grups, 1)}
    ok = cob > 0
    nman, cn = parametres(cob, par)[:2]
    resum = {"conca": nom, "limits": [[n - nfil * pas_y, w], [n, w + ncol * pas_x]], "km2": km2,
             "n_manning_mitja": round(float(nman[ok].mean()), 3) if ok.any() else None,
             "cn_mitja": round(float(cn[ok].mean()), 1) if ok.any() else None,
             "font": f"Mapa de cobertes del sòl de Catalunya ({CAPA[-4:]}), ICGC, CC BY 4.0", "calculat": time.strftime("%Y-%m-%d")}
    (d / "resum.json").write_text(json.dumps(resum, ensure_ascii=False, separators=(",", ":")))
    print(f"{nom}: web {(d / 'cobertes.png').stat().st_size / 1e3:.0f} kB; " + ", ".join(f"{k} {v}" for k, v in km2.items() if v))


def parametres(cob, par=None, grup_sol=None):
    """Per a cada cel·la: n de Manning, número de corba (TR-55, grup hidrològic de sòl) i tipus (0 natural, 1 urbà, 2 aigua).
    Les cel·les sense codi prenen el valor de «sense_dades» del fitxer de paràmetres.
    grup_sol pot ser una lletra (A-D, tota la conca igual) o una matriu com cob amb el grup per cel·la (0=A … 3=D;
    -1 = sense dades → el grup per defecte del fitxer de paràmetres), p. ex. de sols.a_reixa()."""
    par = par or json.loads(PARAM.read_text())
    gdef = "ABCD".index(par["grup_hidrologic_sol"])
    if isinstance(grup_sol, np.ndarray):
        gs = np.where(grup_sol < 0, gdef, grup_sol).astype(np.intp)
    else:
        gs = "ABCD".index(grup_sol) if grup_sol else gdef
    tn, tc, tt = np.zeros(256), np.zeros((256, 4)), np.zeros(256, np.uint8)
    sd = par["classes"][str(par["sense_dades"])]
    tn[:], tc[:], tt[:] = sd["n"], sd["cn"], {"natural": 0, "urba": 1, "aigua": 2}[sd["tipus"]]
    for c, e in par["classes"].items():
        tn[int(c)], tc[int(c)], tt[int(c)] = e["n"], e["cn"], {"natural": 0, "urba": 1, "aigua": 2}[e["tipus"]]
    return tn[cob], tc[cob, gs], tt[cob]


def a_reixa(tr, forma, crs):
    """Codis de coberta a una reixa qualsevol (per a simula_2d.py), a partir dels fitxers de la Release cobertes-cat
    de les conques que hi toquen (moda). Torna None si no n'hi ha cap."""
    import rasterio
    from rasterio.warp import reproject, Resampling, transform_bounds
    from shapely.geometry import shape, box
    caixa = transform_bounds(crs, "EPSG:4326", tr.c, tr.f + forma[0] * tr.e, tr.c + forma[1] * tr.a, tr.f)
    b = box(*caixa)
    out = np.zeros(forma, np.uint8)
    n = 0
    for i, f in conques().items():
        if not shape(f["geometry"]).intersects(b):
            continue
        url = COB_URL.format(f"c{i:03d}")
        try:
            with rasterio.open(url) as src:
                tmp = np.zeros(forma, np.uint8)
                reproject(rasterio.band(src, 1), tmp, dst_transform=tr, dst_crs=crs, src_nodata=0, dst_nodata=0,
                          resampling=Resampling.mode)
        except Exception as e:
            print(f"Cobertes c{i:03d}: {type(e).__name__} {e}"); continue
        out = np.where(out == 0, tmp, out); n += 1
    return out if n else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="ordre", required=True)
    sub.add_parser("llegenda")
    a1 = sub.add_parser("lots"); a1.add_argument("--lots", type=int, default=30); a1.add_argument("--conques", default="pendents")
    a2 = sub.add_parser("conca"); a2.add_argument("--id", type=int, required=True)
    sub.add_parser("index")
    a = ap.parse_args()
    if a.ordre == "llegenda":
        llegenda()
    elif a.ordre == "lots":
        C = conques()
        if a.conques == "pendents":
            a.conques = ",".join(str(i) for i in sorted(C) if not (SORT / f"c{i:03d}" / "resum.json").exists()) or "totes"
        ids = sorted(C) if a.conques == "totes" else [int(x) for x in a.conques.split(",") if x.strip()]
        lots = [[] for _ in range(min(a.lots, len(ids)))]
        pes = [0.0] * len(lots)
        for i in sorted(ids, key=lambda i: -C[i]["properties"]["km2"]):
            k = int(np.argmin(pes)); lots[k].append(i); pes[k] += C[i]["properties"]["km2"] + 20
        print(json.dumps({"lot": [" ".join(map(str, l)) for l in lots]}))
    elif a.ordre == "conca":
        conca(a.id)
    elif a.ordre == "index":
        par = json.loads(PARAM.read_text())
        idx = {}
        for r in sorted(SORT.glob("c*/resum.json")):
            j = json.loads(r.read_text())
            idx[j["conca"]] = {k: j[k] for k in ("limits", "km2", "n_manning_mitja", "cn_mitja")}
        (SORT / "index.json").write_text(json.dumps({"conques": idx, "grups": par["grups"], "font": par["font_cobertes"]},
                                                    ensure_ascii=False, separators=(",", ":")))
        print(f"{len(idx)} conques amb cobertes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
