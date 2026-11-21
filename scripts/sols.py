#!/usr/bin/env python3
"""Mapa de grups hidrològics de sòl (A-D, USDA-NRCS) de tot Catalunya, per al número de corba de simula_2d.py.

Font: SoilGrids 2.0 (ISRIC, CC BY 4.0, 250 m): fracció d'argila i de sorra a 0-5, 5-15, 15-30 i 30-60 cm,
i, si s'hi pot accedir, la fondària a la roca de SoilGrids 2017 (BDTICM, Shangguan et al. 2017).
Res ajustat a mà: el grup surt de les regles del National Engineering Handbook, part 630, capítol 7 (USDA-NRCS 2009),
les mateixes que fa servir HYSOGs250m (Ross et al. 2018):
   A: argila < 10 % i sorra > 90 %
   B: argila 10-20 % i sorra 50-90 %   (i la resta amb argila ≤ 20 % i sorra ≥ 50 %)
   C: argila 20-40 % i sorra < 50 %    (i la resta de textures mitjanes)
   D: argila > 40 %  ·  o roca/capa impermeable a menys de 50 cm
Es desen dues variants (bandes del GeoTIFF), perquè s'avaluïn amb episodis:
   1 «mapa»:     NEH-630: la capa menys permeable dels primers 60 cm, i D si la roca és a < 50 cm.
   2 «mapa-sup»: textura mitjana ponderada dels primers 30 cm (sense regla de roca).
Valors: 1=A, 2=B, 3=C, 4=D, 0 = sense dades (aigua; simula_2d.py hi posa el grup per defecte).
Les cel·les sense dades de SoilGrids a terra (nuclis urbans grans) prenen el grup de la cel·la vàlida més propera.

Sortida: /tmp/sols/grups_sol_cat.tif (EPSG:25831, 250 m, 2 bandes uint8) i /tmp/sols/textura_sol_cat.tif
(argila i sorra 0-30 cm en %, fondària a la roca en cm) per a la Release `sols-cat`; capa web data/catalunya/sols/
(grups.png en EPSG:3857, resum.json amb % de cada grup per conca de 100 km²).

Variant 3 «litologia» (ordre `litologia`): SoilGrids és molt suau a Catalunya (argila 20-32 %, sorra 20-40 % a tot
arreu → gairebé tot C), així que també es fa el grup a partir del Mapa de Permeabilitats d'Espanya 1:200.000 de l'IGME
(fet del Mapa litoestratigràfic, 5 classes), com fa la norma 5.2-IC (2016) amb la litologia. Correspondència ordinal,
a config/sols_litologia.json: molt alta → A, alta → B, mitjana → C, baixa i molt baixa → D.
→ /tmp/sols/grups_sol_litologia_cat.tif (EPSG:25831, 50 m, 1=A … 4=D) i data/catalunya/sols/grups_lit.png.

Ús:  python3 scripts/sols.py construeix
     python3 scripts/sols.py litologia
     python3 scripts/sols.py prova                  # comprova l'accés a les fonts"""
import argparse, json, os, sys, time
from pathlib import Path
import numpy as np

os.environ.setdefault("GDAL_HTTP_MAX_RETRY", "5")
os.environ.setdefault("GDAL_HTTP_RETRY_DELAY", "3")
os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".vrt,.tif")
os.environ.setdefault("GDAL_HTTP_TIMEOUT", "120")

SG = "/vsicurl/https://files.isric.org/soilgrids/latest/data/{p}/{p}_{c}_mean.vrt"
WCS = ("https://maps.isric.org/mapserv?map=/map/{p}.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage"
       "&COVERAGEID={p}_{c}_mean&FORMAT=GEOTIFF_INT16&SUBSET=X({x0},{x1})&SUBSET=Y({y0},{y1})"
       "&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/25831&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/25831"
       "&RESOLUTION=X(250)&RESOLUTION=Y(250)")
ROCA = ["/vsicurl/https://files.isric.org/soilgrids/former/2017-03-10/data/BDTICM_M_250m_ll.tif",
        "/vsicurl/https://files.isric.org/soilgrids/former/2017-03-10/data/BDTICM_M_250m.tif"]
CAPES = [("0-5cm", 5), ("5-15cm", 10), ("15-30cm", 15), ("30-60cm", 30)]
CAT = (258000, 4485000, 530000, 4752000)        # caixa de Catalunya en EPSG:25831 (la mateixa que cobertes.py)
PAS = 250
CRS = "EPSG:25831"
TMP = Path("/tmp/sols")
SORT = Path("data/catalunya/sols")
CONQUES = Path("data/catalunya/conques/conques_100km2.geojson")
REL = "https://github.com/correvents/inundacions/releases/download/sols-cat/"
VARIANTS = {"mapa": ("grups_sol_cat.tif", 1), "mapa-sup": ("grups_sol_cat.tif", 2), "litologia": ("grups_sol_litologia_cat.tif", 1)}
IGME = "https://mapas.igme.es/gis/rest/services/Cartografia_Tematica/IGME_Permeabilidad_200/MapServer/0"
LIT = Path("config/sols_litologia.json")
PAS_LIT = 50
COLORS = {1: (255, 237, 160), 2: (161, 217, 155), 3: (65, 171, 93), 4: (0, 90, 50)}


def reixa():
    from rasterio.transform import from_origin
    x0, y0, x1, y1 = CAT
    return from_origin(x0, y1, PAS, PAS), ((y1 - y0) // PAS, (x1 - x0) // PAS)


def llegeix(p, c):
    """Propietat SoilGrids (g/kg) a la reixa de Catalunya, en %. Primer el VRT (COG) amb /vsicurl, si no el WCS."""
    import rasterio
    from rasterio.warp import reproject, Resampling
    tr, forma = reixa()
    out = np.full(forma, np.nan, np.float32)
    err = []
    for intent in range(3):
        try:
            with rasterio.open(SG.format(p=p, c=c)) as src:
                reproject(rasterio.band(src, 1), out, dst_transform=tr, dst_crs=CRS, src_nodata=src.nodata,
                          dst_nodata=np.nan, resampling=Resampling.bilinear)
            return out / 10
        except Exception as e:
            err.append(f"VRT {type(e).__name__} {e}"); time.sleep(10 * (intent + 1))
    import requests, rasterio.io
    x0, y0, x1, y1 = CAT
    for intent in range(4):
        try:
            r = requests.get(WCS.format(p=p, c=c, x0=x0, x1=x1, y0=y0, y1=y1), timeout=300)
            with rasterio.io.MemoryFile(r.content) as m, m.open() as src:
                reproject(rasterio.band(src, 1), out, dst_transform=tr, dst_crs=CRS, src_nodata=src.nodata,
                          dst_nodata=np.nan, resampling=Resampling.bilinear)
            return out / 10
        except Exception as e:
            err.append(f"WCS {type(e).__name__} {str(e)[:200]}"); time.sleep(10 * (intent + 1))
    raise RuntimeError(f"SoilGrids {p} {c}: " + " | ".join(err))


def roca():
    """Fondària absoluta a la roca (cm) de SoilGrids 2017 a la reixa, o None si no s'hi pot accedir."""
    import rasterio
    from rasterio.warp import reproject, Resampling
    tr, forma = reixa()
    for u in ROCA:
        try:
            with rasterio.open(u) as src:
                out = np.full(forma, np.nan, np.float32)
                reproject(rasterio.band(src, 1), out, dst_transform=tr, dst_crs=CRS, src_nodata=src.nodata,
                          dst_nodata=np.nan, resampling=Resampling.bilinear)
            if np.isfinite(out).mean() > 0.2:
                print(f"Fondària a la roca: {u}")
                return out
        except Exception as e:
            print(f"Roca {u}: {type(e).__name__} {str(e)[:200]}")
    return None


def grup(argila, sorra):
    """NEH-630 cap. 7: textura → grup (1-4); 0 si no hi ha dades."""
    g = np.full(argila.shape, 3, np.uint8)
    g[(argila <= 20) & (sorra >= 50)] = 2
    g[(argila < 10) & (sorra > 90)] = 1
    g[argila > 40] = 4
    g[~(np.isfinite(argila) & np.isfinite(sorra))] = 0
    return g


def omple(g, terra):
    """Cel·les de terra sense dades (nuclis urbans de SoilGrids) → grup de la cel·la vàlida més propera."""
    from scipy import ndimage
    buit = (g == 0) & terra
    if buit.any() and (g > 0).any():
        _, (ii, jj) = ndimage.distance_transform_edt(g == 0, return_indices=True)
        g = np.where(buit, g[ii, jj], g)
    return g


def construeix():
    import rasterio
    TMP.mkdir(parents=True, exist_ok=True)
    tr, forma = reixa()
    cl, sa = {}, {}
    for c, _ in CAPES:
        cl[c], sa[c] = llegeix("clay", c), llegeix("sand", c)
        print(f"SoilGrids {c}: argila {np.nanmean(cl[c]):.1f} %, sorra {np.nanmean(sa[c]):.1f} %, "
              f"vàlid {np.isfinite(cl[c]).mean() * 100:.0f} %")
    # variant 1 (NEH-630): la capa menys permeable dels primers 60 cm (el grup més alt), i D amb roca a < 50 cm
    gs = np.stack([grup(cl[c], sa[c]) for c, _ in CAPES])
    g1 = np.where((gs > 0).all(0), gs.max(0), 0).astype(np.uint8)
    bdr = roca()
    if bdr is not None:
        g1[(g1 > 0) & (bdr < 50)] = 4
    # variant 2: textura mitjana dels primers 30 cm
    w = np.array([e for _, e in CAPES[:3]], float)
    a30 = sum(cl[c] * e for (c, _), e in zip(CAPES[:3], w)) / w.sum()
    s30 = sum(sa[c] * e for (c, _), e in zip(CAPES[:3], w)) / w.sum()
    g2 = grup(a30, s30)
    terra = dins_catalunya(tr, forma)
    g1, g2 = omple(g1, terra) * terra, omple(g2, terra) * terra
    perfil = dict(driver="GTiff", width=forma[1], height=forma[0], crs=CRS, transform=tr, compress="deflate")
    with rasterio.open(TMP / "grups_sol_cat.tif", "w", count=2, dtype="uint8", nodata=0, **perfil) as f:
        f.write(g1, 1); f.write(g2, 2)
        f.set_band_description(1, "mapa: NEH-630, capa menys permeable 0-60 cm + roca < 50 cm")
        f.set_band_description(2, "mapa-sup: textura mitjana 0-30 cm")
        f.update_tags(valors="1=A 2=B 3=C 4=D 0=sense dades", font="SoilGrids 2.0 (ISRIC, CC BY 4.0); BDTICM SoilGrids 2017",
                      metode="USDA-NRCS NEH-630 cap. 7 (2009)")
    with rasterio.open(TMP / "textura_sol_cat.tif", "w", count=3, dtype="float32", nodata=np.nan, **perfil) as f:
        f.write(a30.astype(np.float32), 1); f.write(s30.astype(np.float32), 2)
        f.write((bdr if bdr is not None else np.full(forma, np.nan)).astype(np.float32), 3)
        for i, d in enumerate(["argila 0-30 cm (%)", "sorra 0-30 cm (%)", "fondària a la roca (cm)"], 1):
            f.set_band_description(i, d)
    for nom, g in (("mapa", g1), ("mapa-sup", g2)):
        v = g[g > 0]
        print(f"{nom}: " + ", ".join(f"{'ABCD'[k - 1]} {(v == k).mean() * 100:.0f} %" for k in range(1, 5)))
    web(g1, g2, tr, forma, roca_ok=bdr is not None)


def dins_catalunya(tr, forma):
    """Màscara de les cel·les dins d'alguna conca del tall de 100 km² (Catalunya terra)."""
    from rasterio.features import rasterize
    from shapely.geometry import shape
    from shapely.ops import transform as stransform
    from pyproj import Transformer
    t = Transformer.from_crs(4326, CRS, always_xy=True).transform
    fc = json.loads(CONQUES.read_text())
    geoms = [stransform(t, shape(f["geometry"])) for f in fc["features"]]
    return rasterize([(g, 1) for g in geoms], out_shape=forma, transform=tr, fill=0, dtype="uint8").astype(bool)


def web(g1, g2, tr, forma, roca_ok):
    """grups.png i grups_sup.png en EPSG:3857 (com el relleu, per Leaflet) i resum.json amb % per conca."""
    from rasterio.warp import reproject, Resampling, calculate_default_transform, transform_bounds
    from rasterio.features import rasterize
    from shapely.geometry import shape
    from shapely.ops import transform as stransform
    from pyproj import Transformer
    from PIL import Image
    SORT.mkdir(parents=True, exist_ok=True)
    x0, y0, x1, y1 = CAT
    t3, w3, h3 = calculate_default_transform(CRS, "EPSG:3857", forma[1], forma[0], x0, y0, x1, y1)
    for nom, g in (("grups", g1), ("grups_sup", g2)):
        o = np.zeros((h3, w3), np.uint8)
        reproject(g, o, src_transform=tr, src_crs=CRS, dst_transform=t3, dst_crs="EPSG:3857", src_nodata=0,
                  dst_nodata=0, resampling=Resampling.nearest)
        im = Image.fromarray(o, "P")
        pal = [0, 0, 0] + sum((list(COLORS[k]) for k in range(1, 5)), [])
        im.putpalette(pal + [0] * (768 - len(pal)))
        im.save(SORT / f"{nom}.png", optimize=True, transparency=bytes([0, 200, 200, 200, 200]))
    lo0, la0, lo1, la1 = transform_bounds("EPSG:3857", 4326, t3.c, t3.f + h3 * t3.e, t3.c + w3 * t3.a, t3.f)
    t = Transformer.from_crs(4326, CRS, always_xy=True).transform
    fc = json.loads(CONQUES.read_text())
    conq = {}
    for f in fc["features"]:
        m = rasterize([(stransform(t, shape(f["geometry"])), 1)], out_shape=forma, transform=tr, fill=0,
                      dtype="uint8").astype(bool)
        r = {}
        for nom, g in (("mapa", g1), ("mapa-sup", g2)):
            v = g[m & (g > 0)]
            r[nom] = [round(float((v == k).mean() * 100), 1) for k in range(1, 5)] if v.size else None
        conq["c%03d" % f["properties"]["id"]] = r
    resum = {"limits": [[la0, lo0], [la1, lo1]], "colors": {"ABCD"[k - 1]: list(COLORS[k]) for k in range(1, 5)},
             "font": "SoilGrids 2.0 (ISRIC, CC BY 4.0), 250 m" + ("; fondària a la roca SoilGrids 2017" if roca_ok else ""),
             "metode": "USDA-NRCS NEH-630 cap. 7: textura → grup A-D; «mapa» = capa menys permeable 0-60 cm"
                       + (" i D amb roca a < 50 cm" if roca_ok else "") + "; «mapa-sup» = textura mitjana 0-30 cm",
             "roca": roca_ok, "conques": conq, "calculat": time.strftime("%Y-%m-%d")}
    (SORT / "resum.json").write_text(json.dumps(resum, ensure_ascii=False, separators=(",", ":")))
    print(f"Web: {sum(p.stat().st_size for p in SORT.iterdir()) / 1e3:.0f} kB")


def a_reixa(tr, forma, crs, variant="mapa"):
    """Grup hidrològic per cel·la (0=A … 3=D; -1 sense dades) a una reixa qualsevol (simula_2d.py), per veí més proper
    des de la Release sols-cat. Torna None si no s'hi pot accedir."""
    import rasterio
    from rasterio.warp import reproject, Resampling
    fitxer, banda = VARIANTS[variant]
    try:
        with rasterio.open(REL + fitxer) as src:
            out = np.zeros(forma, np.uint8)
            reproject(rasterio.band(src, banda), out, dst_transform=tr, dst_crs=crs, src_nodata=0,
                      dst_nodata=0, resampling=Resampling.nearest)
    except Exception as e:
        print(f"Grups de sòl: {type(e).__name__} {e}")
        return None
    return out.astype(np.int8) - 1


def igme_poligons():
    """Polígons del Mapa de Permeabilitats 1:200.000 (IGME, ArcGIS REST) sobre Catalunya, en EPSG:25831, a trossos de
    20 km i paginats. Torna (llista de (geometria, valor del camp de permeabilitat), nom del camp)."""
    import requests
    from shapely.geometry import shape
    meta = requests.get(IGME, params={"f": "json"}, timeout=120).json()
    camps = [f["name"] for f in meta.get("fields", [])]
    print("Camps IGME:", camps)
    camp = next((c for c in camps if "perm" in c.lower() and "cod" not in c.lower()), None) or \
        next(c for c in camps if "perm" in c.lower())
    maxim = int(meta.get("maxRecordCount", 1000))
    x0, y0, x1, y1 = CAT
    vistos, out = set(), []
    for ya in range(y0, y1, 20000):
        for xa in range(x0, x1, 20000):
            off = 0
            while True:
                p = {"f": "geojson", "where": "1=1", "geometry": f"{xa},{ya},{xa + 20000},{ya + 20000}",
                     "geometryType": "esriGeometryEnvelope", "inSR": 25831, "outSR": 25831,
                     "spatialRel": "esriSpatialRelIntersects", "outFields": f"OBJECTID,{camp}",
                     "returnGeometry": "true", "resultOffset": off, "resultRecordCount": maxim}
                for intent in range(5):
                    try:
                        j = requests.get(IGME + "/query", params=p, timeout=300).json(); break
                    except Exception as e:
                        print(f"IGME {xa},{ya}: {type(e).__name__}"); time.sleep(10 * (intent + 1))
                else:
                    raise RuntimeError(f"IGME {xa},{ya} no respon")
                fs = j.get("features", [])
                for f in fs:
                    oid = f["properties"].get("OBJECTID", f.get("id"))
                    if oid in vistos or not f.get("geometry"):
                        continue
                    vistos.add(oid); out.append((shape(f["geometry"]), f["properties"].get(camp)))
                if len(fs) < maxim and not j.get("exceededTransferLimit"):
                    break
                off += len(fs)
    return out, camp


def litologia():
    import rasterio, unicodedata
    from collections import Counter
    from rasterio.features import rasterize
    from rasterio.transform import from_origin
    TMP.mkdir(parents=True, exist_ok=True)
    cfg = json.loads(LIT.read_text())
    pol, camp = igme_poligons()
    norm = lambda t: unicodedata.normalize("NFKD", str(t)).encode("ascii", "ignore").decode().strip().lower()
    vals = Counter(norm(v) for _, v in pol)
    print(f"IGME: {len(pol)} polígons; valors de {camp}: {dict(vals)}")
    corr = {norm(k): "ABCD".index(g) + 1 for k, g in cfg["correspondencia"].items()}
    sense = [v for v in vals if v not in corr]
    if sense:
        print("Valors sense correspondència (queden sense dades):", sense)
    x0, y0, x1, y1 = CAT
    tr = from_origin(x0, y1, PAS_LIT, PAS_LIT)
    forma = ((y1 - y0) // PAS_LIT, (x1 - x0) // PAS_LIT)
    g = rasterize([(geo, corr[norm(v)]) for geo, v in pol if norm(v) in corr], out_shape=forma, transform=tr,
                  fill=0, dtype="uint8")
    terra = dins_catalunya(tr, forma)
    g = omple(g, terra) * terra
    perfil = dict(driver="GTiff", width=forma[1], height=forma[0], crs=CRS, transform=tr, compress="deflate")
    with rasterio.open(TMP / "grups_sol_litologia_cat.tif", "w", count=1, dtype="uint8", nodata=0, **perfil) as f:
        f.write(g, 1)
        f.set_band_description(1, "litologia: permeabilitat IGME 1:200.000 → grup")
        f.update_tags(valors="1=A 2=B 3=C 4=D 0=sense dades", font=cfg["font"], correspondencia=json.dumps(cfg["correspondencia"], ensure_ascii=False))
    v = g[g > 0]
    print("litologia: " + ", ".join(f"{'ABCD'[k - 1]} {(v == k).mean() * 100:.0f} %" for k in range(1, 5)))
    web_variant(g, tr, forma, "grups_lit", "litologia", cfg["font"])


def web_variant(g, tr, forma, png, clau, font):
    """Afegeix una variant a la capa web: PNG en EPSG:3857 i % per conca a resum.json."""
    from rasterio.warp import reproject, Resampling, calculate_default_transform
    from rasterio.features import rasterize
    from shapely.geometry import shape
    from shapely.ops import transform as stransform
    from pyproj import Transformer
    from PIL import Image
    SORT.mkdir(parents=True, exist_ok=True)
    x0, y0, x1, y1 = CAT
    t3, w3, h3 = calculate_default_transform(CRS, "EPSG:3857", forma[1], forma[0], x0, y0, x1, y1)
    esc = max(1, int(np.ceil(max(w3, h3) / 2400)))              # PNG de com a molt ~2400 px
    t3, w3, h3 = t3 * t3.scale(esc, esc), w3 // esc, h3 // esc
    o = np.zeros((h3, w3), np.uint8)
    reproject(g, o, src_transform=tr, src_crs=CRS, dst_transform=t3, dst_crs="EPSG:3857", src_nodata=0, dst_nodata=0,
              resampling=Resampling.mode)
    im = Image.fromarray(o, "P")
    pal = [0, 0, 0] + sum((list(COLORS[k]) for k in range(1, 5)), [])
    im.putpalette(pal + [0] * (768 - len(pal)))
    im.save(SORT / f"{png}.png", optimize=True, transparency=bytes([0, 200, 200, 200, 200]))
    from rasterio.warp import transform_bounds
    lo0, la0, lo1, la1 = transform_bounds("EPSG:3857", 4326, t3.c, t3.f + h3 * t3.e, t3.c + w3 * t3.a, t3.f)
    rp = SORT / "resum.json"
    resum = json.loads(rp.read_text()) if rp.exists() else {"conques": {}, "colors": {"ABCD"[k - 1]: list(COLORS[k]) for k in range(1, 5)}}
    resum.setdefault("variants", {})[clau] = {"png": f"{png}.png", "limits": [[la0, lo0], [la1, lo1]], "font": font}
    t = Transformer.from_crs(4326, CRS, always_xy=True).transform
    for f in json.loads(CONQUES.read_text())["features"]:
        m = rasterize([(stransform(t, shape(f["geometry"])), 1)], out_shape=forma, transform=tr, fill=0,
                      dtype="uint8").astype(bool)
        v = g[m & (g > 0)]
        resum["conques"].setdefault("c%03d" % f["properties"]["id"], {})[clau] = \
            [round(float((v == k).mean() * 100), 1) for k in range(1, 5)] if v.size else None
    rp.write_text(json.dumps(resum, ensure_ascii=False, separators=(",", ":")))


def prova():
    import rasterio
    for p, c in (("clay", "0-5cm"),):
        try:
            with rasterio.open(SG.format(p=p, c=c)) as s:
                print("VRT ok", s.crs, s.res)
        except Exception as e:
            print("VRT", type(e).__name__, e)
    for u in ROCA:
        try:
            with rasterio.open(u) as s:
                print("Roca ok", u, s.crs, s.res, s.profile.get("tiled"), s.block_shapes)
        except Exception as e:
            print("Roca", u, type(e).__name__, e)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ordre", choices=["construeix", "litologia", "prova"])
    a = ap.parse_args()
    {"construeix": construeix, "litologia": litologia, "prova": prova}[a.ordre]()


if __name__ == "__main__":
    main()
