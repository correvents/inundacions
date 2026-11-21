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

Ús:  python3 scripts/sols.py construeix
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
URL = "https://github.com/correvents/inundacions/releases/download/sols-cat/grups_sol_cat.tif"
BANDES = {"mapa": 1, "mapa-sup": 2}
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


def a_reixa(tr, forma, crs, variant="mapa", url=URL):
    """Grup hidrològic per cel·la (0=A … 3=D; -1 sense dades) a una reixa qualsevol (simula_2d.py), per veí més proper
    des de la Release sols-cat. Torna None si no s'hi pot accedir."""
    import rasterio
    from rasterio.warp import reproject, Resampling
    try:
        with rasterio.open(url) as src:
            out = np.zeros(forma, np.uint8)
            reproject(rasterio.band(src, BANDES[variant]), out, dst_transform=tr, dst_crs=crs, src_nodata=0,
                      dst_nodata=0, resampling=Resampling.nearest)
    except Exception as e:
        print(f"Grups de sòl: {type(e).__name__} {e}")
        return None
    return out.astype(np.int8) - 1


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
    ap.add_argument("ordre", choices=["construeix", "prova"])
    a = ap.parse_args()
    {"construeix": construeix, "prova": prova}[a.ordre]()


if __name__ == "__main__":
    main()
