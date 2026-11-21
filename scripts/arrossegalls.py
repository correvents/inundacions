#!/usr/bin/env python3
"""Risc d'arrossegalls (troncs, canyes, branques, restes cremades, residus urbans) que poden tapar els passos de
torrentets sota vies (scripts/passos_cat.py) i els ponts, conca per conca a tot Catalunya. És una OPCIÓ del model:
es precalcula aquí i el model l'aplica només si s'activa, per comparar amb episodis si ajuda o no.

Mètode: el d'Australian Rainfall and Runoff 2019 (ARR, llibre 6, cap. 6, «Blockage of hydraulic structures»), que
classifica el potencial d'arrossegalls d'una obra amb tres factors de la seva conca i en treu un bloqueig de disseny:

1. Disponibilitat (quant material hi ha a prop de la llera): a la franja de ±FRANJA m dels torrentets aigües amunt
   de l'obra, fracció de
   - llenyós: boscos densos i clars (cobertes ICGC 2024, classes 7-9, 11-13) i bosc de ribera (15),
   - canyar probable: el mapa de cobertes no té classe de canyar; s'aproxima amb zones humides (20) i prats i
     herbassars (14) a tocar de la llera, que és on creix la canya,
   - cremat: zones cremades de les cobertes (17) i perímetres d'incendis de dades obertes dels últims ANYS_FOC anys
     (els incendis deixen fusta morta i sòl nu: més arrossegalls durant anys),
   - urbà: zones construïdes (21-35): cotxes, contenidors, tanques.
   Pes de cada tipus: llenyós i cremat 1, canyar 1, urbà 0,5 (ARR considera els residus urbans però més petits).
2. Mobilitat (quant fàcil és que es mogui): pendent mitjà de la franja.
3. Transportabilitat (si la riera el pot portar fins a l'obra): pendent mitjà del recorregut principal (com a
   5.2-IC, Témez) i conca drenada.
Cada factor es classifica Alt / Mitjà / Baix pels terços de tot Catalunya (no es tria cap llindar a mà); el
potencial és la suma (A=3, M=2, B=1): ≥ 8 Alt, 6-7 Mitjà, ≤ 5 Baix (a ARR: AAA i AAM = Alt; MMM i AMB = Mitjà;
MMB i BBB = Baix).

Bloqueig de disseny (ARR, taula de bloqueig d'entrada; L10 = llargada del 10 % d'arrossegalls més llargs, 1,5 m per
defecte, W = amplada de l'obertura):
                       Alt   Mitjà  Baix
   W < L10            100 %   50 %   25 %
   L10 ≤ W ≤ 3·L10     20 %   10 %    0 %
   W > 3·L10           10 %    0 %    0 %
- Ponts: W = llargada del pont a OSM.
- Passos (terraplè): l'obra no es veu a l'MDT ni sol ser a OSM; es desen els dos casos, obertura petita (tub
  < 1,5 m, el més habitual en camins i carreteres antigues) i mitjana, i el model tria.
- Guals: l'aigua ja passa per sobre; el bloqueig només empitjora la via (no es fa servir).
Escenari que en surt per al model (scripts/passos_cat.py): bloqueig ≥ 75 % → tapat, 25-75 % → mig tapat, < 25 % → obert.

Entrades: MDT 5 m (Release mdt5-cat), cobertes (Release cobertes-cat), passos (data/catalunya/passos/cNNN.geojson),
vies d'OSM (passos_cat.py osm), perímetres d'incendis (ordre `incendis`).
Sortides: data/catalunya/arrossegalls/cNNN.geojson (un punt per pas i pont) + index.json.

Ús:  python3 scripts/arrossegalls.py incendis --sortida /tmp/incendis.geojson
     python3 scripts/arrossegalls.py conca --id 125 --mdt mdt5_c125.tif --cobertes cobertes5_c125.tif \\
             --vies /tmp/vies.json.gz --incendis /tmp/incendis.geojson
     python3 scripts/arrossegalls.py index
Al model: from arrossegalls import escenari; escenari(pas, "petita") → "tapat" | "mig_tapat" | "obert"."""
import argparse, gzip, io, json, re, sys, tempfile, zipfile
from pathlib import Path
import numpy as np

CONQUES = Path("data/catalunya/conques/conques_100km2.geojson")
PASSOS = Path("data/catalunya/passos")
SORT = Path("data/catalunya/arrossegalls")
CAT_LL = (0.15, 40.5, 3.35, 42.9)                     # caixa de Catalunya (lon/lat)
FRANJA = 10                                           # m a cada costat de la llera (2 cel·les de l'MDT 5 m)
CONCA_MIN = 0.05                                      # km²: torrentet (el mateix que passos_cat.py)
ANYS_FOC = 10
LLENYOS = (7, 8, 9, 11, 12, 13, 15)
CANYAR = (14, 20)
CREMAT = (17,)
URBA = tuple(range(21, 36))
PES = {"llenyos": 1.0, "canyar": 1.0, "cremat": 1.0, "urba": 0.5}
L10 = 1.5
BLOQ = {"petita": {"A": 1.0, "M": 0.5, "B": 0.25}, "mitjana": {"A": 0.2, "M": 0.1, "B": 0.0},
        "gran": {"A": 0.1, "M": 0.0, "B": 0.0}}


# ---------- perímetres d'incendis ----------
def incendis(sortida):
    """Perímetres d'incendis forestals de Catalunya. Prova, per ordre: dades obertes de la Generalitat (cerca al
    catàleg, conjunts amb geometria), i l'arxiu d'àrees cremades d'EFFIS (Copernicus, ≥ 30 ha, des del 2000).
    Desa un GeoJSON (lon/lat) amb la propietat `any`."""
    import requests
    fts, fonts = [], []
    # 1. Generalitat (Socrata): conjunts d'incendis amb geometria, o «Base cartogràfica» amb fitxers adjunts (SHP)
    try:
        cat = requests.get("https://api.eu.socrata.com/api/catalog/v1", timeout=60,
                           params={"domains": "analisi.transparenciacatalunya.cat", "q": "incendis forestals", "limit": 50}).json()
        for r in cat.get("results", []):
            res = r["resource"]
            cols = [c.lower() for c in res.get("columns_datatype", [])]
            print(f"  catàleg: {res['id']} «{res['name']}» {set(cols)}")
            n0 = len(fts)
            if any(c in ("polygon", "multipolygon") for c in cols):
                g = requests.get(f"https://analisi.transparenciacatalunya.cat/resource/{res['id']}.geojson",
                                 params={"$limit": 200000}, timeout=300).json()
                for f in g.get("features", []):
                    a = _any(f.get("properties") or {})
                    if f.get("geometry") and a:
                        fts.append({"type": "Feature", "geometry": f["geometry"], "properties": {"any": a, "font": res["id"]}})
            elif "cartogr" in res["name"].lower() or "perímetr" in res["name"].lower():
                v = requests.get(f"https://analisi.transparenciacatalunya.cat/api/views/{res['id']}.json", timeout=60).json()
                urls = set(re.findall(r'https?://[^"\\\s]+', json.dumps(v)))
                urls |= {f"https://analisi.transparenciacatalunya.cat/api/views/{res['id']}/files/{x['assetId']}?filename={x.get('filename', '')}"
                         for x in (v.get("metadata", {}).get("attachments") or []) if x.get("assetId")}
                for u in sorted(urls):
                    print(f"    enllaç: {u[:200]}")
                for u in [u for u in urls if "format-shp" in u]:         # pàgina d'Agricultura amb un zip per any
                    try:
                        html = requests.get(u.replace("http://", "https://"), timeout=120).text
                        from urllib.parse import urljoin
                        urls |= {urljoin(u, h) for h in re.findall(r'href="([^"]+\.zip)"', html, re.I)}
                    except Exception as e:                               # noqa: BLE001
                        print(f"    {u}: {e}")
                for u in sorted(urls):
                    if not re.search(r"\.zip|/files/|shp|download", u, re.I):
                        continue
                    try:
                        c = requests.get(u, timeout=600).content
                        if c[:2] == b"PK":
                            k = _shp_zip(c, fts, f"Generalitat {res['id']}", Path(u).name)
                            print(f"    {u[:120]} → {k} perímetres")
                    except Exception as e:                               # noqa: BLE001
                        print(f"    {u[:120]}: {e}")
            if len(fts) > n0:
                fonts.append(f"Generalitat {res['id']} ({len(fts) - n0})")
    except Exception as e:                                               # noqa: BLE001
        print(f"  Generalitat: {e}")
    # 2. EFFIS (WFS, Shapefile en zip, ≥ 30 ha): completa els anys que la Generalitat encara no ha publicat
    if True:
        for bbox in (f"{CAT_LL[0]},{CAT_LL[1]},{CAT_LL[2]},{CAT_LL[3]}", f"{CAT_LL[1]},{CAT_LL[0]},{CAT_LL[3]},{CAT_LL[2]}",
                     f"{CAT_LL[0]},{CAT_LL[1]},{CAT_LL[2]},{CAT_LL[3]},EPSG:4326"):
            try:
                r = requests.get("https://maps.effis.emergency.copernicus.eu/effis", timeout=600, params={
                    "service": "WFS", "request": "getfeature", "typename": "ms:modis.ba.poly", "version": "1.1.0",
                    "outputformat": "SHAPEZIP", "bbox": bbox})
                print(f"  EFFIS bbox {bbox}: {r.status_code} {len(r.content)} B {r.content[:80]!r}")
                if r.content[:2] == b"PK":
                    k = _shp_zip(r.content, fts, "EFFIS", "")
                    print(f"    → {k} perímetres")
                    if k:
                        fonts.append(f"EFFIS ({k})"); break
            except Exception as e:                                       # noqa: BLE001
                print(f"  EFFIS: {e}")
    Path(sortida).write_text(json.dumps({"type": "FeatureCollection", "fonts": fonts, "features": fts}, separators=(",", ":")))
    print(f"Incendis: {len(fts)} perímetres; fonts: {fonts or 'cap (només les zones cremades de les cobertes)'}")


def _any(p, nom=""):
    p = {k.lower(): v for k, v in p.items()}
    for k in ("any", "any_incen", "year", "data_inici", "data_incen", "data", "dat_inici", "data_incendi", "firedate", "fecha"):
        v = str(p.get(k) or "")
        m = re.search(r"(19[6-9]\d|20[0-4]\d)", v)
        if m:
            return int(m.group(1))
    m = re.search(r"(19[6-9]\d|20[0-4]\d)", nom)
    return int(m.group(1)) if m else None


def _shp_zip(contingut, fts, font, nom):
    """Llegeix tots els Shapefiles d'un zip (també zips dins de zips) i hi afegeix els polígons de Catalunya."""
    import shapefile
    from pyproj import CRS, Transformer
    from shapely.geometry import shape, mapping
    from shapely.ops import transform as stransform
    z = zipfile.ZipFile(io.BytesIO(contingut)); n0 = len(fts); caixa = _caixa()
    for n in z.namelist():
        if n.lower().endswith(".zip"):
            _shp_zip(z.read(n), fts, font, n)
        if not n.lower().endswith(".shp"):
            continue
        base = n[:-4]
        noms = {x.lower(): x for x in z.namelist()}
        prj = z.read(noms[(base + ".prj").lower()]).decode(errors="ignore") if (base + ".prj").lower() in noms else ""
        sf = shapefile.Reader(shp=io.BytesIO(z.read(n)), dbf=io.BytesIO(z.read(noms[(base + ".dbf").lower()])), encoding="latin-1")
        camps = [c[0] for c in sf.fields[1:]]
        crs = CRS.from_wkt(prj) if prj else CRS.from_epsg(25831)
        a_ll = Transformer.from_crs(crs, 4326, always_xy=True).transform
        print(f"      {n}: {len(sf)} elements, camps {camps[:12]}, {crs.name}")
        for sr in sf.iterShapeRecords():
            a = _any(dict(zip(camps, sr.record)), n + nom)
            if not a or sr.shape.shapeType not in (5, 15, 25):
                continue
            g = stransform(lambda x, y, *_: a_ll(x, y), shape(sr.shape.__geo_interface__))  # noqa
            if g.is_valid and g.intersects(caixa):
                fts.append({"type": "Feature", "geometry": mapping(g.simplify(0.0002)), "properties": {"any": a, "font": font}})
    return len(fts) - n0


def _caixa():
    from shapely.geometry import box
    return box(*CAT_LL)


# ---------- acumulació ponderada aigües avall ----------
def acumula(ordre, fd, pesos):
    """Suma cada capa de `pesos` (k, nf, nc) de totes les cel·les aigües amunt (incloent-hi la mateixa)."""
    import numba

    @numba.njit(cache=True)
    def passa(ordre, fd, w, nf, nc):
        dr = {64: -1, 128: -1, 1: 0, 2: 1, 4: 1, 8: 1, 16: 0, 32: -1}
        dc = {64: 0, 128: 1, 1: 1, 2: 1, 4: 0, 8: -1, 16: -1, 32: -1}
        for k in ordre:
            r = k // nc; c = k % nc; d = fd[r, c]
            if d not in dr:
                continue
            r2 = r + dr[d]; c2 = c + dc[d]
            if r2 < 0 or r2 >= nf or c2 < 0 or c2 >= nc:
                continue
            for i in range(w.shape[0]):
                w[i, r2, c2] += w[i, r, c]

    w = pesos.astype("float64").copy()
    passa(ordre, fd, w, fd.shape[0], fd.shape[1])
    return w


# ---------- una conca ----------
def conca(cid, cami_mdt, cami_cob, vies_gz, cami_inc, any_ref):
    import rasterio
    from rasterio import features
    from pyproj import Transformer
    from scipy import ndimage
    from shapely.geometry import shape, mapping, LineString, Point
    from shapely.ops import transform as stransform
    sys.path.insert(0, str(Path(__file__).parent))
    np.in1d = np.isin
    from passos_cat import recorreguts

    C = {f["properties"]["id"]: f for f in json.loads(CONQUES.read_text())["features"]}
    geom = shape(C[cid]["geometry"])
    with rasterio.open(cami_mdt) as f:
        z = f.read(1).astype("float32"); nd = f.nodata; tr = f.transform; crs = f.crs
    z[z == nd] = np.nan; z /= 100
    z[z <= 0.2] = np.nan
    with rasterio.open(cami_cob) as f:
        cob = f.read(1)
    assert cob.shape == z.shape, "cobertes i MDT han de tenir la mateixa reixa"
    nf, nc = z.shape; cel = abs(tr.a * tr.e)
    a_utm = Transformer.from_crs(4326, crs, always_xy=True).transform
    a_ll = Transformer.from_crs(crs, 4326, always_xy=True).transform
    w, s, e, n = geom.bounds
    print(f"== Conca {cid}: {C[cid]['properties']['km2']} km², reixa {nc}×{nf}")

    # drenatge (igual que passos_cat.py, perquè els passos caiguin a les mateixes cel·les)
    from pysheds.grid import Grid
    tmp = Path(tempfile.mkdtemp()) / "z.tif"
    with rasterio.open(tmp, "w", driver="GTiff", height=nf, width=nc, count=1, dtype="float32", crs=crs, transform=tr,
                       nodata=np.nan) as f:
        f.write(z, 1)
    gr = Grid.from_raster(str(tmp)); dem = gr.read_raster(str(tmp))
    ple = gr.fill_depressions(gr.fill_pits(dem))
    fdr = gr.flowdir(gr.resolve_flats(ple))
    acc = np.asarray(gr.accumulation(fdr), dtype="float64")
    fd = np.asarray(fdr).astype("int64")
    del ple
    km2 = acc * cel / 1e6
    ordre = np.argsort(acc, axis=None, kind="stable")
    lup, zmax = recorreguts(fd, acc, z)

    # franja de la llera i pendent
    llera = km2 >= CONCA_MIN
    franja = ndimage.binary_dilation(llera, iterations=max(int(round(FRANJA / tr.a)), 1)) & np.isfinite(z)
    gy, gx = np.gradient(np.where(np.isfinite(z), z, np.nanmean(z)), tr.a)
    pend = np.hypot(gx, gy).astype("float32")

    # incendis recents rasteritzats
    cremat = np.isin(cob, CREMAT)
    n_foc = 0
    if cami_inc and Path(cami_inc).exists():
        inc = json.loads(Path(cami_inc).read_text())["features"]
        caixa = shape({"type": "Polygon", "coordinates": [[[w, s], [e, s], [e, n], [w, n], [w, s]]]})
        gs = []
        for f in inc:
            if any_ref - f["properties"]["any"] > ANYS_FOC:
                continue
            g = shape(f["geometry"])
            if g.intersects(caixa):
                gs.append(mapping(stransform(lambda x, y: a_utm(x, y), g)))  # noqa
        if gs:
            cremat |= features.rasterize(((g, 1) for g in gs), out_shape=z.shape, transform=tr, dtype="uint8") > 0
        n_foc = len(gs)
    capes = {"franja": franja, "llenyos": franja & np.isin(cob, LLENYOS), "canyar": franja & np.isin(cob, CANYAR),
             "cremat": franja & cremat, "urba": franja & np.isin(cob, URBA)}
    pes = np.stack([capes[k].astype("float32") for k in capes] + [np.where(franja, pend, 0)])
    W = acumula(ordre, fd, pes)
    noms = list(capes) + ["pend_franja"]
    print(f"Franja: {franja.sum() * cel / 1e6:.2f} km²; incendis ≤ {ANYS_FOC} anys a la caixa: {n_foc}")

    def mesura(r, c):
        v = {k: float(W[i, r, c]) for i, k in enumerate(noms)}
        fr = max(v["franja"], 1.0)
        l_km = float(lup[r, c]) * tr.a / 1000
        d = {"conca_km2": round(float(km2[r, c]), 3),
             "franja_ha": round(v["franja"] * cel / 1e4, 2),
             **{f"f_{k}": round(v[k] / fr, 3) for k in PES},
             "pendent_franja": round(v["pend_franja"] / fr, 3),
             "pendent_riera": round((float(zmax[r, c]) - float(z[r, c])) / max(l_km * 1000, 1), 4)}
        d["disponibilitat"] = round(sum(PES[k] * d[f"f_{k}"] for k in PES), 3)
        return d

    def cella(x, y, radi=4, km2_pas=None):
        """Cel·la on el torrentet travessa la via: la de més acumulació a prop del punt o, si passos_cat.py en dona
        la conca, la que hi té la conca més semblant (el punt del pas és al fons del llac, no a la sortida)."""
        cc, rr = ~tr * (x, y); rr, cc = int(rr), int(cc)
        r0, r1, c0, c1 = max(rr - radi, 0), min(rr + radi + 1, nf), max(cc - radi, 0), min(cc + radi + 1, nc)
        if r0 >= r1 or c0 >= c1:
            return None
        sub = acc[r0:r1, c0:c1] if km2_pas is None else -np.abs(np.log((km2[r0:r1, c0:c1] + 1e-6) / km2_pas))
        rs, cs = np.unravel_index(np.argmax(sub), sub.shape)
        return r0 + rs, c0 + cs

    fts = []
    # passos de passos_cat.py
    pj = PASSOS / f"c{cid:03d}.geojson"
    passos = [f for f in json.loads(pj.read_text())["features"] if f["geometry"]["type"] == "Point"] if pj.exists() else []
    for f in passos:
        p = f["properties"]
        t = p["tipus"] == "terraplè"
        rc = cella(*a_utm(*f["geometry"]["coordinates"]), radi=12 if t else 4, km2_pas=p["conca_km2"] if t else None)
        if rc is None:
            continue
        fts.append({"type": "Feature", "geometry": f["geometry"],
                    "properties": {"id": p["id"], "tipus": p["tipus"], "via": p.get("via"), "nom": p.get("nom"), **mesura(*rc)}})
    # ponts d'OSM sobre torrentets
    n_ponts = 0
    if vies_gz:
        with gzip.open(vies_gz, "rt") as fh:
            osmj = json.load(fh)
        for vid, tipus, nom, fl, xy in osmj["vies"]:
            if fl != "p":
                continue
            xs, ys = xy[0::2], xy[1::2]
            if max(xs) < w or min(xs) > e or max(ys) < s or min(ys) > n:
                continue
            g = LineString([a_utm(xy[k], xy[k + 1]) for k in range(0, len(xy), 2)])
            m = features.rasterize([(mapping(g.buffer(tr.a)), 1)], out_shape=z.shape, transform=tr, dtype="uint8") > 0
            m &= llera
            if not m.any():
                continue
            r, c = np.unravel_index(np.argmax(np.where(m, acc, -1)), m.shape)
            x, y = tr * (c + 0.5, r + 0.5)
            n_ponts += 1
            fts.append({"type": "Feature", "geometry": {"type": "Point", "coordinates": [round(q, 5) for q in a_ll(x, y)]},
                        "properties": {"id": f"{cid}-p{n_ponts}", "tipus": "pont", "via": tipus, "nom": nom, "osm_via": vid,
                                       "llum_m": round(g.length, 1), **mesura(r, c)}})
    # un pont pot tenir diversos trams OSM (calçades): es queda el de més llum per cel·la
    vist = {}
    for f in fts:
        if f["properties"]["tipus"] == "pont":
            k = tuple(f["geometry"]["coordinates"])
            if k not in vist or f["properties"]["llum_m"] > vist[k]["properties"]["llum_m"]:
                vist[k] = f
    fts = [f for f in fts if f["properties"]["tipus"] != "pont"] + list(vist.values())
    res = {"km2": C[cid]["properties"]["km2"], "passos": sum(f["properties"]["tipus"] == "terraplè" for f in fts),
           "guals": sum(f["properties"]["tipus"] == "gual" for f in fts), "ponts": len(vist), "incendis": n_foc,
           "franja_km2": round(float(franja.sum()) * cel / 1e6, 2)}
    print(json.dumps(res, ensure_ascii=False))
    SORT.mkdir(parents=True, exist_ok=True)
    (SORT / f"c{cid:03d}.geojson").write_text(json.dumps({"type": "FeatureCollection", "resum": res, "features": fts},
                                                         ensure_ascii=False, separators=(",", ":")))


# ---------- classes per terços de tot Catalunya ----------
def classe(v, t1, t2):
    return "B" if v <= t1 else ("M" if v <= t2 else "A")


def index():
    fitxers = sorted(SORT.glob("c[0-9]*.geojson"))
    tots = [json.loads(f.read_text()) for f in fitxers]
    obres = [f["properties"] for d in tots for f in d["features"] if f["properties"]["tipus"] != "gual"]
    if not obres:
        print("Cap obra"); return
    factors = {"disponibilitat": lambda p: p["disponibilitat"], "mobilitat": lambda p: p["pendent_franja"],
               "transportabilitat": lambda p: p["pendent_riera"] * np.sqrt(p["conca_km2"])}
    talls = {k: [float(np.quantile([fn(p) for p in obres], q)) for q in (1 / 3, 2 / 3)] for k, fn in factors.items()}
    idx, tot = {}, {"A": 0, "M": 0, "B": 0}
    for f, d in zip(fitxers, tots):
        for ft in d["features"]:
            p = ft["properties"]
            cl = {k: classe(fn(p), *talls[k]) for k, fn in factors.items()}
            punts = sum({"A": 3, "M": 2, "B": 1}[c] for c in cl.values())
            pot = "A" if punts >= 8 else ("M" if punts >= 6 else "B")
            p["classes"] = "".join(cl[k] for k in factors); p["potencial"] = pot
            if p["tipus"] == "pont":
                ob = "petita" if p["llum_m"] < L10 else ("mitjana" if p["llum_m"] <= 3 * L10 else "gran")
                p["bloqueig"] = {ob: BLOQ[ob][pot]}
            elif p["tipus"] == "terraplè":
                p["bloqueig"] = {"petita": BLOQ["petita"][pot], "mitjana": BLOQ["mitjana"][pot]}
            if p["tipus"] != "gual":
                tot[pot] += 1
        d["resum"]["potencial"] = {c: sum(ft["properties"].get("potencial") == c and ft["properties"]["tipus"] != "gual"
                                          for ft in d["features"]) for c in "AMB"}
        f.write_text(json.dumps(d, ensure_ascii=False, separators=(",", ":")))
        idx[int(f.stem[1:])] = d["resum"]
    (SORT / "index.json").write_text(json.dumps({"metode": "ARR 2019 llibre 6 cap. 6", "talls_terços": talls,
                                                 "conques": idx, "total": tot}, ensure_ascii=False, separators=(",", ":")))
    print(f"{len(idx)} conques, {len(obres)} obres; potencial {tot}; talls {talls}")


def escenari(pas, obertura="petita"):
    """Escenari d'un pas (properties de data/catalunya/arrossegalls) per a l'opció d'arrossegalls del model."""
    b = (pas.get("bloqueig") or {}).get(obertura)
    if b is None:
        b = next(iter((pas.get("bloqueig") or {"x": 0.0}).values()))
    return "tapat" if b >= 0.75 else ("mig_tapat" if b >= 0.25 else "obert")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="ordre", required=True)
    i = sub.add_parser("incendis"); i.add_argument("--sortida", required=True)
    c = sub.add_parser("conca"); c.add_argument("--id", type=int, required=True); c.add_argument("--mdt", required=True)
    c.add_argument("--cobertes", required=True); c.add_argument("--vies"); c.add_argument("--incendis")
    c.add_argument("--any", type=int, default=2026, help="any de referència per als incendis recents")
    sub.add_parser("index")
    a = ap.parse_args()
    if a.ordre == "incendis":
        return incendis(a.sortida)
    if a.ordre == "index":
        return index()
    return conca(a.id, a.mdt, a.cobertes, a.vies, a.incendis, a.any)


if __name__ == "__main__":
    sys.exit(main())
