#!/usr/bin/env python3
"""Clavegueram pluvial estimat, conca per conca a tot Catalunya, com a OPCIÓ per a simula_2d.py (escenaris funciona / antic / col·lapsat).

Ara la simulació fa servir una capacitat fixa (--clavegueram 30 mm/h) a totes les zones construïdes. Aquí se n'estima una per
cada zona urbana amb el mateix criteri amb què es dimensiona una xarxa (mètode racional), sense cap número triat per a cada poble:

1. Pluja de disseny (subordre `idf`): màxims anuals de la XEMA (Meteocat, dades obertes, des de 2009) de cada estació,
   en 30 min i en 24 h. Gumbel regional (índex: corba de creixement comuna de tot Catalunya × mitjana de cada estació),
   factor 1,13 d'interval fix (Weiss). D'aquí P24(T) i I30(T), i el quocient I1/Id de la fórmula de Témez de cada
   estació, que dona la intensitat I(T, t) per a qualsevol durada t.
2. Zones urbanes (subordre `conca`): cobertes del sòl de l'ICGC (Release cobertes-cat, classes urbanes de
   config/cobertes_parametres.json) sobre l'MDT 5 m (Release mdt5-cat), agregats a 10 m.
3. Sortides de la xarxa: seguint el drenatge de l'MDT des de cada cel·la urbana, el primer punt on l'aigua surt de la
   zona urbana, arriba a una llera (conca ≥ --llera-km2) o al mar. Totes les cel·les que hi van a parar formen una
   subconca de clavegueram amb un col·lector de sortida.
4. Temps de concentració de cada subconca: temps d'entrada (--t-entrada min) + recorregut més llarg / velocitat a la xarxa
   (--v-xarxa m/s), valors habituals de dimensionament.
5. Capacitat (mm/h sobre la superfície urbana) = C·I(T, tc), amb C el coeficient d'escorrentiu mitjà de la subconca segons
   les cobertes (mètode racional, ASCE/WEF 1992) i I la intensitat de la pluja de disseny a la sortida (IDW de les 4 estacions
   més properes): T = 10 anys (criteri habitual de les xarxes noves) i T = 2 anys (xarxes antigues o infradimensionades).
   L'escenari col·lapsat (embornals tapats per fulles i arrossegalls) s'aplica a simula_2d.py com a fracció de la de 10 anys.
6. OSM (subordre `osm`, opcional): embornals, pous i col·lectors cartografiats a OpenStreetMap. A Catalunya n'hi ha
   pocs; es compten per subconca per veure on n'hi ha prou per comprovar l'estimació. No canvien la capacitat.

Sortides:  /tmp/clavegueram/clav10_cNNN.tif (Release clavegueram-cat; EPSG de l'MDT, 10 m, uint32):
               banda 1 capacitat T10 (dècimes de mm/h), banda 2 capacitat T2, banda 3 número de sortida (0 = no urbà)
           data/catalunya/clavegueram/cNNN.json (sortides: x, y, tipus, ha, tc, capacitats, embornals OSM)
           data/catalunya/clavegueram/index.json, idf_estacions.json

Ús:  python3 scripts/clavegueram_cat.py idf
     python3 scripts/clavegueram_cat.py osm --pbf cataluna-latest.osm.pbf --sortida /tmp/clav_osm.json.gz
     python3 scripts/clavegueram_cat.py conca --id 188 --mdt mdt5_c188.tif --cobertes cobertes5_c188.tif [--osm /tmp/clav_osm.json.gz]
     python3 scripts/clavegueram_cat.py index"""
import argparse, gzip, json, math, sys, time
from pathlib import Path
import numpy as np

SORT = Path("data/catalunya/clavegueram")
IDF = SORT / "idf_estacions.json"
TMP = Path("/tmp/clavegueram")
PARAM = Path("config/cobertes_parametres.json")
API = "https://analisi.transparenciacatalunya.cat/resource/"
URL = "https://github.com/correvents/inundacions/releases/download/clavegueram-cat/clav10_{}.tif"
RETORNS = (2, 5, 10, 25)
# coeficient d'escorrentiu del mètode racional per coberta urbana (ASCE/WEF 1992, Design and construction of urban stormwater
# management systems, taula 5.5: valor central de cada interval). La xarxa es dimensiona per a Q = C·I·A.
C_RACIONAL = {21: 0.85, 22: 0.85, 23: 0.65, 25: 0.40, 27: 0.75, 31: 0.85, 33: 0.75, 35: 0.75}
C_ALTRES = 0.25                                                        # parcs i solars dins la trama (gespa, 0,10-0,35)
EXP = lambda t_h: (28 ** 0.1 - t_h ** 0.1) / (28 ** 0.1 - 1)          # exponent de la corba IDF de Témez (t en hores)


# ---------- 1. pluja de disseny amb la XEMA ----------
def socrata(ds, params):
    import requests
    for i in range(5):
        try:
            r = requests.get(API + ds + ".json", params=params, timeout=300); r.raise_for_status(); return r.json()
        except Exception as e:
            print("reintent", ds, type(e).__name__, str(e)[:200]); time.sleep(15 * (i + 1))
    return None


def gumbel_k(norm):
    """Factors de creixement K(T) d'una mostra normalitzada (mitjana 1) pel mètode dels moments."""
    x = np.array(norm); beta = x.std(ddof=1) * math.sqrt(6) / math.pi; mu = x.mean() - 0.5772 * beta
    return {T: float(mu - beta * math.log(-math.log(1 - 1 / T))) for T in RETORNS}


def idf():
    est = socrata("yqwd-vj5e", {"$limit": 5000}) or []
    print(f"XEMA: {len(est)} estacions")
    maxims = {}
    for e in est:
        codi = e["codi_estacio"]
        r = socrata("nzvn-apee", {"$select": "date_trunc_ymd(data_lectura) as d, sum(valor_lectura) as s, count(valor_lectura) as n, max(valor_lectura) as m",
                                  "$where": f"codi_estacio='{codi}' AND codi_variable='35'", "$group": "d", "$limit": 20000})
        if not r:
            continue
        anys = {}
        for x in r:
            n = int(x.get("n", 0))
            if n < 44:                                                   # dies amb ≥ 22 h de dades
                continue
            a = anys.setdefault(x["d"][:4], [0.0, 0.0, 0])
            a[0] = max(a[0], float(x["s"])); a[1] = max(a[1], float(x["m"])); a[2] += 1
        bons = {y: v for y, v in anys.items() if v[2] >= 330}               # anys gairebé complets
        if len(bons) >= 5:
            maxims[codi] = {"nom": e.get("nom_estacio", codi), "lat": float(e["latitud"]), "lon": float(e["longitud"]),
                            "p24": [v[0] for v in bons.values()], "m30": [v[1] for v in bons.values()], "anys": sorted(bons)}
            print(f"{codi} {maxims[codi]['nom'][:30]:30s} {len(bons)} anys: P24 mitjana {np.mean(maxims[codi]['p24']):.0f} mm, 30 min {np.mean(maxims[codi]['m30']):.1f} mm")
    if len(maxims) < 10:
        raise SystemExit(f"Només {len(maxims)} estacions amb ≥ 5 anys complets")
    k24 = gumbel_k([v / np.mean(m["p24"]) for m in maxims.values() for v in m["p24"]])
    k30 = gumbel_k([v / np.mean(m["m30"]) for m in maxims.values() for v in m["m30"] if np.mean(m["m30"]) > 0])
    out = []
    for codi, m in maxims.items():
        p24m, m30m = float(np.mean(m["p24"])), float(np.mean(m["m30"]))
        fila = {"codi": codi, "nom": m["nom"], "lat": m["lat"], "lon": m["lon"], "anys": len(m["anys"])}
        for T in RETORNS:
            id_ = 1.13 * p24m * k24[T] / 24                               # mm/h mitjana en 24 h
            i30 = 1.13 * m30m * k30[T] * 2                                # mm/h en 30 min
            fila[f"p24_t{T}"] = round(id_ * 24, 1); fila[f"i30_t{T}"] = round(i30, 1)
            fila[f"i1id_t{T}"] = round(float(np.clip(math.exp(math.log(max(i30 / id_, 1.01)) / EXP(0.5)), 6, 20)), 2)
        out.append(fila)
    SORT.mkdir(parents=True, exist_ok=True)
    IDF.write_text(json.dumps({"font": "XEMA, Meteocat (dades obertes, nzvn-apee, variable 35)", "metode":
        "Gumbel regional (corba de creixement comuna × mitjana de l'estació), 1,13 interval fix, IDF de Témez amb I1/Id de cada estació",
        "k_p24": k24, "k_30min": k30, "estacions": out}, ensure_ascii=False, indent=1))
    print(f"IDF: {len(out)} estacions; K(10) P24 {k24[10]:.2f}, 30 min {k30[10]:.2f}")


def intensitat(idf_est, lon, lat, T, t_h):
    """I(T, t) en mm/h a (lon, lat): IDW (1/d²) de les 4 estacions més properes, cadascuna amb la seva IDF de Témez."""
    d = np.array([math.hypot((e["lat"] - lat) * 111, (e["lon"] - lon) * 111 * math.cos(math.radians(lat))) for e in idf_est])
    k = np.argsort(d)[:4]; w = 1 / np.maximum(d[k], 1) ** 2
    vals = [idf_est[i][f"p24_t{T}"] / 24 * idf_est[i][f"i1id_t{T}"] ** EXP(t_h) for i in k]
    return float(np.sum(w * vals) / w.sum())


# ---------- 6. OSM ----------
def osm(pbf, sortida):
    import osmium

    class H(osmium.SimpleHandler):
        def __init__(self):
            super().__init__(); self.punts = []; self.linies = []

        def node(self, n):
            t = n.tags
            mh, mm = t.get("manhole", ""), t.get("man_made", "")
            if mm == "storm_drain" or t.get("highway") == "storm_drain" or mh in ("drain", "storm", "gully") or \
               (mm == "manhole" and mh in ("sewer", "drain", "storm", "")) or t.get("drain") == "storm":
                self.punts.append([round(n.location.lon, 6), round(n.location.lat, 6), "embornal" if mh in ("drain", "gully") or "storm" in (mm + mh) else "pou"])

        def way(self, w):
            t = w.tags
            sub = (t.get("substance", "") + t.get("type", "")).lower()
            col = (t.get("man_made") == "pipeline" and ("sewage" in sub or "rainwater" in sub or "water" == sub and t.get("usage") == "drain")) or \
                  (t.get("waterway") in ("drain", "stream", "ditch", "canal", "river") and t.get("tunnel") in ("culvert", "yes", "flooded"))
            if col:
                try:
                    self.linies.append([[round(p.lon, 6), round(p.lat, 6)] for p in w.nodes])
                except osmium.InvalidLocationError:
                    pass

    h = H(); h.apply_file(pbf, locations=True)
    with gzip.open(sortida, "wt") as f:
        json.dump({"punts": h.punts, "colectors": h.linies}, f, separators=(",", ":"))
    print(f"OSM: {len(h.punts)} embornals/pous, {len(h.linies)} trams de col·lector o riera soterrada")


# ---------- 2-5. una conca ----------
def sortides_xarxa(fd, ordre, urba, llera, mar, fora_max, nf, nc):
    """Per a cada cel·la: índex lineal de la cel·la on desguassa la xarxa (-1 si surt de la reixa) i recorregut fins allà
    (cel·les, diagonal √2). La xarxa segueix el drenatge de l'MDT fins a una llera o el mar; pot travessar trams no urbans
    (parcs, solars) de fins a fora_max cel·les. Si n'ha de travessar més, l'aigua surt per la superfície al primer punt no
    urbà (sortida «terreny»). Es recorren les cel·les de més a menys acumulació: la de sota sempre es fa abans."""
    import numba

    @numba.njit(cache=True)
    def passa(ordre, fd, urba, llera, mar, fora_max, sort, dist, fr, nf, nc):
        for k in ordre:
            r = k // nc; c = k % nc
            if mar[r, c]:
                continue
            d = fd[r, c]
            if d == 64: r2, c2 = r - 1, c
            elif d == 128: r2, c2 = r - 1, c + 1
            elif d == 1: r2, c2 = r, c + 1
            elif d == 2: r2, c2 = r + 1, c + 1
            elif d == 4: r2, c2 = r + 1, c
            elif d == 8: r2, c2 = r + 1, c - 1
            elif d == 16: r2, c2 = r, c - 1
            elif d == 32: r2, c2 = r - 1, c - 1
            else:
                sort[r, c] = k; dist[r, c] = 0.5; continue                  # pou: desguassa a la mateixa cel·la
            pas = 1.4142135 if r2 != r and c2 != c else 1.0
            if r2 < 0 or r2 >= nf or c2 < 0 or c2 >= nc:
                sort[r, c] = -1; dist[r, c] = pas; fr[r, c] = 0.0 if urba[r, c] else pas; continue
            k2 = r2 * nc + c2
            if llera[r2, c2] or mar[r2, c2]:
                sort[r, c] = k2; dist[r, c] = pas; fr[r, c] = 0.0 if urba[r, c] else pas
            elif urba[r2, c2] or (urba[r, c] and fr[r2, c2] <= fora_max) or not urba[r, c]:
                sort[r, c] = sort[r2, c2]; dist[r, c] = dist[r2, c2] + pas
                fr[r, c] = 0.0 if urba[r, c] else pas + (0.0 if urba[r2, c2] else fr[r2, c2])
            else:                                                          # massa camp fins a la llera: surt per la superfície
                sort[r, c] = k2; dist[r, c] = pas; fr[r, c] = 0.0

    sort = np.full((nf, nc), -2, np.int64); dist = np.zeros((nf, nc), np.float32); fr = np.zeros((nf, nc), np.float32)
    passa(ordre, fd, urba, llera, mar, fora_max, sort, dist, fr, nf, nc)
    return sort, dist


def conca(cid, cami_mdt, cami_cob, cami_osm, a):
    import rasterio, tempfile
    from rasterio.warp import reproject, Resampling
    from pyproj import Transformer
    from scipy import ndimage
    np.in1d = np.isin                                                  # pysheds encara fa servir np.in1d (numpy 2)
    from pysheds.grid import Grid
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import cobertes as COB

    idf_est = json.loads(IDF.read_text())["estacions"]
    with rasterio.open(cami_mdt) as f:
        z5 = f.read(1).astype("float32"); nd = f.nodata; tr5 = f.transform; crs = f.crs
    z5[z5 == nd] = np.nan; z5 /= 100; z5[z5 <= 0.2] = np.nan                  # cm → m; mar
    fac = max(int(round(a.pas / tr5.a)), 1)
    nf, nc = z5.shape[0] // fac, z5.shape[1] // fac
    with np.errstate(all="ignore"):
        z = np.nanmean(z5[: nf * fac, : nc * fac].reshape(nf, fac, nc, fac), axis=(1, 3)).astype("float32")
    tr = tr5 * tr5.scale(fac, fac); dx = tr.a
    cob = np.zeros((nf, nc), np.uint8)
    with rasterio.open(cami_cob) as f:
        reproject(rasterio.band(f, 1), cob, dst_transform=tr, dst_crs=crs, src_nodata=0, dst_nodata=0, resampling=Resampling.mode)
    _, _, tipus = COB.parametres(cob)
    urba = (tipus == 1) & np.isfinite(z)
    urba = ndimage.binary_closing(urba, iterations=2) & np.isfinite(z)        # carrers i places petites que el mapa parteix
    print(f"== Conca {cid}: reixa {nc}×{nf} de {dx:.0f} m, urbà {urba.sum() * dx * dx / 1e6:.2f} km²")
    if urba.sum() == 0:
        return desa(cid, None, tr, crs, [], {"urba_km2": 0})

    tmp = Path(tempfile.mkdtemp()) / "z.tif"
    with rasterio.open(tmp, "w", driver="GTiff", height=nf, width=nc, count=1, dtype="float32", crs=crs, transform=tr, nodata=np.nan) as f:
        f.write(z, 1)
    gr = Grid.from_raster(str(tmp)); dem = gr.read_raster(str(tmp))
    fd = gr.flowdir(gr.resolve_flats(gr.fill_depressions(gr.fill_pits(dem))))
    acc = np.asarray(gr.accumulation(fd), dtype="float64")
    fd = np.asarray(fd).astype("int64")
    llera = acc * dx * dx / 1e6 >= a.llera_km2
    sort, dist = sortides_xarxa(fd, np.argsort(-acc, axis=None, kind="stable"), urba, llera, ~np.isfinite(z), a.fora_max / dx, nf, nc)
    del acc, gr, dem

    # subconques de clavegueram: cel·les urbanes per sortida
    cel = np.flatnonzero(urba); s = sort.ravel()[cel]
    claus, inv = np.unique(s, return_inverse=True)
    area = np.bincount(inv) * dx * dx
    lmax = np.zeros(len(claus)); np.maximum.at(lmax, inv, dist.ravel()[cel] * dx)
    tc_c = np.full(256, C_ALTRES); [tc_c.__setitem__(k, v) for k, v in C_RACIONAL.items()]
    crac = np.bincount(inv, weights=tc_c[cob.ravel()[cel]]) / np.bincount(inv)
    a_ll = Transformer.from_crs(crs, 4326, always_xy=True).transform
    taula, cap10, cap2 = [], np.zeros(len(claus)), np.zeros(len(claus))
    for j, k in enumerate(claus):
        if k >= 0:
            r, c = divmod(int(k), nc)
        else:                                                          # surt de la reixa: la cel·la urbana més baixa de la subconca
            m = cel[inv == j]; r, c = divmod(int(m[np.nanargmin(z.ravel()[m])]), nc)
        x, y = tr.c + (c + 0.5) * dx, tr.f - (r + 0.5) * dx
        lon, lat = a_ll(x, y)
        tc_h = (a.t_entrada * 60 + lmax[j] / a.v_xarxa) / 3600
        cap10[j] = crac[j] * intensitat(idf_est, lon, lat, 10, tc_h); cap2[j] = crac[j] * intensitat(idf_est, lon, lat, 2, tc_h)
        tip = "fora" if k < 0 else "mar" if not np.isfinite(z[r, c]) else "llera" if llera[r, c] else "terreny"
        taula.append({"n": j + 1, "x": round(x, 1), "y": round(y, 1), "tipus": tip, "ha": round(area[j] / 1e4, 2),
                      "recorregut_m": round(float(lmax[j])), "c": round(float(crac[j]), 2), "tc_min": round(tc_h * 60, 1),
                      "cap_t10_mm_h": round(cap10[j], 1), "cap_t2_mm_h": round(cap2[j], 1)})

    # OSM: embornals i col·lectors cartografiats dins de cada subconca (només per comprovar)
    n_osm = {"embornals": 0, "pous": 0, "colectors": 0}
    if cami_osm and Path(cami_osm).exists():
        a_utm = Transformer.from_crs(4326, crs, always_xy=True).transform
        o = json.loads(gzip.open(cami_osm, "rt").read())
        num = np.zeros((nf, nc), np.int64); num.ravel()[cel] = inv + 1
        for lon, lat, t in o["punts"]:
            x, y = a_utm(lon, lat); c, r = int((x - tr.c) // dx), int((tr.f - y) // dx)
            if 0 <= r < nf and 0 <= c < nc and num[r, c]:
                fila = taula[num[r, c] - 1]; clau = "embornals" if t == "embornal" else "pous"
                fila[clau] = fila.get(clau, 0) + 1; n_osm[clau] += 1
        for lin in o["colectors"]:
            for lon, lat in lin[:: max(len(lin) // 5, 1)]:
                x, y = a_utm(lon, lat); c, r = int((x - tr.c) // dx), int((tr.f - y) // dx)
                if 0 <= r < nf and 0 <= c < nc and num[r, c]:
                    n_osm["colectors"] += 1; break

    b1 = np.zeros((nf, nc), np.uint32); b2 = np.zeros((nf, nc), np.uint32); b3 = np.zeros((nf, nc), np.uint32)
    b1.ravel()[cel] = np.rint(cap10[inv] * 10); b2.ravel()[cel] = np.rint(cap2[inv] * 10); b3.ravel()[cel] = inv + 1
    sup = area.sum()
    res = {"urba_km2": round(sup / 1e6, 3), "subconques": len(claus), "cap_t10_mm_h_mitjana": round(float((cap10 * area).sum() / sup), 1),
           "cap_t2_mm_h_mitjana": round(float((cap2 * area).sum() / sup), 1),
           "cap_t10_mm_h_rang": [round(float(cap10.min()), 1), round(float(cap10.max()), 1)],
           "sortides": {t: int(sum(1 for f in taula if f["tipus"] == t)) for t in ("llera", "terreny", "mar", "fora")}, "osm": n_osm}
    print(f"   {len(claus)} subconques; capacitat T10 {res['cap_t10_mm_h_mitjana']} mm/h (rang {res['cap_t10_mm_h_rang']}), "
          f"T2 {res['cap_t2_mm_h_mitjana']} mm/h; OSM {n_osm}")
    return desa(cid, np.stack([b1, b2, b3]), tr, crs, taula, res, a)


def desa(cid, bandes, tr, crs, taula, res, a=None):
    import rasterio
    nom = f"c{cid:03d}"
    if bandes is not None:
        TMP.mkdir(parents=True, exist_ok=True)
        with rasterio.open(TMP / f"clav10_{nom}.tif", "w", driver="GTiff", height=bandes.shape[1], width=bandes.shape[2], count=3,
                           dtype="uint32", crs=crs, transform=tr, nodata=0, compress="deflate", tiled=True) as f:
            f.write(bandes)
            f.update_tags(1, nom="capacitat T10 (dècimes de mm/h)"); f.update_tags(2, nom="capacitat T2 (dècimes de mm/h)")
            f.update_tags(3, nom=f"número de sortida (data/catalunya/clavegueram/{nom}.json)")
    SORT.mkdir(parents=True, exist_ok=True)
    par = None if a is None else {"pas_m": a.pas, "llera_km2": a.llera_km2, "t_entrada_min": a.t_entrada, "fora_max_m": a.fora_max, "v_xarxa_m_s": a.v_xarxa}
    (SORT / f"{nom}.json").write_text(json.dumps({"conca": cid, "crs": str(crs), "resum": res, "parametres": par, "sortides": taula},
                                                 ensure_ascii=False, separators=(",", ":")))
    return 0


def index():
    idx = {}
    for f in sorted(SORT.glob("c[0-9][0-9][0-9].json")):
        j = json.loads(f.read_text()); idx[f.stem] = j["resum"]
    (SORT / "index.json").write_text(json.dumps({"conques": idx, "release": "clavegueram-cat"}, ensure_ascii=False, indent=0))
    print(f"Clavegueram: {len(idx)} conques")


# ---------- per a simula_2d.py ----------
def a_reixa(tr, forma, crs, urba):
    """Capacitats T10 i T2 (mm/h) i índex lineal de la cel·la de sortida a una reixa qualsevol, a partir dels fitxers de la
    Release clavegueram-cat de les conques que hi toquen. Sortida -1 = l'aigua del clavegueram surt del domini.
    Les cel·les urbanes de la reixa sense dada queden a NaN (simula_2d.py hi posa la capacitat fixa). None si no n'hi ha cap."""
    import rasterio
    from rasterio.warp import reproject, Resampling, transform_bounds
    from shapely.geometry import shape, box
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import cobertes as COB
    nf, nc = forma; dx = tr.a
    b = box(*transform_bounds(crs, "EPSG:4326", tr.c, tr.f + nf * tr.e, tr.c + nc * dx, tr.f))
    c10 = np.full(forma, np.nan); c2 = np.full(forma, np.nan); so = np.full(forma, -1, np.int64); n = 0
    for i, f in COB.conques().items():
        if not shape(f["geometry"]).intersects(b):
            continue
        nom = f"c{i:03d}"; taula = SORT / f"{nom}.json"
        if not taula.exists():
            continue
        t = json.loads(taula.read_text())["sortides"]
        if not t:
            continue
        try:
            loc = TMP / f"clav10_{nom}.tif"                             # calculat en aquesta mateixa màquina
            with rasterio.open(loc if loc.exists() else URL.format(nom)) as src:
                bs = []
                for bd, rs in ((1, Resampling.max), (2, Resampling.max), (3, Resampling.mode)):
                    tmp = np.zeros(forma, np.uint32)
                    reproject(rasterio.band(src, bd), tmp, dst_transform=tr, dst_crs=crs, src_nodata=0, dst_nodata=0, resampling=rs)
                    bs.append(tmp)
        except Exception as e:
            print(f"Clavegueram {nom}: {type(e).__name__} {e}"); continue
        nou = (bs[0] > 0) & np.isnan(c10)
        c10[nou] = bs[0][nou] / 10; c2[nou] = bs[1][nou] / 10
        # sortida de cada subconca a la reixa de la simulació (des de la CRS de la conca)
        from pyproj import Transformer
        a_sim = Transformer.from_crs(json.loads(taula.read_text())["crs"], crs, always_xy=True).transform
        lin = np.full(len(t) + 1, -1, np.int64)
        for fila in t:
            if fila["tipus"] == "fora":
                continue
            x, y = a_sim(fila["x"], fila["y"]); c, r = int((x - tr.c) // dx), int((tr.f - y) // dx)
            if 0 <= r < nf and 0 <= c < nc:
                lin[fila["n"]] = r * nc + c
        so[nou] = lin[np.clip(bs[2][nou].astype(np.int64), 0, len(t))]
        n += 1
    if not n:
        return None
    print(f"Clavegueram estimat ({n} conques): {np.isfinite(c10[urba]).mean() * 100:.0f} % de les cel·les urbanes amb dada; "
          f"T10 {np.nanmean(c10[urba]):.0f} mm/h, T2 {np.nanmean(c2[urba]):.0f} mm/h; desguassa dins el domini el "
          f"{(so[urba & np.isfinite(c10)] >= 0).mean() * 100:.0f} %")
    return c10, c2, so


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="ordre", required=True)
    sub.add_parser("idf")
    o = sub.add_parser("osm"); o.add_argument("--pbf", required=True); o.add_argument("--sortida", required=True)
    c = sub.add_parser("conca"); c.add_argument("--id", type=int, required=True)
    c.add_argument("--mdt", required=True); c.add_argument("--cobertes", required=True); c.add_argument("--osm")
    c.add_argument("--pas", type=float, default=10, help="mida de cel·la (m)")
    c.add_argument("--llera-km2", type=float, default=0.5, help="conca a partir de la qual una cel·la és llera on desguassa la xarxa")
    c.add_argument("--t-entrada", type=float, default=10, help="temps d'entrada a la xarxa (min)")
    c.add_argument("--fora-max", type=float, default=200, help="tram no urbà màxim que pot travessar un col·lector fins a la llera (m)")
    c.add_argument("--v-xarxa", type=float, default=1.5, help="velocitat de l'aigua als col·lectors (m/s)")
    sub.add_parser("index")
    a = ap.parse_args()
    if a.ordre == "idf":
        return idf()
    if a.ordre == "osm":
        return osm(a.pbf, a.sortida)
    if a.ordre == "index":
        return index()
    return conca(a.id, a.mdt, a.cobertes, a.osm, a)


if __name__ == "__main__":
    sys.exit(main())
