#!/usr/bin/env python3
"""Passos de torrentets sota carreteres i vies de tren, conca per conca a tot Catalunya, amb escenaris obert / mig tapat / tapat.

Idea: l'MDT de 5 m de l'IGN (LiDAR) és el terreny nu. Una carretera o una via de tren en terraplè hi queda com un dic
massís: l'obra de drenatge (tub, marc, pas inferior) no s'hi veu. Per tant, aigües amunt de cada terraplè que talla un
torrentet, l'MDT té una depressió tancada, i aquesta depressió és exactament el llac que es faria si l'obra es tapés.
Sense cap número triat a mà per a cada pas:

1. Vies d'OpenStreetMap (Geofabrik, extracte de Catalunya): carreteres (motorway → unclassified, residential) i vies de
   tren. Els trams en pont o en túnel no fan de dic. Les rieres amb tunnel=culvert d'OSM són l'inventari conegut
   d'obres de drenatge i serveixen per comprovar la detecció.
2. MDT 5 m de la conca (Release mdt5-cat) → omplir depressions, direccions de drenatge i acumulació (pysheds).
3. Cada depressió (omplert − terreny) que toca una via i té conca ≥ --conca-min km², fondària ≥ --fons-min m i volum ≥
   --volum-min m³ és un pas probable. El pas es posa al punt de la via més proper al fons de la depressió. Es desa:
   - conca drenada (km²; només la de dins la conca del tall de 100 km²),
   - llac si es tapa: superfície, volum, fondària màxima i contorn,
   - punt de vessament (per on passa l'aigua quan el llac és ple) i si vessa per sobre de la mateixa via,
   - si coincideix amb una obra de drenatge dibuixada a OSM.
4. Torrentets que travessen una via sense depressió ni pont a OSM: pas a nivell (gual). L'aigua hi passa per sobre amb
   qualsevol cabal.
5. Escenaris de cada pas:
   - tapat: l'obra no desguassa. L'aigua omple el llac i vessa per sobre quan l'escolament de la conca supera
     volum/conca (mm).
   - mig tapat: la meitat de la secció. Desguassa la meitat del cabal de disseny.
   - obert: desguassa el cabal de disseny de la norma de carreteres 5.2-IC (crescuda de 100 anys de la conca de l'obra),
     calculat com a scripts/passos.py: Gumbel regional dels màxims de pluja diària de la XEMA, Témez, Kt, Ka, i P0 =
     absorció inicial de l'episodi amb el terreny més moll del model (la norma dimensiona amb terreny humit). Només on hi ha XEMA i episodis (conques pilot
     i les veïnes que cauen dins la caixa de simulació d'una pilot, amb la pluja d'aquesta); a la resta, null.
   Qualsevol cabal entrant per sobre del de l'escenari va omplint el llac; quan n'entra el volum, l'aigua passa per
   sobre de la via.

Sortides:  data/catalunya/passos/cNNN.geojson (punts dels passos, contorns dels llacs, guals)
           data/catalunya/passos/index.json (resum per conca)

Ús:  python3 scripts/passos_cat.py osm --pbf cataluna-latest.osm.pbf --sortida /tmp/vies.json.gz
     python3 scripts/passos_cat.py conca --id 188 --mdt mdt5_c188.tif --vies /tmp/vies.json.gz
     python3 scripts/passos_cat.py index"""
import argparse, gzip, json, sys, tempfile
from pathlib import Path
import numpy as np
np.in1d = np.isin                      # pysheds encara fa servir np.in1d (numpy 2)
import rasterio
from rasterio import features
from pyproj import Transformer
from scipy import ndimage
from shapely.geometry import shape, mapping, LineString, Point
from shapely.ops import transform as stransform
from shapely.strtree import STRtree

CONQUES = Path("data/catalunya/conques/conques_100km2.geojson")
SORT = Path("data/catalunya/passos")
VIALS = ("motorway", "trunk", "primary", "secondary", "tertiary", "unclassified", "residential", "living_street",
         "motorway_link", "trunk_link", "primary_link", "secondary_link", "tertiary_link")
TRENS = ("rail", "light_rail", "narrow_gauge")
RIERES = ("river", "stream", "ditch", "drain", "canal", "brook")
# conques del tall de 100 km² amb XEMA i episodis del model (per al cabal de disseny de l'escenari obert)
PILOTS = {124: "ridaura", 125: "ridaura", 188: "granollers", 199: "cabrera"}
CAIXA_RIDAURA = [2.86, 41.75, 3.09, 41.89]       # la de relleu.py (config/ridaura.json no en porta)


def pilot_de(cid, geom_ll=None):
    """Conca pilot que dona la pluja de disseny: la mateixa conca, o la pilot el domini de simulació de la qual (caixa de
    config/<pilot>.json) toca la conca. Així els passos de les conques veïnes que cauen dins la simulació 2D d'una
    pilot també tenen cabal de disseny, amb la pluja regional de la pilot."""
    if cid in PILOTS:
        return PILOTS[cid]
    if geom_ll is None:
        return None
    from shapely.geometry import box
    for nom in dict.fromkeys(PILOTS.values()):
        cfg = Path("config") / f"{nom}.json"
        caixa = json.loads(cfg.read_text()).get("caixa", CAIXA_RIDAURA if nom == "ridaura" else None) if cfg.exists() else None
        if caixa and geom_ll.intersects(box(*caixa)):
            return nom
    return None


# ---------- 1. OSM ----------
def osm(pbf, sortida):
    import osmium

    class H(osmium.SimpleHandler):
        def __init__(self):
            super().__init__(); self.vies = []; self.obres = []

        def way(self, w):
            t = w.tags
            hw, rw, ww = t.get("highway"), t.get("railway"), t.get("waterway")
            via = hw in VIALS or rw in TRENS
            obra = ww in RIERES and (t.get("tunnel") in ("culvert", "yes") or t.get("culvert") == "yes")
            if not (via or obra):
                return
            try:
                xy = [v for n in w.nodes for v in (round(n.lon, 6), round(n.lat, 6))]
            except osmium.InvalidLocationError:
                return
            if len(xy) < 4:
                return
            if via:
                pont = t.get("bridge", "no") != "no"
                tunel = t.get("tunnel", "no") != "no" or t.get("layer", "0").startswith("-")
                self.vies.append([w.id, hw or ("tren" if rw != "light_rail" else "tramvia"), t.get("ref") or t.get("name") or "",
                                  "p" if pont else ("t" if tunel else ""), xy])
            else:
                self.obres.append([w.id, ww, t.get("name") or "", xy])

    h = H()
    h.apply_file(pbf, locations=True, idx="flex_mem")
    with gzip.open(sortida, "wt") as f:
        json.dump({"vies": h.vies, "obres": h.obres}, f, separators=(",", ":"))
    print(f"OSM: {len(h.vies)} trams de via ({sum(v[3] == 'p' for v in h.vies)} ponts), {len(h.obres)} obres de drenatge")


# ---------- 5. cabal de disseny (5.2-IC), com a scripts/passos.py ----------
def p24_t100(conca):
    """Pluja de 100 anys en 24 h: Gumbel regional (índex de crescuda) amb els màxims anuals de pluja diària de la XEMA."""
    import pandas as pd
    norm, mitj = [], []
    for f in sorted((Path("data") / conca / "xema").glob("*.csv")):
        if f.name.startswith("episodis"):
            continue
        s = pd.read_csv(f, parse_dates=["temps"]).set_index("temps").mm
        dia = s.resample("1D").agg(["sum", "count"])
        dia = dia[dia["count"] >= 44]                                  # dies amb ≥ 22 h de dades
        anual = dia.groupby(dia.index.year)["sum"].agg(["max", "count"])
        anual = anual[anual["count"] >= 330]["max"]                    # anys gairebé complets
        if len(anual) >= 5:
            mitj.append(anual.mean()); norm += list(anual / anual.mean())
    if len(norm) < 20:
        return None, len(norm)
    x = np.array(norm); beta = x.std(ddof=1) * np.sqrt(6) / np.pi; mu = x.mean() - 0.5772 * beta
    return 1.13 * float(np.mean(mitj)) * (mu - beta * np.log(-np.log(1 - 1 / 100))), len(norm)


def q_disseny(area_km2, l_km, j, p24, p0):
    tc = max(0.3 * (l_km / max(j, 1e-3) ** 0.25) ** 0.76, 0.25)       # Témez, h
    ka = 1 - np.log10(area_km2) / 15 if area_km2 > 1 else 1.0
    pd_ = p24 * ka; id_ = pd_ / 24
    i = id_ * 11 ** ((28 ** 0.1 - tc ** 0.1) / (28 ** 0.1 - 1))
    kt = 1 + tc ** 1.25 / (tc ** 1.25 + 14)
    r = pd_ / p0
    c = (r - 1) * (r + 23) / (r + 11) ** 2 if r > 1 else 0.0
    return c * i * area_km2 * kt / 3.6, tc


def pluja_pilot(cid, geom_ll=None):
    nom = pilot_de(cid, geom_ll)
    if not nom:
        return None
    p24, n = p24_t100(nom)
    idx = Path("data") / nom / "episodis" / "index.json"
    ias = [e["ia_mm"] for e in json.loads(idx.read_text())["episodis"] if e.get("ia_mm")] if idx.exists() else []
    if not p24 or not ias:
        print(f"{nom}: sense cabal de disseny ({n} anys-estació XEMA, {len(ias)} episodis)")
        return None
    p0 = float(min(ias))                                               # terreny moll: l'episodi amb menys absorció
    print(f"{nom}: pluja de 100 anys en 24 h {p24:.0f} mm ({n} anys-estació XEMA); P0 = {p0:.0f} mm")
    return {"p24_t100_mm": round(p24), "p0_mm": round(p0), "anys_estacio": n, "_p24": p24, "_p0": p0}


# ---------- 2-4. una conca ----------
def conca(cid, cami_mdt, vies_gz, a):
    C = {f["properties"]["id"]: f for f in json.loads(CONQUES.read_text())["features"]}
    geom = shape(C[cid]["geometry"])
    with rasterio.open(cami_mdt) as f:
        z = f.read(1).astype("float32"); nd = f.nodata; tr = f.transform; crs = f.crs
    z[z == nd] = np.nan; z /= 100
    z[z <= 0.2] = np.nan                                                # mar
    nf, nc = z.shape; cel = abs(tr.a * tr.e)
    a_utm = Transformer.from_crs(4326, crs, always_xy=True).transform
    a_ll = Transformer.from_crs(crs, 4326, always_xy=True).transform
    w, s, e, n = geom.bounds                                           # lon/lat
    print(f"== Conca {cid}: {C[cid]['properties']['km2']} km², reixa {nc}×{nf}")

    # vies de la caixa de la conca
    with gzip.open(vies_gz, "rt") as f:
        osmj = json.load(f)
    def dins(xy):
        xs, ys = xy[0::2], xy[1::2]
        return max(xs) >= w and min(xs) <= e and max(ys) >= s and min(ys) <= n
    vies = []
    for vid, tipus, nom, fl, xy in osmj["vies"]:
        if dins(xy):
            g = LineString([a_utm(xy[k], xy[k + 1]) for k in range(0, len(xy), 2)])
            vies.append({"id": vid, "tipus": tipus, "nom": nom, "fl": fl, "g": g})
    obres = [LineString([a_utm(xy[k], xy[k + 1]) for k in range(0, len(xy), 2)]) for _, _, _, xy in osmj["obres"] if dins(xy)]
    massisses = [v for v in vies if not v["fl"]]
    ponts = [v for v in vies if v["fl"] == "p"]
    print(f"OSM: {len(vies)} trams de via ({len(ponts)} ponts), {len(obres)} obres de drenatge dibuixades")
    if not massisses:
        return desa(cid, [], {"km2": C[cid]["properties"]["km2"], "vies": 0})

    # rasters de vies: índex (1-based) del tram que passa per cada cel·la
    via_r = features.rasterize(((mapping(v["g"]), k + 1) for k, v in enumerate(massisses)), out_shape=z.shape,
                               transform=tr, all_touched=True, dtype="int32")
    pont_r = features.rasterize(((mapping(v["g"].buffer(10)), 1) for v in ponts), out_shape=z.shape, transform=tr,
                                dtype="uint8") if ponts else np.zeros(z.shape, "uint8")

    # 2. pysheds
    from pysheds.grid import Grid
    tmp = Path(tempfile.mkdtemp()) / "z.tif"
    with rasterio.open(tmp, "w", driver="GTiff", height=nf, width=nc, count=1, dtype="float32", crs=crs, transform=tr,
                       nodata=np.nan) as f:
        f.write(z, 1)
    gr = Grid.from_raster(str(tmp)); dem = gr.read_raster(str(tmp))
    ple = gr.fill_depressions(gr.fill_pits(dem))
    fd = gr.flowdir(gr.resolve_flats(ple))
    acc = np.asarray(gr.accumulation(fd), dtype="float64")
    zp = np.asarray(ple, dtype="float32")
    fons = np.where(np.isfinite(z), zp - z, 0).astype("float32")
    del ple
    print("Drenatge fet")

    # 3. depressions que toquen una via
    lab, nlab = ndimage.label(fons > 0.05, structure=np.ones((3, 3)))
    idx = np.arange(1, nlab + 1)
    prof = ndimage.maximum(fons, lab, idx); vol = ndimage.sum(fons, lab, idx) * cel
    candidats = idx[(prof >= a.fons_min) & (vol >= a.volum_min)]
    sl = ndimage.find_objects(lab)
    via_dil = ndimage.binary_dilation(via_r > 0, iterations=a.peu)
    arbre_o = STRtree(obres) if obres else None
    pilot = pluja_pilot(cid, shape(C[cid]["geometry"]))
    rec = recorreguts(fd, acc, z) if pilot else None
    fts, passos, ocupat = [], [], np.zeros(z.shape, bool)
    for k in candidats:
        r0, r1 = max(sl[k - 1][0].start - a.peu - 1, 0), min(sl[k - 1][0].stop + a.peu + 1, nf)
        c0, c1 = max(sl[k - 1][1].start - a.peu - 1, 0), min(sl[k - 1][1].stop + a.peu + 1, nc)
        m = lab[r0:r1, c0:c1] == k
        anell = ndimage.binary_dilation(m, iterations=a.peu) & ~m
        tocs = via_r[r0:r1, c0:c1][anell]
        tocs = tocs[tocs > 0]
        if not len(tocs):
            continue
        vora = ndimage.binary_dilation(m) & ~m
        ac = np.where(vora, acc[r0:r1, c0:c1], -1)
        rs, cs = np.unravel_index(np.argmax(ac), ac.shape)
        area = float(ac[rs, cs]) * cel / 1e6
        if area < a.conca_min:
            continue
        # punt del pas: el tram que més toca, al punt més proper al fons de la depressió
        v = massisses[int(np.bincount(tocs).argmax()) - 1]
        f_ = np.where(m, fons[r0:r1, c0:c1], -1); rf, cf = np.unravel_index(np.argmax(f_), f_.shape)
        xf, yf = tr * (c0 + cf + 0.5, r0 + rf + 0.5)
        pp = v["g"].interpolate(v["g"].project(Point(xf, yf)))
        xs, ys = tr * (c0 + cs + 0.5, r0 + rs + 0.5)
        vessa_via = bool(via_dil[r0 + rs, c0 + cs])
        osm_obra = bool(arbre_o is not None and len(arbre_o.query(pp.buffer(a.osm_dist), predicate="intersects")))
        cp, rp = ~tr * (pp.x, pp.y)
        if 0 <= rp < nf and 0 <= cp < nc and pont_r[int(rp), int(cp)]:
            continue                                                    # pont dibuixat a OSM: no fa de dic
        ocupat[r0:r1, c0:c1] |= ndimage.binary_dilation(m, iterations=3)
        V, A, H = float(vol[k - 1]), float(m.sum() * cel), float(prof[k - 1])
        esc = {"tapat": {"q_m3s": 0, "vessa_amb_escolament_mm": round(V / (area * 1e6) * 1000, 1)}}
        if pilot:
            q, tc = cabal_obra(rec, z, tr, r0 + rs, c0 + cs, area, pilot)
            if q is not None:
                esc["mig_tapat"] = {"q_m3s": round(q / 2, 2)}
                esc["obert"] = {"q_m3s": round(q, 2), "tc_h": round(tc, 2)}
        pid = f"{cid}-{len(passos) + 1}"
        prop = {"id": pid, "tipus": "terraplè", "via": v["tipus"], "nom": v["nom"], "osm_via": v["id"],
                "conca_km2": round(area, 3), "llac_m3": round(V), "llac_m2": round(A), "prof_m": round(H, 1),
                "vessa_per_la_via": vessa_via, "vessament": [round(c, 5) for c in a_ll(xs, ys)],
                "obra_osm": osm_obra, "escenaris": esc}
        passos.append(prop)
        fts.append({"type": "Feature", "geometry": {"type": "Point", "coordinates": [round(c, 5) for c in a_ll(pp.x, pp.y)]},
                    "properties": prop})
        # contorn del llac (tapat)
        sub_tr = tr * tr.translation(c0, r0)
        for gj, val in features.shapes(m.astype("uint8"), mask=m, transform=sub_tr):
            pol = shape(gj).simplify(cel * 0.6)
            if pol.area >= 4 * cel:
                pol = stransform(lambda x, y: tuple(np.round(a_ll(x, y), 5)), pol)  # noqa
                fts.append({"type": "Feature", "geometry": mapping(pol), "properties": {"llac_de": pid}})

    # 4. guals: torrentet que travessa una via a nivell, sense depressió ni pont
    rierol = acc * cel / 1e6 >= a.conca_min
    gual = rierol & (via_r > 0) & ~ocupat & ~(pont_r > 0)
    lg, ng = ndimage.label(gual, structure=np.ones((3, 3)))
    if ng:
        pos = ndimage.maximum_position(acc, lg, np.arange(1, ng + 1))
        for (r, c) in pos:
            area = float(acc[r, c]) * cel / 1e6
            if area < a.gual_min:
                continue
            v = massisses[int(via_r[r, c]) - 1]
            x, y = tr * (c + 0.5, r + 0.5)
            pid = f"{cid}-g{sum(p['tipus'] == 'gual' for p in passos) + 1}"
            prop = {"id": pid, "tipus": "gual", "via": v["tipus"], "nom": v["nom"], "osm_via": v["id"], "conca_km2": round(area, 3),
                    "obra_osm": bool(arbre_o is not None and len(arbre_o.query(Point(x, y).buffer(a.osm_dist), predicate="intersects")))}
            passos.append(prop)
            fts.append({"type": "Feature", "geometry": {"type": "Point", "coordinates": [round(q, 5) for q in a_ll(x, y)]},
                        "properties": prop})

    # comprovació amb l'inventari d'OSM: obres dibuixades que creuen una via i que hem trobat
    creuen = [o for o in obres if any(o.intersects(massisses[j]["g"]) for j in STRtree([v["g"] for v in massisses]).query(o))]
    pts = [Point(f["geometry"]["coordinates"]) for f in fts if f["geometry"]["type"] == "Point" and f["properties"]["tipus"] == "terraplè"]
    pts = [Point(a_utm(p.x, p.y)) for p in pts]
    trobades = sum(1 for o in creuen if any(o.distance(p) <= a.osm_dist for p in pts)) if pts else 0
    t = [p for p in passos if p["tipus"] == "terraplè"]
    res = {"km2": C[cid]["properties"]["km2"], "vies": len(massisses), "passos": len(t),
           "guals": len(passos) - len(t), "vessen_per_la_via": sum(p["vessa_per_la_via"] for p in t),
           "obres_osm_sota_via": len(creuen), "obres_osm_trobades": trobades, "pilot": pilot and {k: v for k, v in pilot.items() if not k.startswith("_")}}
    print(json.dumps(res, ensure_ascii=False))
    for p in sorted(t, key=lambda p: -p["llac_m3"])[:12]:
        print(f"  {p['id']:>8} {p['via']:12s} {p['nom'][:28]:28s} conca {p['conca_km2']:7.2f} km²  llac {p['llac_m3']:8.0f} m³ "
              f"{p['prof_m']:4.1f} m  vessa {'per la via' if p['vessa_per_la_via'] else 'per un altre lloc'}  "
              f"tapat: {p['escenaris']['tapat']['vessa_amb_escolament_mm']} mm  obert: {p['escenaris'].get('obert', {}).get('q_m3s')}")
    return desa(cid, fts, res)


def recorreguts(fd, acc, z):
    """Per a cada cel·la: recorregut d'aigua més llarg que hi arriba (cel·les; diagonal √2) i cota màxima de la seva conca.
    Es recorren les cel·les de menys a més acumulació: una cel·la sempre en té menys que la de sota."""
    import numba

    @numba.njit(cache=True)
    def passa(ordre, fd, lup, zmax, nf, nc):
        dr = {64: -1, 128: -1, 1: 0, 2: 1, 4: 1, 8: 1, 16: 0, 32: -1}
        dc = {64: 0, 128: 1, 1: 1, 2: 1, 4: 0, 8: -1, 16: -1, 32: -1}
        for k in ordre:
            r = k // nc; c = k % nc; d = fd[r, c]
            if d not in dr:
                continue
            r2 = r + dr[d]; c2 = c + dc[d]
            if r2 < 0 or r2 >= nf or c2 < 0 or c2 >= nc:
                continue
            pas = 1.4142135 if dr[d] != 0 and dc[d] != 0 else 1.0
            if lup[r, c] + pas > lup[r2, c2]:
                lup[r2, c2] = lup[r, c] + pas
            if zmax[r, c] > zmax[r2, c2]:
                zmax[r2, c2] = zmax[r, c]

    nf, nc = z.shape
    lup = np.zeros(z.shape, "float32"); zmax = np.where(np.isfinite(z), z, -1e9).astype("float32")
    passa(np.argsort(acc, axis=None, kind="stable"), np.asarray(fd).astype("int64"), lup, zmax, nf, nc)
    return lup, zmax


def cabal_obra(rec, z, tr, r, c, area, pilot):
    """Cabal de disseny 5.2-IC de l'obra amb sortida a (r, c): recorregut més llarg i pendent de la seva conca."""
    lup, zmax = rec
    l_km = float(lup[r, c]) * tr.a / 1000
    j = (float(zmax[r, c]) - float(z[r, c])) / max(l_km * 1000, 1)
    return q_disseny(area, max(l_km, 0.05), j, pilot["_p24"], pilot["_p0"])


def desa(cid, fts, res):
    SORT.mkdir(parents=True, exist_ok=True)
    (SORT / f"c{cid:03d}.geojson").write_text(json.dumps({"type": "FeatureCollection", "resum": res, "features": fts},
                                                         ensure_ascii=False, separators=(",", ":")))
    return 0


def index():
    idx = {}
    for f in sorted(SORT.glob("c*.geojson")):
        idx[int(f.stem[1:])] = json.loads(f.read_text())["resum"]
    tot = {k: sum((r.get(k) or 0) for r in idx.values()) for k in ("passos", "guals", "vessen_per_la_via", "obres_osm_sota_via", "obres_osm_trobades")}
    (SORT / "index.json").write_text(json.dumps({"conques": idx, "total": tot}, ensure_ascii=False, separators=(",", ":")))
    print(f"{len(idx)} conques: {tot}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="ordre", required=True)
    o = sub.add_parser("osm"); o.add_argument("--pbf", required=True); o.add_argument("--sortida", required=True)
    c = sub.add_parser("conca"); c.add_argument("--id", type=int, required=True); c.add_argument("--mdt", required=True)
    c.add_argument("--vies", required=True)
    c.add_argument("--conca-min", type=float, default=0.05, help="conca mínima d'un torrentet (km², la de les rieretes)")
    c.add_argument("--gual-min", type=float, default=0.25, help="conca mínima d'un gual per desar-lo (km²)")
    c.add_argument("--fons-min", type=float, default=0.5, help="fondària mínima del llac si es tapa (m)")
    c.add_argument("--volum-min", type=float, default=100, help="volum mínim del llac si es tapa (m³)")
    c.add_argument("--peu", type=int, default=4, help="cel·les entre l'eix de la via i el peu del terraplè (4 = 20 m)")
    c.add_argument("--osm-dist", type=float, default=30, help="distància per casar amb una obra dibuixada a OSM (m)")
    sub.add_parser("index")
    a = ap.parse_args()
    if a.ordre == "osm":
        return osm(a.pbf, a.sortida)
    if a.ordre == "index":
        return index()
    return conca(a.id, a.mdt, a.vies, a)


if __name__ == "__main__":
    sys.exit(main())
