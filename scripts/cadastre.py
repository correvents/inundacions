#!/usr/bin/env python3
"""Edificis del Cadastre (INSPIRE Buildings, Dirección General del Catastro) dins la caixa de la conca.

De cada part d'edifici (BuildingPart) en treu la planta, les plantes sobre rasant i les plantes sota rasant:
  - Les parts amb alguna planta sobre rasant són obstacles per a l'aigua a la simulació 2D (simula_2d.py --edificis).
    L'alçada aproximada és plantes × 3 m.
  - Les parts amb plantes sota rasant (soterranis, aparcaments) són les que es poden inundar si l'aigua arriba a tocar.
Els municipis es trien pel seu requadre al feed ATOM del Cadastre; si el feed no el porta, per nom (config "municipis_cadastre").

Sortides (data/<conca>/cadastre/):
  edificis.json.gz   parts en UTM 31N (EPSG:25831): [ref, plantes sobre, plantes sota, [anells [[x, y], ...]]]
  soterranis.geojson un element per edifici amb soterrani (lat/lon): ref, plantes sota/sobre, ús, habitatges
  edificis.png       edificis pintats per alçada, per al visor (amb resum.json: límits i recomptes)

Ús:  python3 scripts/cadastre.py --config config/ridaura.json"""
import argparse, gzip, io, json, re, sys, unicodedata, zipfile
from pathlib import Path
import numpy as np
import requests
from lxml import etree
from pyproj import Transformer
from shapely.geometry import Polygon, MultiPolygon, box, mapping
from shapely.ops import unary_union
from PIL import Image
import rasterio.features
from rasterio.transform import from_origin

ATOM = ["https://www.catastro.hacienda.gob.es/INSPIRE/buildings/ES.SDGC.BU.atom.xml",
        "http://www.catastro.minhap.es/INSPIRE/buildings/ES.SDGC.BU.atom.xml"]
PROV = {"ridaura": [17], "granollers": [8], "cabrera": [8]}
PAS_PNG = 4                                      # m per píxel de edificis.png
H_PLANTA = 3.0
CAIXES = {"ridaura": [2.86, 41.75, 3.09, 41.89]}  # la de relleu.py quan la config no en porta


def ln(el):
    return etree.QName(el).localname if isinstance(el.tag, str) else ""


def normal(s):
    s = unicodedata.normalize("NFD", s).encode("ascii", "ignore").decode().upper()
    return re.sub(r"[^A-Z]", "", s)


def get(url, **kw):
    for k in range(4):
        try:
            r = requests.get(requests.utils.requote_uri(url), timeout=180, **kw)
            r.raise_for_status()
            return r
        except Exception as e:
            print(f"  {url}: {e}"); err = e
    raise err


def xml_de(url):
    r = get(url)
    try:
        return etree.fromstring(r.content)
    except etree.XMLSyntaxError:
        print(f"{url} no és XML ({r.status_code}, {r.headers.get('content-type')}):\n{r.text[:1500]}"); raise


def feed_provincial(prov):
    """El feed general enllaça un feed per gerència (província)."""
    for top in ATOM:
        try:
            xml = xml_de(top)
        except Exception:
            continue
        hrefs = [x.get("href") for x in xml.iter() if ln(x) == "link" and x.get("href")]
        bo = [h for h in hrefs if re.search(rf"atom_{prov:02d}\.xml$", h, re.I)]
        if bo:
            return bo[0]
        print(f"{top}: cap enllaç a la gerència {prov:02d}. Enllaços: {hrefs[:20]}")
    raise RuntimeError("No trobo el feed ATOM del Cadastre")


def municipis(prov, caixa, noms):
    """(nom, url del zip) dels municipis del feed provincial que toquen la caixa."""
    xml = xml_de(feed_provincial(prov))
    cx = box(*caixa); tria = []
    entrades = [e for e in xml.iter() if ln(e) == "entry"]
    print(f"Feed {prov:02d}: {len(entrades)} entrades")
    for e in entrades:
        titol = next((x.text for x in e if ln(x) == "title"), "") or ""
        enllac = [x.get("href") for x in e if ln(x) == "link" and (x.get("href") or "").lower().endswith((".zip", ".xml"))]
        pol = next((x.text for x in e.iter() if ln(x) in ("polygon", "box")), None)
        toca = None
        if pol:
            v = [float(t) for t in pol.split()]
            if len(v) == 4:                        # georss:box lat_s lon_w lat_n lon_e
                g = box(v[1], v[0], v[3], v[2])
            else:                                  # georss:polygon lat lon lat lon ...
                g = Polygon(list(zip(v[1::2], v[0::2])))
            toca = g.intersects(cx)
        if toca is None:
            toca = any(normal(n) in normal(titol) for n in noms)
        if toca and enllac:
            tria.append((titol, enllac[0]))
    if not tria:
        print("Cap municipi trobat. Mostra del feed:\n" + etree.tostring(entrades[1] if len(entrades) > 1 else xml, pretty_print=True).decode()[:3000])
    return tria


def zip_de(url):
    """El feed provincial pot apuntar al zip o a un feed municipal que hi apunta."""
    if url.lower().endswith(".zip"):
        return url
    xml = xml_de(url)
    return next(x.get("href") for x in xml.iter() if ln(x) == "link" and (x.get("href") or "").lower().endswith(".zip"))


def anells(el, a_utm):
    """Anells (exterior primer) dels polígons d'un element GML, en UTM 31N."""
    out = []
    for pl in el.iter():
        if ln(pl) not in ("PolygonPatch", "Polygon"):
            continue
        for anell in pl.iter():
            if ln(anell) != "posList":
                continue
            v = np.array(anell.text.split(), dtype=float)
            dim = int(anell.get("srsDimension") or 2)
            xy = v.reshape(-1, dim)[:, :2]
            out.append(a_utm(xy))
    return out


def transformador(srs):
    epsg = re.findall(r"(\d{4,5})\s*$", srs or "")
    epsg = int(epsg[0]) if epsg else 25831
    if epsg == 25831:
        return lambda xy: xy
    t = Transformer.from_crs(epsg, 25831, always_xy=epsg not in (4258, 4326))
    return lambda xy: np.column_stack(t.transform(xy[:, 0], xy[:, 1]))


def llegeix(zbytes, caixa_utm):
    z = zipfile.ZipFile(io.BytesIO(zbytes))
    noms = z.namelist()
    fpart = next(n for n in noms if n.lower().endswith("buildingpart.gml"))
    fedif = next(n for n in noms if n.lower().endswith("building.gml"))
    srs = re.search(rb'srsName="([^"]+)"', z.read(fpart)[:20000])
    a_utm = transformador(srs.group(1).decode() if srs else "")
    parts = []
    for _, el in etree.iterparse(z.open(fpart), tag="{*}BuildingPart"):
        gid = el.get("{http://www.opengis.net/gml/3.2}id") or ""
        ref = gid.replace("ES.SDGC.BU.", "").split("_part")[0]
        v = {ln(x): x.text for x in el.iter() if ln(x) in ("numberOfFloorsAboveGround", "numberOfFloorsBelowGround")}
        sobre, sota = int(v.get("numberOfFloorsAboveGround") or 0), int(v.get("numberOfFloorsBelowGround") or 0)
        an = anells(el, a_utm)
        el.clear()
        if not an:
            continue
        x0, y0 = an[0].min(axis=0); x1, y1 = an[0].max(axis=0)
        if x1 < caixa_utm[0] or x0 > caixa_utm[2] or y1 < caixa_utm[1] or y0 > caixa_utm[3]:
            continue
        parts.append([ref, sobre, sota, [np.round(a, 1).tolist() for a in an]])
    edif = {}
    for _, el in etree.iterparse(z.open(fedif), tag="{*}Building"):
        gid = el.get("{http://www.opengis.net/gml/3.2}id") or ""
        ref = gid.replace("ES.SDGC.BU.", "")
        v = {ln(x): (x.text or "").strip() for x in el.iter() if ln(x) in ("currentUse", "numberOfDwellings", "numberOfBuildingUnits")}
        edif[ref] = {"us": v.get("currentUse", ""), "habitatges": int(v.get("numberOfDwellings") or 0)}
        el.clear()
    return parts, edif


US = {"1_residential": "residencial", "2_agriculture": "agrari", "3_industrial": "industrial", "4_1_office": "oficines",
      "4_2_retail": "comercial", "4_3_publicServices": "serveis públics"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    conca = cfg["conca"]
    caixa = cfg.get("caixa") or CAIXES[conca]
    a_utm = Transformer.from_crs(4326, 25831, always_xy=True)
    xs, ys = a_utm.transform([caixa[0], caixa[2], caixa[0], caixa[2]], [caixa[1], caixa[1], caixa[3], caixa[3]])
    caixa_utm = [min(xs), min(ys), max(xs), max(ys)]
    noms = cfg.get("municipis_cadastre", [])

    parts, edif = [], {}
    for prov in PROV.get(conca, cfg.get("provincies_cadastre", [])):
        for titol, url in municipis(prov, caixa, noms):
            try:
                zb = get(zip_de(url)).content
                p, e = llegeix(zb, caixa_utm)
            except Exception as ex:
                print(f"{titol}: {type(ex).__name__} {ex}"); continue
            print(f"{titol}: {len(p)} parts dins la caixa, {sum(1 for x in p if x[2] > 0)} amb plantes sota rasant")
            parts += p; edif.update(e)
    if not parts:
        print("Cap edifici."); return 1

    d = Path("data") / conca / "cadastre"; d.mkdir(parents=True, exist_ok=True)
    with gzip.open(d / "edificis.json.gz", "wt") as f:
        json.dump({"crs": "EPSG:25831", "h_planta_m": H_PLANTA, "camps": ["ref", "plantes_sobre", "plantes_sota", "anells"], "parts": parts},
                  f, separators=(",", ":"))

    # soterranis: un element per edifici, amb les parts que tenen plantes sota rasant
    a_ll = Transformer.from_crs(25831, 4326, always_xy=True)
    per_ref = {}
    for ref, sobre, sota, an in parts:
        r = per_ref.setdefault(ref, {"sota": 0, "sobre": 0, "pol": []})
        r["sobre"] = max(r["sobre"], sobre)
        if sota > 0:
            r["sota"] = max(r["sota"], sota)
            ll = [np.column_stack(a_ll.transform(np.array(x)[:, 0], np.array(x)[:, 1])) for x in an]
            r["pol"].append(Polygon(ll[0], ll[1:]).buffer(0))
    fts = []
    for ref, r in per_ref.items():
        if not r["sota"]:
            continue
        g = unary_union(r["pol"])
        if g.geom_type == "GeometryCollection":
            g = unary_union([x for x in g.geoms if x.geom_type in ("Polygon", "MultiPolygon")])
        g = g.simplify(0.000005)
        if g.is_empty or g.geom_type not in ("Polygon", "MultiPolygon"):
            continue
        e = edif.get(ref, {})
        fts.append({"type": "Feature", "geometry": json.loads(json.dumps(mapping(g), default=float)),
                    "properties": {"ref": ref, "plantes_sota": r["sota"], "plantes_sobre": r["sobre"],
                                   "us": US.get(e.get("us"), e.get("us", "")), "habitatges": e.get("habitatges", 0)}})
    def arrodoneix(c):
        return [arrodoneix(x) for x in c] if isinstance(c[0], (list, tuple)) else [round(c[0], 6), round(c[1], 6)]
    for f in fts:
        f["geometry"]["coordinates"] = arrodoneix(f["geometry"]["coordinates"])
    (d / "soterranis.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": fts}, ensure_ascii=False, separators=(",", ":")))

    # imatge dels edificis per alçada (lat/lon), per al visor
    pas_g = PAS_PNG / 111320
    w, s, e_, n = caixa
    nc_, nf_ = int((e_ - w) / pas_g * np.cos(np.radians((s + n) / 2))), int((n - s) / pas_g)
    tr = from_origin(w, n, (e_ - w) / nc_, (n - s) / nf_)
    formes = []
    for ref, sobre, sota, an in parts:
        if sobre < 1:
            continue
        ll = [np.column_stack(a_ll.transform(np.array(x)[:, 0], np.array(x)[:, 1])) for x in an]
        formes.append((Polygon(ll[0], ll[1:]), min(sobre, 250)))
    formes.sort(key=lambda f: f[1])                                     # les més altes per sobre
    pl = rasterio.features.rasterize(formes, out_shape=(nf_, nc_), transform=tr, fill=0, dtype="uint8")
    rgba = np.zeros((nf_, nc_, 4), dtype=np.uint8)
    for llindar, c in ((1, (189, 189, 189, 200)), (3, (150, 150, 150, 215)), (6, (99, 99, 99, 230)), (10, (37, 37, 37, 240))):
        rgba[pl >= llindar] = c
    Image.fromarray(rgba, "RGBA").save(d / "edificis.png", optimize=True)
    sot = [f["properties"]["plantes_sota"] for f in fts]
    (d / "resum.json").write_text(json.dumps({
        "font": "Dirección General del Catastro, INSPIRE Buildings", "limits": [[s, w], [n, e_]], "png": "edificis.png",
        "llegenda_plantes": [1, 3, 6, 10], "h_planta_m": H_PLANTA, "parts": len(parts), "parts_obstacle": len(formes),
        "edificis_amb_soterrani": len(fts), "soterranis_per_plantes": {str(k): sot.count(k) for k in sorted(set(sot))}}, ensure_ascii=False, indent=1))
    print(f"{len(parts)} parts ({len(formes)} obstacles), {len(fts)} edificis amb soterrani → {d}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
