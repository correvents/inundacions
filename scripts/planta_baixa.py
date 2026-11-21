#!/usr/bin/env python3
"""Ús de la planta baixa (i del soterrani) de cada edifici del Cadastre, per saber QUÈ hi ha on arriba l'aigua.

Dues fonts, de més a menys exacta:
  1. Consulta per parcel·la del Cadastre (OVC Consulta_DNPRC, oberta): llista cada immoble de la parcel·la amb la seva planta
     i el seu ús (Residencial, Comercial, Almacén-Estacionamiento...). Dona l'ús real de la planta baixa. Una petició per
     parcel·la: es fa només per a les conques triades (per defecte les pilot) i es desa a la Release `planta-baixa-cat`
     (pb_ovc_cNNN.json.gz, resposta resumida [[planta, ús], ...] per parcel·la) per no repetir-la.
  2. INSPIRE Buildings (tot Catalunya): ús principal de l'edifici, habitatges i unitats. Si no hi ha la consulta per parcel·la,
     la planta baixa s'estima: unifamiliar → habitatge; plurifamiliar amb més unitats que habitatges → «habitatge o local».
     Es desa a la Release `planta-baixa-cat` com us_MMMMM.json.gz {ref: [ús, habitatges, unitats]}.

Sortida per conca (tall de 100 km²): data/catalunya/planta_baixa/cNNN.json
  {"classes": [...], "camps": ["font", "classe planta baixa", "habitatges a planta baixa", "immobles a planta baixa",
   "classe soterrani"], "edificis": {ref: [...]}}   font 0 = consulta per parcel·la, 1 = estimat amb INSPIRE.
La «classe planta baixa» és la més sensible dels immobles de la planta (habitatge > equipament > comerç > ...).
simula_2d.py --planta-baixa la fa servir (opcional) per comptar quines plantes baixes tenen aigua a tocar i de quin ús.

Ús:  python3 scripts/planta_baixa.py inspire --lot 3 --lots 16 [--salta noms_ja_pujats.txt]
     python3 scripts/planta_baixa.py ovc --conca 125 [--max 50] [--previ pb_ovc_c125.json.gz]
     python3 scripts/planta_baixa.py conca --conca 125"""
import argparse, gzip, io, json, re, sys, time, zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import requests
from lxml import etree
from shapely.geometry import shape, Point
from shapely.prepared import prep
from pyproj import Transformer

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cadastre as C

TMP = Path("/tmp/planta_baixa")
SORT = Path("data/catalunya/planta_baixa")
CAD_IDX = Path("data/catalunya/cadastre/municipis.json")
CONQUES = Path("data/catalunya/conques/conques_100km2.geojson")
REL = "https://github.com/correvents/inundacions/releases/download/cadastre-cat/"       # edificis_MMMMM.json.gz
REL_PB = "https://github.com/correvents/inundacions/releases/download/planta-baixa-cat/"  # us_*, pb_ovc_* (la de cadastre-cat és plena: 1000 fitxers)
OVC = "https://ovc.catastro.meh.es/OVCServWeb/OVCWcfCallejero/COVCCallejero.svc/json/Consulta_DNPRC"

# classes, de més a menys sensible (la planta baixa pren la més sensible dels seus immobles)
CLASSES = ["habitatge", "equipament", "comerç", "oficines", "habitatge o local", "industrial", "garatge o magatzem", "agrari", "altres"]
K = {c: i for i, c in enumerate(CLASSES)}
INSPIRE = {"1_residential": "habitatge", "2_agriculture": "agrari", "3_industrial": "industrial", "4_1_office": "oficines",
           "4_2_retail": "comerç", "4_3_publicServices": "equipament"}
# paraules de l'ús del Cadastre (luso de l'immoble o lcd de la construcció) → classe
PARAULES = [("VIVIENDA", "habitatge"), ("RESIDENCIAL", "habitatge"),
            ("APARCAMIENTO", "garatge o magatzem"), ("ESTACIONAMIENTO", "garatge o magatzem"), ("ALMACEN", "garatge o magatzem"),
            ("TRASTERO", "garatge o magatzem"), ("GARAJE", "garatge o magatzem"),
            ("COMERC", "comerç"), ("OCIO", "comerç"), ("HOSTEL", "comerç"), ("OFICINA", "oficines"), ("INDUSTRI", "industrial"),
            ("CULTURAL", "equipament"), ("DEPORTIV", "equipament"), ("SANIDAD", "equipament"), ("BENEFICENCIA", "equipament"),
            ("RELIGIOS", "equipament"), ("ESPECTACULO", "equipament"), ("SINGULAR", "equipament"), ("ENSE", "equipament"),
            ("AGRARI", "agrari"), ("AGRICOL", "agrari"), ("AGRAR", "agrari")]
NO_COMPTA = ("ELEMENTOS COMUNES", "PORCHE", "PISCINA", "DEPORTIVO DESCUBIERTO", "SUELO", "URBANIZACION", "JARDIN", "PATIO", "TERRAZA")


def classe_us(txt):
    t = C.normal(txt or "")
    if not t or any(C.normal(n) in t for n in NO_COMPTA):
        return None
    return next((c for p, c in PARAULES if C.normal(p) in t), "altres")


def nivell(pt):
    """Planta del Cadastre → 0 (baixa), -1 (sota rasant), 1 (pis) o None (desconeguda)."""
    p = (pt or "").strip().upper().lstrip("+")
    if not p:
        return None
    if p in ("00", "0", "000", "BJ", "PB", "BA", "BX", "B", "BAJ"):
        return 0
    if p.startswith("-") or p in ("ST", "SS", "SO", "SM", "S1", "S2", "S3", "SB", "SÓ", "SE"):
        return -1
    if re.fullmatch(r"\d+", p) or p in ("EN", "PR", "AT", "SA", "BC", "PL", "EP"):
        return 1
    return None


# ---------- 1. INSPIRE: ús principal, habitatges i unitats de cada edifici ----------
def inspire_municipi(url):
    z = zipfile.ZipFile(io.BytesIO(C.get(C.zip_de(url)).content))
    fedif = next(n for n in z.namelist() if n.lower().endswith("building.gml"))
    out = {}
    for _, el in etree.iterparse(z.open(fedif), tag="{*}Building"):
        ref = (el.get("{http://www.opengis.net/gml/3.2}id") or "").replace("ES.SDGC.BU.", "")
        v = {C.ln(x): (x.text or "").strip() for x in el.iter() if C.ln(x) in ("currentUse", "numberOfDwellings", "numberOfBuildingUnits")}
        el.clear()
        if ref:
            out[ref] = [v.get("currentUse", ""), int(v.get("numberOfDwellings") or 0), int(v.get("numberOfBuildingUnits") or 0)]
    return out


def ordre_inspire(a):
    TMP.mkdir(parents=True, exist_ok=True)
    ja = set(Path(a.salta).read_text().split()) if a.salta and Path(a.salta).exists() else set()
    fallats, fets = [], 0
    for prov, codi, nom, url in C_tots()[a.lot::a.lots]:
        f = TMP / f"us_{codi}.json.gz"
        if f.name in ja:
            continue
        try:
            d = inspire_municipi(url)
        except Exception as ex:
            print(f"{codi} {nom}: {type(ex).__name__} {ex}"); fallats.append(codi); continue
        with gzip.open(f, "wt", compresslevel=9) as g:
            json.dump({"camps": ["ús INSPIRE", "habitatges", "unitats"], "edificis": d}, g, separators=(",", ":"))
        fets += 1
        n = {}
        for u, *_ in d.values():
            n[u] = n.get(u, 0) + 1
        print(f"{codi} {nom}: {len(d)} edificis {n}")
    (TMP / f"fallats_{a.lot}.txt").write_text(" ".join(fallats))
    print(f"Lot {a.lot}: {fets} municipis nous, fallats: {' '.join(fallats) or 'cap'}")
    return 0


def C_tots():
    import cadastre_cat
    return cadastre_cat.tots()


# ---------- edificis d'una conca ----------
def baixa(nom, rel=REL):
    TMP.mkdir(parents=True, exist_ok=True)
    f = TMP / nom
    if f.exists():
        return f
    for k in range(4):
        try:
            r = requests.get(rel + nom, timeout=300)
            if r.status_code == 404:
                return None
            r.raise_for_status(); f.write_bytes(r.content); return f
        except requests.RequestException as e:
            print(f"  {nom}: {e}"); time.sleep(10 * (k + 1))
    raise RuntimeError(f"no puc baixar {nom}")


def edificis_conca(cid):
    """{ref: (codi municipi, plantes sobre, plantes sota)} dels edificis amb el centre dins la conca."""
    fc = json.loads(CONQUES.read_text())
    g = shape(next(f for f in fc["features"] if f["properties"]["id"] == cid)["geometry"])
    t = Transformer.from_crs(4326, 25831, always_xy=True)
    import shapely
    g31 = shapely.transform(g, lambda xy: np.column_stack(t.transform(xy[:, 0], xy[:, 1])))
    pg = prep(g31)
    x0, y0, x1, y1 = g31.bounds
    out = {}
    for codi, m in json.loads(CAD_IDX.read_text()).items():
        b = m["bbox"]
        if b[2] < x0 or b[0] > x1 or b[3] < y0 or b[1] > y1:
            continue
        f = baixa(f"edificis_{codi}.json.gz")
        if f is None:
            continue
        for ref, sobre, sota, an in json.loads(gzip.open(f, "rt").read())["parts"]:
            if sobre < 1 and sota < 1:
                continue
            xa, ya = np.mean(an[0], axis=0)
            if not (x0 <= xa <= x1 and y0 <= ya <= y1) or not pg.contains(Point(xa, ya)):
                continue
            o = out.get(ref)
            out[ref] = (codi, max(sobre, o[1] if o else 0), max(sota, o[2] if o else 0))
    print(f"Conca {cid}: {len(out)} edificis")
    return out


# ---------- 2. consulta per parcel·la ----------
def troba(d, clau):
    """Primer valor de `clau` en un JSON niat."""
    if isinstance(d, dict):
        if clau in d and not isinstance(d[clau], (dict, list)):
            return d[clau]
        for v in d.values():
            r = troba(v, clau)
            if r is not None:
                return r
    elif isinstance(d, list):
        for v in d:
            r = troba(v, clau)
            if r is not None:
                return r
    return None


def llista(x):
    return x if isinstance(x, list) else ([] if x is None else [x])


def resumeix(j):
    """Resposta de Consulta_DNPRC → [[planta, ús], ...] (un per immoble o, si és un sol immoble, per construcció)."""
    r = j.get("consulta_dnprcResult", j)
    if r.get("lerr"):
        return {"error": str(troba(r["lerr"], "des"))}
    out = []
    for u in llista(troba_obj(r, "lrcdnp", "rcdnp")):
        out.append([troba(u, "pt") or "", troba(u, "luso") or ""])
    bico = r.get("bico")
    if bico:
        luso = troba(bico.get("bi", {}), "luso") or ""
        cons = llista(troba_obj(bico, "lcons", None))
        for c in cons:
            out.append([troba(c, "pt") or "", c.get("lcd") or luso])
        if not cons:
            out.append([troba(bico.get("bi", {}), "pt") or "", luso])
    return out


def troba_obj(d, k1, k2):
    v = d.get(k1) if isinstance(d, dict) else None
    if k2 and isinstance(v, dict):
        v = v.get(k2)
    return v


MOSTRA = [0]


def consulta(ref):
    for k in range(5):
        try:
            r = requests.get(OVC, params={"RefCat": ref[:14]}, timeout=60)
            if r.status_code in (403, 429, 503):          # el Cadastre talla si anem massa de pressa: esperem
                time.sleep(30 * (k + 1)); continue
            r.raise_for_status()
            j = r.json()
            if MOSTRA[0] < 3:                       # unes quantes respostes senceres al registre, per revisar el format
                MOSTRA[0] += 1; print(f"Mostra {ref}: {json.dumps(j, ensure_ascii=False)[:1500]}")
            return ref, resumeix(j)
        except Exception as e:
            err = f"{type(e).__name__} {e}"; time.sleep(10 * (k + 1))
    return ref, {"error": err}


def ordre_ovc(a):
    ed = edificis_conca(a.conca)
    previ = {}
    if a.previ and Path(a.previ).exists():
        previ = json.loads(gzip.open(a.previ, "rt").read())
    refs = [r for r in sorted(ed) if not isinstance(previ.get(r), list)]
    if a.max:
        refs = refs[:a.max]
    print(f"{len(previ)} parcel·les ja consultades; en consulto {len(refs)}")
    t0, fets, errors, seguits = time.time(), 0, 0, 0
    with ThreadPoolExecutor(a.fils) as ex:
        for ref, res in ex.map(consulta, refs):
            previ[ref] = res; fets += 1
            if isinstance(res, dict):
                errors += 1; seguits += 1
                if errors <= 5: print(f"  Error {ref}: {res.get('error')}")
                if seguits >= 40:
                    print("40 errors seguits: el Cadastre ens talla. Ho deixo per a la propera passada"); break
            else:
                seguits = 0
            if fets % 500 == 0:
                print(f"  {fets}/{len(refs)} ({errors} errors) en {time.time() - t0:.0f} s")
            if time.time() - t0 > a.minuts * 60:
                print("Temps esgotat: la resta a la propera passada"); break
    TMP.mkdir(parents=True, exist_ok=True)
    f = TMP / f"pb_ovc_c{a.conca:03d}.json.gz"
    with gzip.open(f, "wt", compresslevel=9) as g:
        json.dump(previ, g, separators=(",", ":"), ensure_ascii=False)
    ok = sum(1 for v in previ.values() if isinstance(v, list))
    print(f"Conca {a.conca}: {ok}/{len(ed)} parcel·les amb resposta, {errors} errors en aquesta passada → {f}")
    return 0


# ---------- 3. capa per conca ----------
def classifica_ovc(units):
    pb, sota, hab, n = [], [], 0, 0
    for pt, us in units:
        c = classe_us(us)
        if c is None:
            continue
        nv = nivell(pt)
        if nv == 0:
            pb.append(c); n += 1; hab += c == "habitatge"
        elif nv == -1:
            sota.append(c)
    if not pb and not sota:
        return None
    return [0, min((K[c] for c in pb), default=-1), hab, n, min((K[c] for c in sota), default=-1)]


def estima_inspire(us, hab, unit, sobre, sota):
    c = INSPIRE.get(us)
    if c is None:
        return None
    if c == "habitatge":
        if hab <= 1 or sobre <= 1:
            pb, n_hab = "habitatge", min(hab, 1) or 1
        elif unit > hab:
            pb, n_hab = "habitatge o local", 0
        else:
            pb, n_hab = "habitatge", 1
    else:
        pb, n_hab = c, 0
    s = ("garatge o magatzem" if c == "habitatge" and hab > 1 else c) if sota > 0 else None
    return [1, K[pb] if sobre >= 1 else -1, n_hab, 1 if sobre >= 1 else 0, K[s] if s else -1]


def ordre_conca(a):
    ed = edificis_conca(a.conca)
    ovc = {}
    f = baixa(f"pb_ovc_c{a.conca:03d}.json.gz", REL_PB)
    if f is not None:
        ovc = json.loads(gzip.open(f, "rt").read())
    ins = {}
    for codi in sorted({v[0] for v in ed.values()}):
        fu = baixa(f"us_{codi}.json.gz", REL_PB)
        if fu is not None:
            ins.update(json.loads(gzip.open(fu, "rt").read())["edificis"])
    out, fonts = {}, [0, 0, 0]
    for ref, (codi, sobre, sota) in ed.items():
        v = classifica_ovc(ovc[ref]) if isinstance(ovc.get(ref), list) else None
        if v is None and ref in ins:
            v = estima_inspire(*ins[ref], sobre, sota)
        if v is None:
            fonts[2] += 1; continue
        fonts[v[0]] += 1
        out[ref] = v
    per = {c: sum(1 for v in out.values() if v[1] == i) for i, c in enumerate(CLASSES)}
    SORT.mkdir(parents=True, exist_ok=True)
    (SORT / f"c{a.conca:03d}.json").write_text(json.dumps({
        "font": "Dirección General del Catastro: consulta per parcel·la (OVC Consulta_DNPRC) i INSPIRE Buildings",
        "classes": CLASSES, "camps": ["font (0 consulta per parcel·la, 1 estimat amb INSPIRE)", "classe planta baixa (-1 cap)",
                                      "habitatges a planta baixa", "immobles a planta baixa", "classe soterrani (-1 cap)"],
        "edificis": out}, ensure_ascii=False, separators=(",", ":")))
    resum = {"edificis": len(ed), "amb_us": len(out), "consulta_parcela": fonts[0], "estimat_inspire": fonts[1], "sense_dades": fonts[2],
             "planta_baixa_per_classe": per, "habitatges_planta_baixa": sum(v[2] for v in out.values())}
    (SORT / f"c{a.conca:03d}.resum.json").write_text(json.dumps(resum, ensure_ascii=False))
    ordre_index(a)
    print(f"Conca {a.conca}: {json.dumps(resum, ensure_ascii=False)}")
    return 0


def ordre_index(a=None):
    """index.json a partir de tots els cNNN.resum.json (així dues execucions no es trepitgen)."""
    idx = {f.name[1:4].lstrip("0"): json.loads(f.read_text()) for f in sorted(SORT.glob("c*.resum.json"))}
    (SORT / "index.json").write_text(json.dumps(idx, ensure_ascii=False, indent=1))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="ordre", required=True)
    a1 = sub.add_parser("inspire"); a1.add_argument("--lot", type=int, required=True); a1.add_argument("--lots", type=int, required=True)
    a1.add_argument("--salta")
    a2 = sub.add_parser("ovc"); a2.add_argument("--conca", type=int, required=True); a2.add_argument("--max", type=int, default=0)
    a2.add_argument("--previ"); a2.add_argument("--fils", type=int, default=2); a2.add_argument("--minuts", type=float, default=300)
    a3 = sub.add_parser("conca"); a3.add_argument("--conca", type=int, required=True)
    sub.add_parser("index")
    a = ap.parse_args()
    return {"inspire": ordre_inspire, "ovc": ordre_ovc, "conca": ordre_conca, "index": ordre_index}[a.ordre](a)


if __name__ == "__main__":
    sys.exit(main())
