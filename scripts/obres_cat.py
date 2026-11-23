#!/usr/bin/env python3
"""Secció real i capacitat hidràulica de cada pas de torrentet sota via (data/catalunya/passos), a tot Catalunya.

El LiDAR no veu l'obra de drenatge (tub, calaix, pont) que hi ha sota un terraplè. Aquí se li posa la secció i la
capacitat de desguàs a partir de fonts obertes, per ordre de confiança, i cada valor diu d'on surt:

  1. ACA, Punts Singulars Estructurals (PSE) del Districte de Conca Fluvial de Catalunya: 2.619 ponts, obres de
     drenatge ≥ 2 m i cobertes de llera als trams ARPSI, amb fitxa (secció, gàlib, llum, obertures, cotes de rasant i
     de clau) i cabal màxim de desguàs calculat per l'ACA amb model 2D (IBER). Capa WFS AIGUA_PS_ESTRUCTURALS_PUBLIC
     (sig.gencat.cat/ows/AIGUA/wfs) + fitxes PDF (acanet.gencat.cat/sig/fitxes/PSE/<codi>.pdf). El servei WFS no
     accepta clients automàtics, per això les fitxes es baixen des d'un ordinador i se'n desa una instantània al repo
     (data/catalunya/obres/aca_pse.json, ordre «aca»). Si el pas és un PSE, la capacitat és la de l'ACA.
  1b. Catàleg de ponts i estructures de la xarxa de carreteres de la Generalitat (Territori, v2 2018, CC BY-NC,
     terra.repositori.gencat.cat/handle/20.500.13045/344): 3.594 ponts, pontons (llum 3-10 m) i viaductes amb
     coordenades, longitud i nombre de trams. Només els que el nom diu que passen sobre un curs d'aigua. No dona
     gàlib: l'alçada es pren com l'alçada del terraplè (MDT) i es diu. Instantània gencat_ponts.json (ordre «gencat»).
  2. OpenStreetMap: waterway amb tunnel=culvert i etiquetes diameter / width / height (ordre «osm», Geofabrik).
  3. Relació empírica secció ~ conca drenada, AJUSTADA AMB LES OBRES DE LES FONTS 1 I 2 de tot Catalunya
     (log A = a + b·log S, mínims quadrats). No hi ha cap coeficient triat a mà: si no hi ha prou obres per ajustar,
     no s'estima (null).

Capacitat (si no la dona l'ACA): control d'entrada, la norma habitual en obres curtes de drenatge, amb l'aigua que
pot pujar abans de passar per sobre de la via (HW = fondària del llac si es tapa, prof_m de passos_cat.py, MDT 5 m):
  - HW ≤ H: entrada en lámina lliure, règim crític a la boca:  Q = n·B·√g·(2/3·HW)^1.5
  - HW > H: boca en càrrega (orifici, aresta viva Cd = 0,6):  Q = n·Cd·B·H·√(2g·(HW − H/2))  (mai menys que amb HW = H)
  Un tub de diàmetre D es tracta com un calaix de la mateixa àrea (B = H = D·√π/2).
  Escenaris: obert (secció sencera), mig tapat (meitat de l'alçada útil), tapat (0, el de passos_cat.py).
  Quan hi hagi el MET 1 m de l'ICGC, HW (alçada del terraplè) es podrà refer amb decímetres.

Sortides:  data/catalunya/obres/aca_pse.json   (instantània ACA, tots els PSE, amb coordenades)
           data/catalunya/obres/cNNN.json      ({id_pas: {seccio, font, capacitat...}})
           data/catalunya/obres/index.json     (resum per font i ajust empíric)
És una OPCIÓ: passos_cat.py i la simulació no canvien si no es llegeixen aquests fitxers.

Ús:  python3 scripts/obres_cat.py aca --txt ~/Downloads/inundacions/obres/txt --wfs aca_pse_wfs.json(opcional)
     python3 scripts/obres_cat.py gencat --txt <carpeta amb les fitxes en text (pdftotext -layout)>
     python3 scripts/obres_cat.py osm --pbf cataluna-latest.osm.pbf --sortida /tmp/culverts.json.gz
     python3 scripts/obres_cat.py tot --culverts /tmp/culverts.json.gz"""
import argparse, gzip, json, math, re, sys
from pathlib import Path

PASSOS = Path("data/catalunya/passos")
SORT = Path("data/catalunya/obres")
G = 9.81
CD = 0.6                 # coeficient d'orifici d'una boca d'aresta viva (manuals d'hidràulica, no ajustat)
RIERES = ("river", "stream", "ditch", "drain", "canal", "brook")


# ---------- utilitats ----------
def num(s):
    """'2,35' / '35.4' / '> Q500_est' / '-' → float o None"""
    if s is None:
        return None
    m = re.search(r"-?\d[\d.,]*", str(s))
    if not m:
        return None
    t = m.group().rstrip(".,")
    if "," in t:                      # format de les fitxes: 1.153,22
        t = t.replace(".", "").replace(",", ".")
    try:
        return float(t)
    except ValueError:
        return None


def mida_osm(v, mm_si_gran=True):
    """Etiqueta OSM de mida → metres. '600' (sense unitat, > 20) i '600 mm' són mil·límetres; '1.2', '1,2 m' són metres."""
    if not v:
        return None
    s = str(v).strip().lower().replace(",", ".")
    m = re.match(r"^(\d+(?:\.\d+)?)\s*(mm|cm|m)?\b", s)
    if not m:
        return None
    x, u = float(m.group(1)), m.group(2)
    if u == "mm" or (u is None and mm_si_gran and x > 20):
        x /= 1000
    elif u == "cm":
        x /= 100
    return x if 0.1 <= x <= 30 else None


def capacitat(B, H, n, hw):
    """Control d'entrada (m³/s): lámina lliure si hw ≤ H, orifici si hw > H."""
    if not (B and H and hw) or hw <= 0:
        return None
    lliure = lambda h: n * B * math.sqrt(G) * (2 / 3 * h) ** 1.5
    if hw <= H:
        return lliure(hw)
    return max(n * CD * B * H * math.sqrt(2 * G * (hw - H / 2)), lliure(H))


def lonlat(x, y):
    from pyproj import Transformer
    return Transformer.from_crs(25831, 4326, always_xy=True).transform(x, y)


# ---------- 1. ACA PSE ----------
CAMPS = {  # etiqueta de la fitxa → clau
    "Nom estructura": "nom", "Tipus estructura": "tipus", "Secció de l'obertura": "seccio", "Gàlib": "galib_m",
    "Longitud del pont": "llum_m", "Amplada tauler": "amplada_tauler_m", "Nombre d'obertures": "obertures",
    "Nombre de piles": "piles", "Cota rasant via": "cota_rasant_m", "Cota clau d'obra": "cota_clau_m",
    "Carretera (via)": "via", "Tipus de via": "tipus_via", "Titularitat": "titular", "Nom riu": "riu", "Municipi": "municipi",
    "Cabal màxim de desguàs (m³/s)": "q_max_m3s", "Cabal d’afecció a rasant (m³/s)": "q_afeccio_rasant",
    "Cabal d'afecció a rasant (m³/s)": "q_afeccio_rasant", "Grau d’afecció hidràulica": "grau_afeccio",
    "Grau d'afecció hidràulica": "grau_afeccio", "Metodologia determinació afecció": "metode",
    "Necessitat de gestió del risc d'inundació": "necessitat_gestio", "Actualització": "actualitzacio",
    "Comarca": "comarca", "Codi riu": "codi_riu", "PK (km)": "pk", "Codi conca": "codi_conca", "Conca": "conca",
    "Categoria de la via": "categoria_via", "Àmbit de la via": "ambit_via", "ID_ES model hidràulic": "model_id",
    "Codi TRI": "codi_tri", "Nucli urbà": "nucli_urba", "Via única d’evacuació": "via_unica_evacuacio",
    "Diàmetre": "diametre_m", "Amplada": "amplada_m", "Alçada": "alcada_m", "Llum": "llum_m", "Longitud": "longitud_m",
}
NUMS = {"galib_m", "llum_m", "amplada_tauler_m", "obertures", "piles", "cota_rasant_m", "cota_clau_m", "q_max_m3s",
        "diametre_m", "amplada_m", "alcada_m", "longitud_m"}


def llegeix_fitxa(txt):
    t = re.sub(r"[ \t]+", " ", txt)
    etiq = sorted(CAMPS, key=len, reverse=True)
    pat = re.compile("(" + "|".join(re.escape(e) for e in etiq) + r")\s*:\s*")
    r = {}
    for linia in t.splitlines():
        ms = list(pat.finditer(linia))
        for i, m in enumerate(ms):
            v = linia[m.end(): ms[i + 1].start() if i + 1 < len(ms) else len(linia)].strip()
            k = CAMPS[m.group(1)]
            if k in r or not v:
                continue
            if k in NUMS and "Q" in v:          # «> Q500_est»: més que el cabal T500 que hi arriba, no un número
                r[k + "_supera"] = v.replace(" ", ""); continue
            r[k] = num(v) if k in NUMS else v
    if "nom" in r:  # el nom pot continuar a la línia següent: es pren de la descripció del codi
        m = re.search(r"Nom estructura:\s*(.+?)\((\d{3}E\d{4})\)", t, re.S)
        if m:
            r["nom"] = " ".join(m.group(1).split())
    m = re.search(r"Coordenades UTM:[^(\n]*\((\d+)[, ]+(\d+)\)", t)
    if m:
        r["utm"] = [int(m.group(1)), int(m.group(2))]
    m = re.search(r"Descripció:\s*(.+?)(?:\n\s*\n|Codificació)", t, re.S)
    if m:
        r["descripcio"] = " ".join(m.group(1).split())
    # taula T10/T100/T500: avinguda, estructura, làmina, capacitat
    for T in (10, 100, 500):
        m = re.search(rf"\b{T} anys\s+([\d.,-]+)\s+([\d.,-]+)\s+([\d.,-]+)\s+([A-ZÀ-Ú ]+)", t)
        if m:
            r[f"t{T}"] = {"q_avinguda": num(m.group(1)), "q_estructura": num(m.group(2)), "lamina_m": num(m.group(3)),
                          "capacitat": m.group(4).strip()}
    return r


def aca(a):
    d = Path(a.txt).expanduser()
    wfs = {}
    if a.wfs and Path(a.wfs).exists():
        for f in json.loads(Path(a.wfs).read_text())["features"]:
            wfs[f["properties"]["CODI_PSE"]] = f["properties"]
    codis = sorted(p.stem for p in d.glob("*.txt") if not p.name.startswith("._"))
    out, buits = [], []
    for c in codis:
        r = llegeix_fitxa((d / f"{c}.txt").read_text(errors="replace"))
        w = wfs.get(c, {})
        if "utm" not in r and w.get("X"):
            r["utm"] = [w["X"], w["Y"]]
        if "utm" not in r:
            buits.append(c); print("  sense coordenades:", c); continue
        r["lonlat"] = [round(v, 6) for v in lonlat(*r["utm"])]
        for k, kw in (("q_max_m3s", "Q_MAX_DESG"), ("tipus", "TIPUS_ESTR"), ("grau_afeccio", "G_AFEC_HID")):
            if r.get(k) in (None, "") and w.get(kw) not in (None, ""):
                r[k] = num(w[kw]) if k == "q_max_m3s" else w[kw]
        if r.get("q_max_m3s") is None and r.get("q_max_m3s_supera") and r.get("t500", {}).get("q_estructura"):
            r["q_max_min_m3s"] = r["t500"]["q_estructura"]      # cota inferior: passa tot el T500 que hi arriba
        r["codi"] = c
        r["fitxa"] = f"https://acanet.gencat.cat/sig/fitxes/PSE/{c}.pdf"
        out.append(r)
    SORT.mkdir(parents=True, exist_ok=True)
    (SORT / "aca_pse.json").write_text(json.dumps({"font": "ACA, Punts Singulars Estructurals (WFS AIGUA_PS_ESTRUCTURALS_PUBLIC i fitxes PSE)",
                                                    "llicencia": "Dades obertes Generalitat de Catalunya (citar ACA)",
                                                    "n": len(out), "pse": out}, ensure_ascii=False, separators=(",", ":")))
    from collections import Counter
    print(f"ACA PSE: {len(out)} fitxes llegides, {len(buits)} sense coordenades")
    print(" tipus:", Counter(r.get("tipus") for r in out).most_common())
    print(" secció:", Counter(r.get("seccio") for r in out).most_common(12))
    print(" amb gàlib:", sum(r.get("galib_m") is not None for r in out), " amb llum:", sum(r.get("llum_m") is not None for r in out),
          " amb q_max:", sum(r.get("q_max_m3s") is not None for r in out))


# ---------- 1b. Catàleg de ponts i estructures de la xarxa de carreteres de la Generalitat ----------
AIGUA = re.compile(r"\b(riu|riera|rieral|torrent|torrentet|rec|barranc|rial|sèquia|séquia|canal|font|clot|sot|fondo|xaragall|"
                   r"reguer|regatell|llera|riereta|aigua|drenatge|desguàs|bassa|estany|ribera|curs|gorg|rasa|cuneta)\b", re.I)


def gencat(a):
    """Llegeix les fitxes del Catàleg (pdftotext -layout) → data/catalunya/obres/gencat_ponts.json."""
    out = []
    for fp in sorted(Path(a.txt).expanduser().glob("*.txt")):
        t = fp.read_text(errors="replace")
        for bloc in t.split("Nom de l'estructura:")[1:]:
            b = re.sub(r"[ \t]+", " ", bloc)
            linies = b.splitlines()
            r = {"nom": linies[0].strip()}
            def camp(et, nums=False):
                m = re.search(rf"^\s*{et}\s+(.+?)\s*$", b, re.M)
                return (num(m.group(1)) if nums else m.group(1).strip()) if m else None
            r["carretera"] = camp("Carretera"); r["pk"] = camp("PK"); r["municipi"] = camp("Municipi")
            r["titular"] = camp("Titularitat"); r["tipologia"] = camp("Tipologia general")
            r["longitud_m"] = camp(r"Longitud total \(m\)", True); r["superficie_m2"] = camp(r"Superfície total \(m2\)", True)
            r["trams"] = camp("Nombre de trams", True); r["material"] = (camp("Material") or "").replace("m ", "m").strip() or None
            m = re.search(r"(Convencional[^\n]*\n?[^\n]*?(?:bigues|calaix)\)|Convencional)", b)
            if "calaix" in b[:2000]:
                r["tauler"] = "llosa/calaix"
            elif "bigues" in b[:2000]:
                r["tauler"] = "bigues"
            x, y = camp("UTM X", True), camp("UTM Y", True)
            m = re.search(r"Generalitat de Catalunya\s+(\S+-\S+)\s*$", b, re.M)
            r["codi"] = m.group(1) if m else None
            if not (x and y and 250000 < x < 550000 and 4480000 < y < 4760000):
                continue
            r["utm"] = [round(x, 1), round(y, 1)]
            r["lonlat"] = [round(v, 6) for v in lonlat(x, y)]
            r["sobre_aigua"] = bool(AIGUA.search(r["nom"]))
            r["font_fitxer"] = fp.stem
            out.append(r)
    SORT.mkdir(parents=True, exist_ok=True)
    (SORT / "gencat_ponts.json").write_text(json.dumps({
        "font": "Catàleg de ponts i estructures de la xarxa de carreteres de la Generalitat de Catalunya (v2, 2018, Territori)",
        "url": "https://terra.repositori.gencat.cat/handle/20.500.13045/344", "llicencia": "CC BY-NC 4.0",
        "n": len(out), "ponts": out}, ensure_ascii=False, separators=(",", ":")))
    from collections import Counter
    print(f"Catàleg Generalitat: {len(out)} estructures, {sum(r['sobre_aigua'] for r in out)} sobre un curs d'aigua (pel nom)")
    print(" tipologia:", Counter(r.get("tipologia") for r in out).most_common())


# ---------- 2. OSM ----------
def osm(pbf, sortida):
    import osmium

    class H(osmium.SimpleHandler):
        def __init__(self):
            super().__init__(); self.o = []

        def way(self, w):
            t = w.tags
            if t.get("waterway") not in RIERES or not (t.get("tunnel") in ("culvert", "yes") or t.get("culvert") == "yes"):
                return
            et = {k: t.get(k) for k in ("diameter", "width", "height", "culvert:diameter", "culvert:width", "culvert:height",
                                        "est_width", "material", "tunnel", "layer") if t.get(k)}
            try:
                xy = [v for n in w.nodes for v in (round(n.lon, 6), round(n.lat, 6))]
            except osmium.InvalidLocationError:
                return
            if len(xy) >= 4:
                self.o.append([w.id, t.get("waterway"), et, xy])

    h = H(); h.apply_file(pbf, locations=True, idx="flex_mem")
    with gzip.open(sortida, "wt") as f:
        json.dump(h.o, f, separators=(",", ":"))
    amb = sum(1 for o in h.o if any(k in o[2] for k in ("diameter", "width", "culvert:diameter", "culvert:width")))
    print(f"OSM: {len(h.o)} obres de drenatge (tunnel=culvert), {amb} amb mida")


def seccio_osm(et):
    D = mida_osm(et.get("diameter") or et.get("culvert:diameter"))
    if D:
        return {"forma": "circular", "D_m": round(D, 2), "B_m": round(D * math.sqrt(math.pi) / 2, 2),
                "H_m": round(D * math.sqrt(math.pi) / 2, 2), "n": 1}
    B = mida_osm(et.get("culvert:width") or et.get("width"), mm_si_gran=False)
    H = mida_osm(et.get("culvert:height") or et.get("height"), mm_si_gran=False)
    if B:
        return {"forma": "rectangular" if H else "rectangular (alçada desconeguda)", "B_m": round(B, 2),
                "H_m": round(H, 2) if H else None, "n": 1}
    return None


# ---------- 3. assignació i capacitat ----------
def tot(a):
    from pyproj import Transformer
    from shapely.geometry import Point, LineString
    from shapely.strtree import STRtree
    import numpy as np
    utm = Transformer.from_crs(4326, 25831, always_xy=True).transform

    pse = json.loads((SORT / "aca_pse.json").read_text())["pse"] if (SORT / "aca_pse.json").exists() else []
    pse_pt = [Point(r["utm"]) for r in pse]
    arbre_p = STRtree(pse_pt) if pse_pt else None
    gp = []
    if (SORT / "gencat_ponts.json").exists():
        gp = [r for r in json.loads((SORT / "gencat_ponts.json").read_text())["ponts"] if r["sobre_aigua"]]
    gp_pt = [Point(r["utm"]) for r in gp]
    cul = []
    if a.culverts:
        with gzip.open(a.culverts, "rt") as f:
            for wid, ww, et, xy in json.load(f):
                s = seccio_osm(et)
                if s:
                    cul.append((wid, s, LineString([utm(xy[k], xy[k + 1]) for k in range(0, len(xy), 2)])))
    arbre_c = STRtree([c[2] for c in cul]) if cul else None
    print(f"ACA PSE: {len(pse)}  ·  Catàleg Generalitat sobre aigua: {len(gp)}  ·  OSM culverts amb mida: {len(cul)}")

    # tots els passos (punts) de Catalunya
    pas = []                                                   # (cid, props, punt UTM)
    for fp in sorted(PASSOS.glob("c*.geojson")):
        cid = int(fp.stem[1:])
        for f in json.loads(fp.read_text())["features"]:
            if "id" in f["properties"] and f["geometry"]["type"] == "Point":
                pas.append((cid, f["properties"], Point(utm(*f["geometry"]["coordinates"]))))
    arbre_pas = STRtree([x[2] for x in pas])
    # cada obra es casa amb UN sol pas: el PSE de l'ACA amb el pas de més conca a menys de dist-aca (terraplè abans
    # que gual: el PSE és al riu principal, no al torrentet del costat); el culvert d'OSM amb el pas més proper.
    de_pse, de_gp, de_cul = {}, {}, {}
    for j, pt in enumerate(pse_pt):
        ii = arbre_pas.query(pt.buffer(a.dist_aca))
        if len(ii):
            k = max(ii, key=lambda k: (pas[k][1].get("tipus") == "terraplè", pas[k][1].get("conca_km2") or 0))
            if k not in de_pse or pt.distance(pas[k][2]) < pse_pt[de_pse[k]].distance(pas[k][2]):
                de_pse[k] = j
    for j, pt in enumerate(gp_pt):
        ii = arbre_pas.query(pt.buffer(a.dist_aca))
        if len(ii):
            k = max(ii, key=lambda k: (pas[k][1].get("tipus") == "terraplè", pas[k][1].get("conca_km2") or 0))
            if k not in de_gp or pt.distance(pas[k][2]) < gp_pt[de_gp[k]].distance(pas[k][2]):
                de_gp[k] = j
    for j, (wid, s_, g) in enumerate(cul):
        ii = arbre_pas.query(g.buffer(a.dist_osm))
        if len(ii):
            k = min(ii, key=lambda k: g.distance(pas[k][2]))
            if k not in de_cul or g.distance(pas[k][2]) < cul[de_cul[k]][2].distance(pas[k][2]):
                de_cul[k] = j

    res, mostres = {}, []
    for k, (cid, p, pt) in enumerate(pas):
        d = res.setdefault(cid, {})
        hw = p.get("prof_m") if p.get("tipus") == "terraplè" else None
        e = {"tipus_pas": p.get("tipus"), "conca_km2": p.get("conca_km2"), "hw_m": hw,
             "hw_font": "fondària del llac si es tapa (MDT 5 m IGN, metres enters)" if hw else None}
        # 1. ACA
        if k in de_pse:
            j = de_pse[k]; r = pse[j]
            n = int(r.get("obertures") or 1)
            B = (r.get("llum_m") or r.get("amplada_m") or r.get("diametre_m"))
            H = r.get("galib_m") or r.get("alcada_m") or r.get("diametre_m")
            if B and r.get("amplada_tauler_m") and (r.get("tipus") or "").startswith(("Obra", "Cobertura")):
                # en obres i cobertes la fitxa a vegades dona la «longitud» al llarg del curs i l'«amplada» a través
                # (24 de 248): la llum és la més petita de les dues (una obra sota una via és més llarga que ampla)
                B = min(B, r["amplada_tauler_m"])
            if B and n > 1:
                B = B / n                                       # llum total repartida entre obertures (piles no restades)
            e.update({"font": "ACA PSE", "ref": r["codi"], "fitxa": r["fitxa"], "dist_m": round(pse_pt[j].distance(pt)),
                      "estructura": r.get("tipus"), "forma": r.get("seccio"), "B_m": B, "H_m": H, "n": n,
                      "q_obert_m3s": r.get("q_max_m3s") or r.get("q_max_min_m3s"),
                      "q_obert_es_minim": bool(r.get("q_max_min_m3s") and not r.get("q_max_m3s")) or None,
                      "q_font": f"ACA, fitxa PSE ({r.get('metode') or 'model hidràulic'})" if (r.get("q_max_m3s") or r.get("q_max_min_m3s")) else None,
                      "q_afeccio_rasant": r.get("q_afeccio_rasant"), "grau_afeccio": r.get("grau_afeccio")})
            if H and B and r.get("cota_rasant_m") and r.get("cota_clau_m"):
                e["hw_m"] = round(r["cota_rasant_m"] - r["cota_clau_m"] + H, 2)
                e["hw_font"] = "cota rasant − cota clau + gàlib (fitxa ACA)"
        # 1b. Catàleg de ponts de la Generalitat: llum i trams; sense alçada (gàlib) a la fitxa
        elif k in de_gp:
            r = gp[de_gp[k]]
            n = int(r.get("trams") or 1)
            B = r["longitud_m"] / n if r.get("longitud_m") else None
            H = hw                                              # la fitxa no dona gàlib: com a màxim, l'alçada del terraplè
            e.update({"font": "Generalitat, catàleg de ponts", "ref": r.get("codi"), "nom": r.get("nom"), "dist_m": round(gp_pt[de_gp[k]].distance(pt)),
                      "estructura": r.get("tipologia"), "tauler": r.get("tauler"), "material": r.get("material"),
                      "B_m": round(B, 2) if B else None, "H_m": H, "n": n,
                      "H_font": "alçada del terraplè (MDT 5 m): la fitxa no dona gàlib" if H else None})
        # 2. OSM
        elif k in de_cul:
            wid, s_, g = cul[de_cul[k]]
            e.update({"font": "OSM", "ref": f"way/{wid}", "dist_m": round(g.distance(pt)), **s_})
        if e.get("font") and e.get("B_m") and e.get("H_m") and p.get("conca_km2") and p.get("tipus") == "terraplè":
            es_drenatge = e["font"] == "OSM" or (e.get("estructura") or "").startswith(("Obra de drenatge", "Cobertura"))
            if es_drenatge:
                mostres.append((p["conca_km2"], e["B_m"] * e["H_m"] * e.get("n", 1), e["font"]))
        d[p["id"]] = e

    # 3. ajust empíric secció ~ conca, amb les obres de drenatge trobades (no ponts)
    aj = None
    if len(mostres) >= a.min_mostres:
        S = np.log([m[0] for m in mostres]); A = np.log([m[1] for m in mostres])
        b, a0 = np.polyfit(S, A, 1)
        resid = A - (a0 + b * S)
        aj = {"n": len(mostres), "n_aca": sum(m[2] == "ACA PSE" for m in mostres), "n_osm": sum(m[2] == "OSM" for m in mostres), "a": round(float(a0), 4), "b": round(float(b), 4), "r2": round(float(1 - resid.var() / A.var()), 3),
              "p25": round(float(np.percentile(resid, 25)), 3), "p75": round(float(np.percentile(resid, 75)), 3),
              "formula": "A_m2 = exp(a) · S_km2^b"}
        print(f"Ajust empíric: {aj}")
    else:
        print(f"Ajust empíric: només {len(mostres)} obres amb mida i conca (< {a.min_mostres}); no s'estima")

    cnt = {"ACA PSE": 0, "Generalitat, catàleg de ponts": 0, "OSM": 0, "empíric": 0, "sense": 0}
    SORT.mkdir(parents=True, exist_ok=True)
    for cid, d in res.items():
        for pid, e in d.items():
            if "font" not in e:
                if aj and e.get("conca_km2") and e["tipus_pas"] == "terraplè":
                    A = math.exp(aj["a"]) * e["conca_km2"] ** aj["b"]
                    s = math.sqrt(A); H = min(s, e["hw_m"]) if e.get("hw_m") else s
                    e.update({"font": "empíric", "ref": "ajust secció~conca amb obres ACA+OSM (index.json)", "forma": "equivalent quadrada",
                              "B_m": round(A / H, 2), "H_m": round(H, 2), "n": 1,
                              "A_interval_m2": [round(A * math.exp(aj["p25"]), 2), round(A * math.exp(aj["p75"]), 2)],
                              "confianca": "baixa" if aj["r2"] < 0.5 or aj["n_osm"] < 30 else "mitjana",
                              "avis": "ajust amb obres de les fonts 1-2; si gairebé totes són ACA ≥ 2 m, sobreestima les obres petites"
                              if aj["n_osm"] < 30 else None})
                else:
                    e["font"] = None
            if e.get("q_obert_m3s") is None and e.get("B_m") and e.get("H_m") and e.get("hw_m"):
                q = capacitat(e["B_m"], e["H_m"], e.get("n", 1), e["hw_m"])
                e["q_obert_m3s"] = round(q, 2) if q else None
                e["q_font"] = "control d'entrada (boca lliure / orifici Cd 0,6) amb hw"
            if e.get("B_m") and e.get("H_m") and e.get("hw_m"):
                qm = capacitat(e["B_m"], e["H_m"] / 2, e.get("n", 1), e["hw_m"])
                e["q_mig_tapat_m3s"] = round(qm, 2) if qm else None
            e = {k: (round(v, 2) if isinstance(v, float) else v) for k, v in e.items() if v is not None}
            d[pid] = e
            cnt[e.get("font") or "sense"] += 1
        (SORT / f"c{cid:03d}.json").write_text(json.dumps(d, ensure_ascii=False, separators=(",", ":")))
    (SORT / "index.json").write_text(json.dumps({"passos_per_font": cnt, "ajust_empiric": aj, "aca_pse": len(pse), "gencat_ponts_sobre_aigua": len(gp),
                                                 "osm_culverts_amb_mida": len(cul)}, ensure_ascii=False, indent=1))
    print("Passos per font:", cnt)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="ordre", required=True)
    s = sub.add_parser("aca"); s.add_argument("--txt", required=True); s.add_argument("--wfs")
    s = sub.add_parser("gencat"); s.add_argument("--txt", required=True)
    s = sub.add_parser("osm"); s.add_argument("--pbf", required=True); s.add_argument("--sortida", required=True)
    s = sub.add_parser("tot"); s.add_argument("--culverts")
    s.add_argument("--dist-aca", type=float, default=50, help="distància màxima pas ↔ PSE de l'ACA (m)")
    s.add_argument("--dist-osm", type=float, default=30, help="distància màxima pas ↔ culvert d'OSM (m, com passos_cat.py)")
    s.add_argument("--min-mostres", type=int, default=30, help="obres amb mida mínimes per fer l'ajust empíric")
    a = ap.parse_args()
    return {"aca": aca, "gencat": gencat, "tot": tot}[a.ordre](a) if a.ordre != "osm" else osm(a.pbf, a.sortida)


if __name__ == "__main__":
    sys.exit(main())
