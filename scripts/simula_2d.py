#!/usr/bin/env python3
"""Simulació 2D de la pluja sobre el relleu (rain-on-grid): on corre i on s'embassa l'aigua a cada moment.

Esquema inercial local de Bates et al. (2010), el de LISFLOOD-FP: a cada cara entre cel·les el cabal unitari q s'actualitza
amb el pendent de la superfície de l'aigua i la fricció de Manning, i la làmina h amb la divergència de q.
  - Relleu: la reixa de relleu.py (MDT 5 m de l'IGN) agregada a --pas metres.
  - Pluja: graella del radar OPERA de l'episodi (15 min, 1 km), corregida hora a hora perquè la mitjana sobre la conca
    de l'aforador coincideixi amb el pluviòmetre de l'ACA (el radar infravalora els nuclis convectius).
  - Cobertes del sòl: mapa de l'ICGC (scripts/cobertes.py, Release cobertes-cat). Cada coberta té la seva n de Manning
    (Chow 1959) i un número de corba (USDA TR-55); --ia-natural i --perdua-natural són la mitjana del terreny natural
    i es reparteixen per cel·la en proporció a la retenció S = 25400/CN - 254. Amb --sense-cobertes: ESA WorldCover.
  - Pèrdues: llindar per cel·la. Terreny natural: --ia-natural mm de mitjana (calibrat amb l'aforador el 3-4/10/2026).
    Zones construïdes: --ia-urba mm.
    Després del llindar, la làmina d'aigua de cada cel·la perd --perdua-natural mm/h (infiltració) al terreny natural i
    --clavegueram mm/h a les zones construïdes (embornals), excepte on el calat supera --h-llera (llera saturada). Així l'aigua embassada també es buida quan para de ploure.
  - Clavegueram (--escenari-clavegueram): «fix» = --clavegueram mm/h a tot arreu i l'aigua desapareix (com fins ara).
    «funciona» / «antic» / «collapsat»: capacitat estimada per subconca de clavegueram (scripts/clavegueram_cat.py,
    Release clavegueram-cat (rasters i taules): C·I de la pluja de 10 anys, de 2 anys, o la de 10 anys × --factor-collapsat si els
    embornals es tapen). L'aigua que entra a la xarxa surt al seu punt de desguàs (llera, mar o terreny) i hi segueix
    corrent. Amb --clavegueram-ofegat M, una subconca deixa de desguassar mentre el calat al desguàs passa de M metres.
  - Fricció: n de Manning per coberta (o --n-natural, --n-urba amb --sense-cobertes).
  - Mar i vores del domini: l'aigua surt.
  - Passos de torrentets sota vies (--passos, opció): per defecte tapats, com surten a l'MDT (el terraplè fa de dic).
    Amb oberts, cada pas de scripts/passos_cat.py és un conducte amb la seva capacitat; amb arrossegalls, a més, es tapa
    total o parcialment segons el risc d'arrossegalls (ARR 2019, scripts/arrossegalls.py), també als ponts. Vegeu scripts/passos_2d.py.
Sortides (data/<conca>/sim2d/<id>/): un PNG de calat cada --cada minuts en lat/lon, calat màxim, l'hidrograma simulat a
l'aforador (per validar-lo amb l'observat), index.json i, per al model d'inundació, dues graelles numèriques en UTM:
calat màxim (cm) i minut en què cada cel·la supera 15 cm d'aigua.
Opcions activables (per provar amb episodis quines ajuden):
  --edificis-diposit  els edificis són dipòsits: l'aigua hi entra per portals i rampes de garatge (abocador) i omple
                      soterranis i plantes baixes (scripts/diposit.py); surt edificis_diposit.json i el volum retingut.
  --velocitat         velocitat màxima de l'aigua i classe de perill AIDR H1-H6 (calat, velocitat, calat × velocitat).

Ús:  python3 scripts/simula_2d.py --config config/ridaura.json --mdt /tmp/relleu/mdt.tif \\
         --radar data/ridaura/episodis/2026-10-03_radar.json --episodi data/ridaura/episodis/2026-10-04.json"""
import argparse, json, sys, time
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject
from rasterio.transform import from_origin
from pyproj import Transformer
from numba import njit
from PIL import Image

G = 9.81
WORLDCOVER = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N39E{:03d}_Map.tif"


@njit(cache=True)
def pas_temps(z, h, qx, qy, n, dx, dt, actiu, psx, psy, phi):
    """Un pas de l'esquema inercial local. qx a la cara est de cada cel·la, qy a la cara sud.
    Porositat dels edificis: q és el cabal per metre d'amplada lliure; psx/psy, la fracció lliure de cada cara, i phi, la
    fracció lliure de cada cel·la (h és el calat a la part lliure, el carrer). Sense edificis tot val 1."""
    nf, nc = z.shape
    for i in range(nf):
        for j in range(nc - 1):
            if not (actiu[i, j] or actiu[i, j + 1]) or psx[i, j] == 0.0:
                qx[i, j] = 0.0; continue
            eta1, eta2 = z[i, j] + h[i, j], z[i, j + 1] + h[i, j + 1]
            hf = max(eta1, eta2) - max(z[i, j], z[i, j + 1])
            if hf <= 1e-3:
                qx[i, j] = 0.0; continue
            nn = 0.5 * (n[i, j] + n[i, j + 1])
            s = (eta2 - eta1) / dx
            q = (qx[i, j] - G * hf * dt * s) / (1 + G * dt * nn * nn * abs(qx[i, j]) / hf ** (7 / 3))
            qmax = hf * np.sqrt(G * hf)                   # limitador de Froude 1 per estabilitat
            qx[i, j] = min(max(q, -qmax), qmax)
    for i in range(nf - 1):
        for j in range(nc):
            if not (actiu[i, j] or actiu[i + 1, j]) or psy[i, j] == 0.0:
                qy[i, j] = 0.0; continue
            eta1, eta2 = z[i, j] + h[i, j], z[i + 1, j] + h[i + 1, j]
            hf = max(eta1, eta2) - max(z[i, j], z[i + 1, j])
            if hf <= 1e-3:
                qy[i, j] = 0.0; continue
            nn = 0.5 * (n[i, j] + n[i + 1, j])
            s = (eta2 - eta1) / dx
            q = (qy[i, j] - G * hf * dt * s) / (1 + G * dt * nn * nn * abs(qy[i, j]) / hf ** (7 / 3))
            qmax = hf * np.sqrt(G * hf)
            qy[i, j] = min(max(q, -qmax), qmax)
    for i in range(nf):
        for j in range(nc):
            if not actiu[i, j]:
                h[i, j] = 0.0; continue                   # mar o fora del domini: l'aigua surt
            ent = 0.0
            if j > 0: ent += psx[i, j - 1] * qx[i, j - 1]
            if j < nc - 1: ent -= psx[i, j] * qx[i, j]
            if i > 0: ent += psy[i - 1, j] * qy[i - 1, j]
            if i < nf - 1: ent -= psy[i, j] * qy[i, j]
            h[i, j] = max(h[i, j] + dt * ent / (dx * phi[i, j]), 0.0)


@njit(cache=True)
def velocitat(h, qx, qy, actiu, vmax, hvmax, perill, hmin):
    """Velocitat a cada cel·la (mitjana dels cabals unitaris de les quatre cares / calat) i classe de perill AIDR H1-H6
    (Australian Disaster Resilience Handbook 7, Smith et al. 2014) segons calat, velocitat i calat × velocitat. Desa els màxims."""
    nf, nc = h.shape
    for i in range(nf):
        for j in range(nc):
            hh = h[i, j]
            if not actiu[i, j] or hh < hmin:
                continue
            ux = 0.5 * ((qx[i, j - 1] if j > 0 else 0.0) + qx[i, j]) / hh
            uy = 0.5 * ((qy[i - 1, j] if i > 0 else 0.0) + qy[i, j]) / hh
            v = np.sqrt(ux * ux + uy * uy)
            hv = hh * v
            if v > vmax[i, j]: vmax[i, j] = v
            if hv > hvmax[i, j]: hvmax[i, j] = hv
            if hh <= 0.3 and v <= 2.0 and hv <= 0.3: c = 1
            elif hh <= 0.5 and v <= 2.0 and hv <= 0.6: c = 2
            elif hh <= 1.2 and v <= 2.0 and hv <= 0.6: c = 3
            elif hh <= 2.0 and v <= 2.0 and hv <= 1.0: c = 4
            elif hh <= 4.0 and v <= 4.0 and hv <= 4.0: c = 5
            else: c = 6
            if c > perill[i, j]: perill[i, j] = c


@njit(cache=True)
def hmax_actiu(h, actiu):
    m = 0.0
    for i in range(h.shape[0]):
        for j in range(h.shape[1]):
            if actiu[i, j] and h[i, j] > m:
                m = h[i, j]
    return m


def obstacles(fitxer, tr, nf, nc, dx, fi=2.0):
    """Porositat dels edificis del Cadastre a la reixa de la simulació (Sanders et al. 2008, porositat integral).
    Les parts amb plantes sobre rasant es rasteritzen a fi metres; de cada cel·la en surt la fracció lliure (phi) i de cada cara
    la fracció lliure per on pot passar l'aigua (psx, psy), limitada a la de les dues cel·les perquè l'esquema sigui estable.
    Cel·les gairebé tot edifici (< 10 % lliure): no hi entra aigua i la pluja de les teulades va a la cel·la lliure més propera."""
    import gzip, rasterio.features
    from scipy.ndimage import distance_transform_edt
    from shapely.geometry import Polygon
    ed = json.loads(gzip.open(fitxer, "rt").read())
    f = max(int(round(dx / fi)), 1)
    ftr = rasterio.transform.from_origin(tr.c, tr.f, dx / f, dx / f)
    formes = [(Polygon(an[0], an[1:]), 1) for ref, sobre, sota, an in ed["parts"] if sobre >= 1 and len(an[0]) >= 4]
    fina = rasterio.features.rasterize(formes, out_shape=(nf * f, nc * f), transform=ftr, fill=0, dtype="uint8").astype(bool)
    phi = 1 - fina.reshape(nf, f, nc, f).mean(axis=(1, 3))
    vx = fina[:, f - 1:-1:f][:, :nc - 1] | fina[:, f::f][:, :nc - 1]          # columnes fines a banda i banda de cada cara est
    psx = np.ones((nf, nc)); psx[:, :nc - 1] = 1 - vx.reshape(nf, f, nc - 1).mean(axis=1)
    vy = fina[f - 1:-1:f, :][:nf - 1] | fina[f::f, :][:nf - 1]
    psy = np.ones((nf, nc)); psy[:nf - 1, :] = 1 - vy.reshape(nf - 1, nc, f).mean(axis=2)
    del fina, vx, vy
    bloc = phi < 0.1
    psx[:, :-1] = np.minimum(psx[:, :-1], np.minimum(phi[:, :-1], phi[:, 1:])); psy[:-1] = np.minimum(psy[:-1], np.minimum(phi[:-1], phi[1:]))
    psx[:, :-1][bloc[:, :-1] | bloc[:, 1:]] = 0; psy[:-1][bloc[:-1] | bloc[1:]] = 0
    phi = np.where(bloc, 1.0, phi)                 # inert: no hi entra ni hi plou
    _, (ri, ci) = distance_transform_edt(bloc, return_indices=True)
    desti = (ri * nc + ci).ravel()
    ed["ocupacio_pct"] = round(100 * float(1 - np.where(bloc, 0, phi).mean()), 2)
    print(f"Edificis: {len(formes)} parts, ocupen {ed['ocupacio_pct']} % del domini; {bloc.sum()} cel·les tancades")
    return phi, psx, psy, desti, ed


def soterranis(ed, tr, nf, nc, dx, hmax, t15, durada, z, z_mdt, tr_mdt, toca=5.0, tria=None):
    """Edificis amb plantes sota rasant que tenen aigua a tocar: a les cel·les a menys de `toca` m de la planta de les parts amb
    soterrani, l'alçada màxima de l'aigua per sobre del terra de l'edifici (mediana del MDT sota la planta), sense passar del
    calat. Així l'aigua d'una riera encaixada al costat, més baixa que l'edifici, no compta. I minut en què hi passa de 15 cm."""
    from shapely.geometry import Polygon, box
    from shapely.ops import unary_union
    tria = tria or (lambda ref, sobre, sota: sota > 0)       # per defecte, les parts amb soterrani
    per_ref = {}
    for ref, sobre, sota, an in ed["parts"]:
        if tria(ref, sobre, sota) and len(an[0]) >= 4:
            p = Polygon(an[0], an[1:]).buffer(0)
            if not p.is_empty:
                per_ref.setdefault(ref, []).append(p)
    out = {}
    for ref, pols in per_ref.items():
        pl = unary_union(pols)
        if pl.is_empty:
            continue
        g = pl.buffer(toca)
        mx = ~tr_mdt
        bx0, by0, bx1, by1 = pl.bounds
        a0, b0 = mx * (bx0, by1); a1, b1 = mx * (bx1, by0)
        tros = z_mdt[max(int(b0), 0):max(int(b1) + 1, int(b0) + 1), max(int(a0), 0):max(int(a1) + 1, int(a0) + 1)]
        terra = float(np.nanmedian(tros)) if tros.size and np.isfinite(tros).any() else -1e9
        x0, y0, x1, y1 = g.bounds
        c0, c1 = max(int((x0 - tr.c) // dx), 0), min(int((x1 - tr.c) // dx), nc - 1)
        r0, r1 = max(int((tr.f - y1) // dx), 0), min(int((tr.f - y0) // dx), nf - 1)
        cm, mn, fi = 0, None, 0
        for r in range(r0, r1 + 1):
            for c in range(c0, c1 + 1):
                cel = box(tr.c + c * dx, tr.f - (r + 1) * dx, tr.c + (c + 1) * dx, tr.f - r * dx)
                if not g.intersects(cel):
                    continue
                sobre = min(hmax[r, c], z[r, c] + hmax[r, c] - terra)
                if sobre <= 0:
                    continue
                cm = max(cm, int(round(sobre * 100)))
                if np.isfinite(t15[r, c]) and sobre >= 0.15:
                    mn = t15[r, c] if mn is None else min(mn, t15[r, c])
                    if np.isfinite(durada[r, c]): fi = max(fi, t15[r, c] + durada[r, c])   # l'última cel·la del voltant que s'asseca
        if cm >= 5:
            out[ref] = [cm, None if mn is None else int(mn), 0 if mn is None else int(max(fi - mn, 0))]
    return len(per_ref), out


def plantes_baixes(fitxers, ed, dest, ini, tr, nf, nc, dx, hmax, t15, durada, z, z_mdt, tr_mdt, llindar=15):
    """OPCIÓ --planta-baixa: plantes baixes amb aigua a tocar (mateix criteri que els soterranis) i de quin ús són."""
    pb, classes = {}, None
    for f in fitxers:
        j = json.loads(Path(f).read_text())
        classes = j["classes"]; pb.update(j["edificis"])
    _, aigua = soterranis(ed, tr, nf, nc, dx, hmax, t15, durada, z, z_mdt, tr_mdt,
                          tria=lambda ref, sobre, sota: sobre >= 1 and ref in pb and pb[ref][1] >= 0)
    out = {ref: v + [pb[ref][1], pb[ref][2], pb[ref][0]] for ref, v in aigua.items()}
    risc = {ref: v for ref, v in out.items() if v[0] >= llindar}
    per = {c: sum(1 for v in risc.values() if v[3] == i) for i, c in enumerate(classes)}
    resum = {"fitxer": "plantes_baixes.json", "llindar_cm": llindar, "edificis_amb_us": len(pb), "amb_aigua_a_tocar": len(out),
             "risc": len(risc), "risc_per_classe": per, "habitatges_planta_baixa_risc": sum(v[4] for v in risc.values()),
             "risc_amb_us_de_consulta_parcela": sum(1 for v in risc.values() if v[5] == 0)}
    (dest / "plantes_baixes.json").write_text(json.dumps({"inici_utc": ini.strftime("%Y-%m-%dT%H:%M"), "toca_m": 5, "llindar_cm": llindar,
        "classes": classes, "camps": ["aigua a tocar per sobre del terra (cm)", "minut en què passa de 15 cm", "minuts fins que baixa de 5 cm",
                                      "classe planta baixa", "habitatges a planta baixa", "font de l'ús (0 consulta per parcel·la, 1 estimat)"],
        "edificis": out}, ensure_ascii=False, separators=(",", ":")))
    print(f"Plantes baixes: {len(risc)} amb ≥{llindar} cm a tocar: {per}")
    return resum


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--mdt", required=True)
    ap.add_argument("--radar", required=True)
    ap.add_argument("--episodi", help="JSON amb pluviòmetre i radar a la conca (episodi.py) per corregir el radar; sense, radar tal qual")
    ap.add_argument("--id", help="nom de la carpeta de sortida (per defecte, la data de l'últim fotograma)")
    ap.add_argument("--pas", type=float, default=20)
    ap.add_argument("--ia-natural", type=float, default=150)
    ap.add_argument("--ia-urba", type=float, default=5)
    ap.add_argument("--perdua-natural", type=float, default=12, help="infiltració un cop superat el llindar (mm/h)")
    ap.add_argument("--clavegueram-desti", choices=["perdut", "riu"], default="perdut",
                    help="on va l'aigua dels embornals: a la riera més propera amb conca ≥ --riu-km2 (riu) o fora del model (perdut)")
    ap.add_argument("--riu-km2", type=float, default=1.0, help="conca mínima de la riera que rep el clavegueram (km²)")
    ap.add_argument("--perdua-sobre", choices=["calat", "pluja"], default="calat",
                    help="d'on es resta --perdua-natural: de l'aigua que corre per cada cel·la (calat) o de la pluja que hi cau (pluja)")
    ap.add_argument("--h-llera", type=float, default=0.1, help="calat a partir del qual la cel·la és llera saturada i no infiltra (m)")
    ap.add_argument("--clavegueram", type=float, default=30, help="capacitat dels embornals a les zones construïdes (mm/h)")
    ap.add_argument("--escenari-clavegueram", default="fix", choices=["fix", "funciona", "antic", "collapsat"],
                    help="fix: --clavegueram a tot arreu; funciona/antic/collapsat: xarxa estimada (clavegueram_cat.py) amb T10, T2 o T10 × --factor-collapsat")
    ap.add_argument("--factor-collapsat", type=float, default=0.25, help="fracció de la capacitat T10 que queda amb els embornals tapats")
    ap.add_argument("--clavegueram-ofegat", type=float, default=0, help="calat al punt de desguàs (m) a partir del qual la xarxa no desguassa (0 = mai)")
    ap.add_argument("--n-natural", type=float, default=0.06, help="n de Manning del terreny natural (només amb --sense-cobertes)")
    ap.add_argument("--n-urba", type=float, default=0.02, help="n de Manning de les zones construïdes (només amb --sense-cobertes)")
    ap.add_argument("--extra-hores", type=float, default=0, help="hores addicionals sense pluja després de l'última dada de radar, per veure com es buida l'aigua")
    ap.add_argument("--cada", type=int, default=10, help="minuts entre fotogrames")
    ap.add_argument("--sense-cobertes", action="store_true", help="sense el mapa de cobertes de l'ICGC: WorldCover i --n-natural/--n-urba únics")
    ap.add_argument("--grup-sol", help="grup hidrològic de sòl per al número de corba: A-D (tota la conca igual), «mapa» "
                    "(scripts/sols.py, SoilGrids NEH-630 0-60 cm + roca), «mapa-sup» (SoilGrids 0-30 cm) o «litologia» (permeabilitat IGME); per defecte el de config/cobertes_parametres.json")
    ap.add_argument("--sol-referencia", choices=["mitjana", "B"], default="mitjana",
                    help="amb el mapa de sòls: «mitjana» conserva la mitjana calibrada (només reparteix); «B» pren la calibració "
                    "com feta amb grup B a tot arreu, i el mapa també en canvia la mitjana")
    ap.add_argument("--planta-baixa", nargs="+", help="OPCIÓ: ús de la planta baixa (planta_baixa.py: data/catalunya/planta_baixa/cNNN.json, una o més conques); amb --edificis, compta les plantes baixes amb aigua a tocar per ús")
    ap.add_argument("--passos", choices=["tapats", "oberts", "arrossegalls"], default="tapats",
                    help="passos de torrentets sota vies (scripts/passos_2d.py): tapats (per defecte, com l'MDT), oberts amb la seva capacitat, "
                         "o tapats en part segons el risc d'arrossegalls (ARR 2019), també als ponts")
    ap.add_argument("--obertura", choices=["petita", "mitjana"], default="petita",
                    help="amb --passos arrossegalls, obertura dels passos de mida desconeguda (taula de bloqueig ARR)")
    ap.add_argument("--edificis", help="edificis del Cadastre (cadastre.py: data/<conca>/cadastre/edificis.json.gz): obstacles i soterranis")
    ap.add_argument("--edificis-diposit", action="store_true", help="OPCIÓ: els edificis són dipòsits; l'aigua hi entra per portals i "
                    "rampes i omple soterranis i plantes baixes (diposit.py). Cal --edificis")
    ap.add_argument("--diposits", help="dipòsits precalculats (diposit.py: diposits.json.gz o Release diposits-cat); sense, es calculen de --edificis")
    ap.add_argument("--velocitat", action="store_true", help="OPCIÓ: calcula i desa la velocitat màxima de l'aigua i el perill AIDR H1-H6")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    rad = json.loads(Path(a.radar).read_text())
    ep = json.loads(Path(a.episodi).read_text()) if a.episodi else None

    # relleu agregat a --pas m
    with rasterio.open(a.mdt) as f:
        fac = max(int(round(a.pas / f.res[0])), 1)
        z0 = f.read(1)
        nf, nc = f.height // fac, f.width // fac
        blocs = z0[: nf * fac, : nc * fac].reshape(nf, fac, nc, fac)
        with np.errstate(all="ignore"):
            z = np.nanmin(np.nanmin(blocs, axis=3), axis=1)            # el mínim de cada bloc conserva els llits
        tr = f.transform * f.transform.scale(fac, fac)
        tr_mdt = f.transform
        crs = f.crs
    mar = ~np.isfinite(z)
    z = np.where(mar, 0.0, z).astype(np.float64)
    dx = float(tr.a)
    print(f"Domini: {nc}×{nf} cel·les de {dx:.0f} m")

    # cobertes del sòl: mapa de l'ICGC (Release cobertes-cat, scripts/cobertes.py) → n de Manning per coberta (Chow 1959)
    # i absorció repartida segons la retenció S del número de corba (TR-55); sense, ESA WorldCover i valors únics.
    cob = None
    if not a.sense_cobertes:
        sys.path.insert(0, str(Path(__file__).parent))
        import cobertes as C
        cob = C.a_reixa(tr, (nf, nc), crs)
    if cob is not None:
        grup = a.grup_sol
        if grup in ("mapa", "mapa-sup", "litologia"):
            import sols
            grup = sols.a_reixa(tr, (nf, nc), crs, variant=a.grup_sol)
            if grup is None:
                sys.exit("No es pot llegir el mapa de grups hidrològics de sòl (Release sols-cat)")
            v = grup[grup >= 0]
            print("Grups de sòl (" + a.grup_sol + "): " + ", ".join(f"{'ABCD'[k]} {(v == k).mean() * 100:.0f} %" for k in range(4)))
        n, cn, tipus = C.parametres(cob, grup_sol=grup)
        urba = tipus == 1
        aigua = tipus == 2
        S = 25400 / cn - 254                                     # retenció potencial màxima (mm), SCS
        nat = ~urba & ~aigua & ~mar
        Sref = S if a.sol_referencia == "mitjana" else 25400 / C.parametres(cob, grup_sol="B")[1] - 254
        rel = S / Sref[nat].mean() if nat.any() else np.ones_like(S)   # la mitjana de referència conserva el valor calibrat
        ia = np.where(urba, a.ia_urba, a.ia_natural * rel)
        sorti = np.where(urba, a.clavegueram, a.perdua_natural * rel) / 1000 / 3600
        ia[aigua] = 0; sorti[aigua] = 0                          # sobre l'aigua la pluja no s'infiltra
        print(f"Cobertes ICGC: construït {urba.mean() * 100:.1f} %, aigua {aigua.mean() * 100:.1f} %, n natural mitjana "
              f"{n[nat].mean():.3f}, absorció natural {a.ia_natural * rel[nat].min():.0f}-{a.ia_natural * rel[nat].max():.0f} mm")
    else:
        urba = np.zeros((nf, nc), dtype=bool)
        for lon0 in (0, 3):
            try:
                with rasterio.open(WORLDCOVER.format(lon0)) as f:
                    dst = np.zeros((nf, nc), dtype=np.uint8)
                    reproject(rasterio.band(f, 1), dst, dst_transform=tr, dst_crs=crs, resampling=Resampling.mode)
                    urba |= dst == 50
            except Exception as e:
                print(f"WorldCover E{lon0:03d}: {type(e).__name__} {e}")
        print(f"Construït (WorldCover): {urba.mean() * 100:.1f} % del domini")
        ia = np.where(urba, a.ia_urba, a.ia_natural)
        sorti = np.where(urba, a.clavegueram, a.perdua_natural) / 1000 / 3600     # m/s que perd la làmina un cop superat el llindar
        n = np.where(urba, a.n_urba, a.n_natural)
    # clavegueram estimat per subconca (opció): capacitat per cel·la i cel·la de desguàs
    clav = None
    if a.escenari_clavegueram != "fix":
        sys.path.insert(0, str(Path(__file__).parent))
        import clavegueram_cat as CL
        r = CL.a_reixa(tr, (nf, nc), crs, urba)
        if r is None:
            print("Clavegueram estimat: cap conca precalculada al domini; es fa servir el fix")
        else:
            c10, c2, desg = r
            cap = {"funciona": c10, "antic": c2, "collapsat": c10 * a.factor_collapsat}[a.escenari_clavegueram]
            amb = urba & np.isfinite(cap)
            sorti = np.where(amb, np.nan_to_num(cap) / 1000 / 3600, sorti)
            ok = amb & (desg >= 0) & (desg != np.arange(nf * nc).reshape(nf, nc))
            clav = {"cel": np.flatnonzero(ok), "desg": desg[ok]}
            clav["ofeg"] = np.flatnonzero(amb)
            clav["ofeg_desg"] = desg[amb]
            print(f"Clavegueram «{a.escenari_clavegueram}»: {np.nanmean(cap[amb]):.0f} mm/h de mitjana a les zones urbanes, "
                  f"{len(clav['cel'])} cel·les desguassen a dins el domini")
    phi, psx, psy, desti, edif = obstacles(a.edificis, tr, nf, nc, dx) if a.edificis else (np.ones((nf, nc)), np.ones((nf, nc)), np.ones((nf, nc)), None, None)
    dip = None
    if a.edificis_diposit:
        if not a.edificis:
            sys.exit("--edificis-diposit necessita --edificis")
        sys.path.insert(0, str(Path(__file__).parent))
        import diposit as D
        if a.diposits:
            files = D.llegeix(a.diposits)
        else:
            files = D.per_edifici(edif["parts"])
        dip = D.Diposits(files, tr, nf, nc, phi, desti, z0, tr_mdt)
    t_pas = None
    if a.passos != "tapats":
        sys.path.insert(0, str(Path(__file__).parent))
        import passos_2d as P2
        t_pas, info_pas, resum_pas = P2.prepara(a.passos, a.obertura, tr, nf, nc, crs, z, mar, psx, psy)
        qmax_pas, vol_pas = np.zeros(len(t_pas)), np.zeros(len(t_pas))
    # clavegueram: l'aigua que entra pels embornals surt a la riera més propera (data/<conca>/rieretes.geojson, relleu.py)
    desti_clav = None
    fr = Path("data") / cfg["conca"] / "rieretes.geojson"
    if a.clavegueram_desti == "riu" and fr.exists():
        from rasterio import features as rfeat
        from scipy import ndimage
        from shapely.geometry import shape as shp
        from shapely.ops import transform as stransf
        to_c = Transformer.from_crs(4326, crs, always_xy=True)
        trams = [stransf(to_c.transform, shp(f["geometry"])) for f in json.loads(fr.read_text())["features"]
                 if f["properties"].get("area_km2", 0) >= a.riu_km2]
        if trams:
            riu = rfeat.rasterize([(g, 1) for g in trams], out_shape=(nf, nc), transform=tr, all_touched=True).astype(bool) & ~mar
            _, (ri_, ci_) = ndimage.distance_transform_edt(~riu, return_indices=True)
            desti_clav = (ri_ * nc + ci_).ravel()
            print(f"Clavegueram: surt a {len(trams)} trams de riera amb conca ≥ {a.riu_km2} km²")

    # pluja: graella del radar a la reixa, amb correcció horària amb el pluviòmetre
    to_ll = Transformer.from_crs(crs, 4326, always_xy=True)
    xs = tr.c + (np.arange(nc) + 0.5) * dx
    ys = tr.f - (np.arange(nf) + 0.5) * dx
    X, Y = np.meshgrid(xs, ys)
    LON, LAT = to_ll.transform(X, Y)
    lons, lats = np.array(rad["lons"]), np.array(rad["lats"])
    ci = np.clip(np.rint((LON - lons[0]) / rad["pas_graus"]).astype(int), 0, len(lons) - 1)
    ri = np.clip(np.rint((LAT - lats[0]) / rad["pas_graus"]).astype(int), 0, len(lats) - 1)
    idx = ri * len(lons) + ci
    frames = np.array(rad["mm"])                                  # (temps, cel·les radar) mm en interval
    temps_rad = [datetime.fromisoformat(t) for t in rad["temps"]]
    # correcció: mm del pluviòmetre / mm del radar a la conca de l'aforador, per hores
    t_ep = [datetime.fromisoformat(t) for t in ep["temps_utc"]] if ep else []
    plv = np.array(ep["pluja_pluviometre_mm"], dtype=float) if ep else None
    rdc = np.array(ep["pluja_radar_conca_mm"], dtype=float) if ep else None
    def factor(t):
        if not ep:
            return 1.0
        sel = [k for k, te in enumerate(t_ep) if t - timedelta(minutes=60) < te <= t]
        p, r = np.nansum(plv[sel]), np.nansum(rdc[sel])
        return float(np.clip(p / r, 0.3, 4)) if r > 0.5 else 1.0

    # punt de control: transsecte N-S d'uns 300 m a l'aforador; s'hi mesura el cabal cap a l'est
    to_utm = Transformer.from_crs(4326, crs, always_xy=True)
    gx, gy = to_utm.transform(cfg["aforador"]["lon"], cfg["aforador"]["lat"])
    gc, gr = int((gx - tr.c) // dx), int((tr.f - gy) // dx)
    tf = range(max(gr - int(150 // dx), 0), min(gr + int(150 // dx) + 1, nf))

    h = np.zeros((nf, nc)); qx = np.zeros((nf, nc)); qy = np.zeros((nf, nc))
    actiu = ~mar
    actiu[0, :] = actiu[-1, :] = actiu[:, 0] = actiu[:, -1] = False   # vores del domini obertes: l'aigua surt
    pluja_acum = np.zeros((nf, nc)); neta_ant = np.zeros((nf, nc))
    hmax = np.zeros((nf, nc))
    t15 = np.full((nf, nc), np.inf)                               # primer moment amb més de 15 cm d'aigua (min des de l'inici)
    t_mullat = np.zeros((nf, nc))                                 # darrer moment amb més de 5 cm d'aigua (min des de l'inici)
    if a.velocitat:
        vmax = np.zeros((nf, nc)); hvmax = np.zeros((nf, nc)); perill = np.zeros((nf, nc), dtype=np.uint8)
    ini = temps_rad[0] - timedelta(minutes=rad["interval_min"])
    fi_radar = temps_rad[-1]; fi = fi_radar + timedelta(hours=a.extra_hores)
    dest = Path("data") / cfg["conca"] / "sim2d" / (a.id or fi.strftime("%Y-%m-%d")); dest.mkdir(parents=True, exist_ok=True)
    # sortida en lat/lon per al visor
    w, s, e, nn = rasterio.warp.transform_bounds(crs, "EPSG:4326", tr.c, tr.f - nf * dx, tr.c + nc * dx, tr.f)
    pas_g = dx / 111320
    ll_tr = from_origin(w, nn, pas_g, pas_g)
    lnc, lnf = int((e - w) / pas_g), int((nn - s) / pas_g)
    def png(hh, nom):
        out = np.zeros((lnf, lnc), dtype=np.float32)
        reproject(hh.astype(np.float32), out, src_transform=tr, src_crs=crs, dst_transform=ll_tr, dst_crs="EPSG:4326",
                  resampling=Resampling.nearest)
        rgba = np.zeros((lnf, lnc, 4), dtype=np.uint8)
        for llindar, c in ((0.05, (189, 215, 238, 120)), (0.15, (107, 174, 214, 170)), (0.3, (33, 113, 181, 200)),
                           (0.6, (8, 69, 148, 220)), (1.2, (63, 0, 125, 235))):
            rgba[out >= llindar] = c
        Image.fromarray(rgba, "RGBA").save(dest / nom, optimize=True)

    t, k_frame, k_taxa, taxa, fotos, hidro = ini, 0, -1, None, [], []
    hidro_caixa, hidro_trans = [], []
    bal = {"pluja": 0.0, "neta": 0.0, "infiltrat": 0.0, "clavegueram": 0.0}   # m³: balanç d'aigua de tot el domini
    t_seg, dt_max, rellotge = 0.0, 5.0, time.time()
    proxim = ini + timedelta(minutes=a.cada)
    while t < fi:
        # pluja de l'interval de radar actual, repartida uniformement
        while k_frame < len(temps_rad) - 1 and temps_rad[k_frame] < t + timedelta(seconds=1):
            k_frame += 1
        if t >= fi_radar:                                              # passat l'últim fotograma de radar no plou més
            if k_taxa != -2:
                taxa = np.zeros_like(taxa); k_taxa = -2
        elif k_frame != k_taxa:                                       # la pluja només canvia amb cada fotograma de radar
            taxa = frames[k_frame][idx] * factor(temps_rad[k_frame]) / 1000 / (rad["interval_min"] * 60)   # m/s
            if desti is not None:                                     # la pluja de les teulades sense carrer a la cel·la baixa al carrer més proper
                taxa = np.bincount(desti, weights=taxa.ravel(), minlength=nf * nc).reshape(nf, nc)
            k_taxa = k_frame
        hm = hmax_actiu(h, actiu)
        dt = min(dt_max, 0.7 * dx / np.sqrt(G * max(hm, 0.01)))
        # pèrdues: només escorre la pluja per sobre del llindar
        pluja_acum += taxa * dt * 1000
        neta = np.maximum(pluja_acum - ia, 0.0)
        dn = neta - neta_ant
        if a.perdua_sobre == "pluja":                             # el terreny natural infiltra de la pluja que hi cau, fins a --perdua-natural mm/h
            dn = np.where(urba, dn, np.maximum(dn - sorti * dt * 1000, 0.0))
        h += dn / 1000 / phi                                      # la pluja de tota la cel·la (teulades incloses) va al carrer
        bal["pluja"] += float((taxa * dt)[~mar].sum()) * dx * dx; bal["neta"] += float(dn[~mar].sum()) / 1000 * dx * dx
        h_abans = h.copy()
        neta_ant = neta
        # al terreny natural, les lleres on ja corre aigua (> --h-llera m) tenen el fons saturat i no hi infiltren; el clavegueram sí que drena
        perd_m = (pluja_acum >= ia) & (urba | ((h < a.h_llera) & (a.perdua_sobre == "calat")))
        perd = np.minimum(h, np.where(perd_m, sorti * dt / np.where(urba, phi, 1.0), 0.0))
        if clav is not None and a.clavegueram_ofegat > 0:         # desguàs sota l'aigua: la subconca no buida
            hd = np.where(clav["ofeg_desg"] >= 0, h.ravel()[np.maximum(clav["ofeg_desg"], 0)], 0.0)
            perd.ravel()[clav["ofeg"][hd > a.clavegueram_ofegat]] = 0.0
        h -= perd
        dperd = (h_abans - h) * phi
        bal["infiltrat"] += float(dperd[~urba].sum()) * dx * dx; bal["clavegueram"] += float(dperd[urba].sum()) * dx * dx
        if clav is not None:                                      # clavegueram estimat: desguassa a la cel·la de sortida de la subconca
            v = perd.ravel()[clav["cel"]] * phi.ravel()[clav["cel"]]      # m³/m² de cel·la que entra a la xarxa
            h.ravel()[:] += np.bincount(clav["desg"], weights=v, minlength=nf * nc) / phi.ravel()
        elif desti_clav is not None:                              # el clavegueram desguassa a la riera (volum per cel·la → calat)
            h += np.bincount(desti_clav, weights=np.where(urba, dperd, 0.0).ravel(), minlength=nf * nc).reshape(nf, nc) / phi
        if dip is not None:                                       # entrada i sortida d'aigua dels edificis
            dip.pas(z, h, phi, dx, dt, t_seg / 60)
        pas_temps(z, h, qx, qy, n, dx, dt, actiu, psx, psy, phi)
        if t_pas is not None and len(t_pas):                      # conductes dels passos sota vies
            P2.transfereix(t_pas, z, h, phi, dx, dt, qmax_pas, vol_pas)
        if a.velocitat:
            velocitat(h, qx, qy, actiu, vmax, hvmax, perill, 0.05)
        np.maximum(hmax, h, out=hmax)
        t15[(h > 0.15) & ~np.isfinite(t15)] = t_seg / 60
        mullat = h > 0.05
        t_mullat[mullat] = t_seg / 60                              # darrer minut amb més de 5 cm d'aigua
        t_seg += dt; t = ini + timedelta(seconds=t_seg)
        if t >= proxim:
            nom = f"{proxim:%d%H%M}.png"; png(h, nom)
            # cabal a l'aforador de dues maneres: el transsecte N-S (cap a l'est) i el que surt d'una caixa de ±150 m al voltant
            # del punt, sigui quina sigui la direcció del curs (si el riu fa un meandre pel transsecte, aquest descompta)
            k = max(int(150 // dx), 1); r0, r1, c0, c1 = max(gr - k, 1), min(gr + k, nf - 2), max(gc - k, 1), min(gc + k, nc - 2)
            Qx, Qy = qx * psx, qy * psy
            q_caixa = float((np.maximum(Qx[r0:r1 + 1, c1], 0).sum() + np.maximum(-Qx[r0:r1 + 1, c0 - 1], 0).sum() +
                             np.maximum(Qy[r1, c0:c1 + 1], 0).sum() + np.maximum(-Qy[r0 - 1, c0:c1 + 1], 0).sum()) * dx)
            q_trans = float(sum(qx[r, gc] * psx[r, gc] for r in tf) * dx)
            q = q_caixa if cfg.get("aforador_mode") == "caixa" else q_trans
            hidro_caixa.append(round(q_caixa, 1)); hidro_trans.append(round(q_trans, 1))
            fotos.append({"temps": proxim.strftime("%Y-%m-%dT%H:%M"), "png": nom}); hidro.append(round(q, 1))
            print(f"{proxim:%d %H:%M} dt {dt:.2f} s, calat màx {hm:.2f} m, cabal aforador {q:.1f} m³/s, inundat >15 cm {(h > 0.15).sum() * dx * dx / 1e6:.2f} km² ({time.time() - rellotge:.0f} s)")
            proxim += timedelta(minutes=a.cada)
    png(hmax, "maxim.png")
    # temps d'assecat: des que supera 15 cm fins l'últim minut amb més de 5 cm; les cel·les encara mullades al final queden marcades
    mullat_final = h > 0.05
    durada = np.where(np.isfinite(t15), np.maximum(t_mullat - np.nan_to_num(t15), 0), np.nan)
    d_h = durada / 60
    def png_assecat(nom):
        out = np.full((lnf, lnc), np.nan, dtype=np.float32)
        reproject(d_h.astype(np.float32), out, src_transform=tr, src_crs=crs, dst_transform=ll_tr, dst_crs="EPSG:4326",
                  resampling=Resampling.nearest, src_nodata=np.nan, dst_nodata=np.nan)
        fin = np.zeros((lnf, lnc), dtype=np.float32)
        reproject(mullat_final.astype(np.float32), fin, src_transform=tr, src_crs=crs, dst_transform=ll_tr, dst_crs="EPSG:4326", resampling=Resampling.nearest)
        rgba = np.zeros((lnf, lnc, 4), dtype=np.uint8)
        ok = np.isfinite(out)
        for llindar, c in ((0, (255, 255, 178, 200)), (1, (254, 204, 92, 210)), (3, (253, 141, 60, 220)), (6, (240, 59, 32, 230)), (12, (189, 0, 38, 240))):
            rgba[ok & (out >= llindar)] = c
        rgba[ok & (fin > 0.5)] = (84, 39, 143, 240)          # encara mullat al final de la simulació
        Image.fromarray(rgba, "RGBA").save(dest / nom, optimize=True)
    png_assecat("assecat.png")
    np.where(np.isfinite(durada), np.minimum(durada, 65534), 65535).astype("<u2").tofile(dest / "durada_min.u2")
    # graelles numèriques per al model d'inundació (UTM, mateixa reixa que la simulació): calat màxim en cm i minut d'inundació
    sot = None
    if edif is not None:
        n_sot, sot = soterranis(edif, tr, nf, nc, dx, hmax, t15, durada, z, z0, tr_mdt)
        (dest / "soterranis.json").write_text(json.dumps({"inici_utc": ini.strftime("%Y-%m-%dT%H:%M"), "toca_m": 5, "llindar_cm": 15,
            "camps": ["aigua a tocar per sobre del terra de l'edifici, màxim (cm)", "minut des de inici_utc en què passa de 15 cm (null = mai)", "minuts des que passa de 15 cm fins que l'aigua del voltant baixa de 5 cm"], "edificis": sot}, separators=(",", ":")))
        print(f"Soterranis: {sum(1 for v in sot.values() if v[0] >= 15)} de {n_sot} edificis amb soterrani tenen ≥15 cm d'aigua a tocar")
    pb_resum = None
    if edif is not None and a.planta_baixa:
        pb_resum = plantes_baixes(a.planta_baixa, edif, dest, ini, tr, nf, nc, dx, hmax, t15, durada, z, z0, tr_mdt)
    if t_pas is not None:
        k = 0
        for reg in info_pas:
            if reg.get("estat") == "conducte":
                (ru, cu), (rd, cd) = reg["entrada"], reg["sortida"]
                reg.update({"q_max_m3s": round(float(qmax_pas[k]), 2), "volum_m3": round(float(vol_pas[k])),
                            "calat_max_entrada_m": round(float(hmax[ru, cu]), 2), "calat_max_sortida_m": round(float(hmax[rd, cd]), 2),
                            "vessa": bool(hmax[ru, cu] >= reg["fondaria_m"])})
                k += 1
        (dest / "passos.json").write_text(json.dumps({"mode": a.passos, "obertura": a.obertura if a.passos == "arrossegalls" else None,
                                                      "resum": resum_pas, "passos": info_pas}, ensure_ascii=False, separators=(",", ":")))
    np.minimum(np.rint(hmax * 100), 65534).astype("<u2").tofile(dest / "calat_max_cm.u2")
    fitxers_ext, vel_info, dip_info = {}, None, None
    if a.velocitat:
        def png_classes(v, llindars, colors, nom):
            out = np.zeros((lnf, lnc), dtype=np.float32)
            reproject(v.astype(np.float32), out, src_transform=tr, src_crs=crs, dst_transform=ll_tr, dst_crs="EPSG:4326", resampling=Resampling.nearest)
            rgba = np.zeros((lnf, lnc, 4), dtype=np.uint8)
            for ll, c in zip(llindars, colors):
                rgba[out >= ll] = c
            Image.fromarray(rgba, "RGBA").save(dest / nom, optimize=True)
        png_classes(vmax, [0.25, 0.5, 1, 2, 3], [(255, 255, 204, 170), (161, 218, 180, 190), (65, 182, 196, 210), (34, 94, 168, 225), (37, 52, 148, 240)], "velocitat.png")
        png_classes(perill, [1, 2, 3, 4, 5, 6], [(255, 255, 178, 150), (254, 217, 118, 180), (254, 178, 76, 200), (253, 141, 60, 215), (240, 59, 32, 230), (177, 0, 38, 240)], "perill.png")
        np.minimum(np.rint(vmax * 100), 65534).astype("<u2").tofile(dest / "velocitat_max_cms.u2")
        np.minimum(np.rint(hvmax * 100), 65534).astype("<u2").tofile(dest / "hv_max_cm2s.u2")
        perill.astype("u1").tofile(dest / "perill_max.u1")
        fitxers_ext = {"velocitat_max_cms.u2": "uint16 velocitat màxima (cm/s) on el calat passa de 5 cm",
                       "hv_max_cm2s.u2": "uint16 màxim de calat × velocitat (m²/s × 100)",
                       "perill_max.u1": "uint8 classe de perill AIDR H1-H6 màxima (0 = sec)"}
        moll = hmax > 0.05
        vel_info = {"png": "velocitat.png", "perill_png": "perill.png", "llegenda_velocitat_ms": [0.25, 0.5, 1, 2, 3],
                    "perill": "AIDR Handbook 7 (Smith et al. 2014): H1 ≤0,3 m i hv ≤0,3; H2 ≤0,5 m, hv ≤0,6 (cotxes petits); H3 ≤1,2 m (vehicles, infants); H4 ≤2 m, hv ≤1 (persones); H5 ≤4 m, v ≤4, hv ≤4 (edificis vulnerables); H6 la resta",
                    "km2_per_classe": {f"H{c}": round(float((perill == c).sum()) * dx * dx / 1e6, 3) for c in range(1, 7)},
                    "velocitat_p99_ms": round(float(np.percentile(vmax[moll], 99)), 2) if moll.any() else 0}
        print(f"Velocitat: p99 {vel_info['velocitat_p99_ms']} m/s; perill km² {vel_info['km2_per_classe']}")
    if dip is not None:
        res = dip.resultat()
        (dest / "edificis_diposit.json").write_text(json.dumps({"inici_utc": ini.strftime("%Y-%m-%dT%H:%M"),
            "camps": ["calat màxim a la planta baixa (cm)", "% del soterrani ple (null = no en té)", "minut des de inici_utc en què hi comença a entrar aigua", "m³ màxims dins l'edifici"],
            "edificis": res}, separators=(",", ":")))
        dip_info = {"fitxer": "edificis_diposit.json", "hipotesis": D.hipotesis(), "edificis_amb_aigua": len(res),
                    "soterranis_plens_50pct": sum(1 for v in res.values() if v[1] is not None and v[1] >= 50),
                    "plantes_baixes_15cm": sum(1 for v in res.values() if v[0] >= 15),
                    "entrat_hm3": round(dip.entrat / 1e6, 4), "tornat_al_carrer_hm3": round(dip.sortit / 1e6, 4),
                    "retingut_final_hm3": round(float(dip.V.sum()) / 1e6, 4)}
        print(f"Dipòsit: aigua a {len(res)} edificis; entra {dip_info['entrat_hm3']} hm³, en torna {dip_info['tornat_al_carrer_hm3']}, "
              f"queda {dip_info['retingut_final_hm3']}; {dip_info['soterranis_plens_50pct']} soterranis ≥50 % plens, {dip_info['plantes_baixes_15cm']} plantes baixes ≥15 cm")
    np.where(np.isfinite(t15), np.minimum(t15, 65534), 65535).astype("<u2").tofile(dest / "minut_15cm.u2")
    (dest / "index.json").write_text(json.dumps({
        "limits": [[s, w], [nn, e]], "pas_m": dx, "inici_utc": ini.strftime("%Y-%m-%dT%H:%M"),
        "graella_utm": {"crs": str(crs), "x0": tr.c, "y0": tr.f, "nf": nf, "nc": nc, "fitxers": {"calat_max_cm.u2": "uint16 cm",
                        "minut_15cm.u2": "uint16 minuts des de inici_utc fins a superar 15 cm; 65535 = mai", "durada_min.u2": "uint16 minuts entre superar 15 cm i l'últim moment amb més de 5 cm; 65535 = mai superat", **fitxers_ext}}, "fotogrames": fotos, "q_aforador_sim": hidro, "q_aforador_caixa": hidro_caixa, "q_aforador_transsecte": hidro_trans, "maxim": "maxim.png", "assecat": "assecat.png", "assecat_final_mullat_km2": round(float(mullat_final.sum()) * dx * dx / 1e6, 3), "llegenda_assecat_h": [0, 1, 3, 6, 12],
        "hipotesis": {"ia_natural_mm": a.ia_natural, "ia_urba_mm": a.ia_urba, "perdua_natural_mm_h": a.perdua_natural, "perdua_sobre": a.perdua_sobre, "clavegueram_mm_h": a.clavegueram,
                      "clavegueram_desti": "subconca" if clav is not None else "riu" if desti_clav is not None else "perdut", "clavegueram_escenari": a.escenari_clavegueram if clav is not None else "fix",
                      "clavegueram_factor_collapsat": a.factor_collapsat if a.escenari_clavegueram == "collapsat" else None, "clavegueram_ofegat_m": a.clavegueram_ofegat or None, "h_llera_m": a.h_llera,
                      "passos": "tapats (MDT)" if t_pas is None else {"mode": a.passos, "obertura": a.obertura if a.passos == "arrossegalls" else None, **resum_pas, "fitxer": "passos.json"}, "manning": ({"natural": a.n_natural, "construit": a.n_urba} if cob is None else "per coberta (config/cobertes_parametres.json, Chow 1959)"),
                      "cobertes": None if cob is None else {"font": "Mapa de cobertes del sòl de Catalunya, ICGC (CC BY 4.0)", "absorcio": "ia_natural × S/S̄ (número de corba TR-55)", "grup_sol": a.grup_sol or json.loads(Path("config/cobertes_parametres.json").read_text())["grup_hidrologic_sol"],
                                                    "sol_referencia": a.sol_referencia},
                      "pluja": "radar OPERA corregit cada hora amb el pluviòmetre de l'ACA" if ep else "radar OPERA sense corregir", "esquema": "inercial local (Bates et al. 2010)"},
        "llegenda_calat_m": [0.05, 0.15, 0.3, 0.6, 1.2],
        "opcions": {"edificis_diposit": dip is not None, "velocitat": bool(a.velocitat)},
        "velocitat": vel_info, "diposit": dip_info,
        "balanc_hm3": {**{k: round(v / 1e6, 3) for k, v in bal.items()}, "retingut_llindar": round((bal["pluja"] - bal["neta"]) / 1e6, 3),
                       "queda_al_final": round(float((h * phi)[~mar].sum()) * dx * dx / 1e6, 3),
                       "queda_en_clots_gt30cm": round(float((h * phi)[(~mar) & (h > 0.3)].sum()) * dx * dx / 1e6, 3),
                       "sortit_pel_mar_i_vores": round((bal["neta"] - bal["infiltrat"] - (0 if (desti_clav is not None or clav is not None) else bal["clavegueram"]) - float((h * phi)[~mar].sum()) * dx * dx) / 1e6, 3),
                       "volum_aforador_caixa": round(sum(hidro_caixa) * a.cada * 60 / 1e6, 3)},
        "edificis": None if edif is None else {"font": "Cadastre (INSPIRE Buildings)", "metode": "porositat: els edificis treuen espai i amplada de pas a cada cel·la; la pluja de les teulades va al carrer",
                                               "ocupacio_pct": edif["ocupacio_pct"], "soterranis": "soterranis.json",
                                               "soterranis_risc": sum(1 for v in sot.values() if v[0] >= 15)},
        "plantes_baixes": pb_resum}, ensure_ascii=False))
    print(f"Fet en {time.time() - rellotge:.0f} s: {len(fotos)} fotogrames")
    return 0


if __name__ == "__main__":
    sys.exit(main())
