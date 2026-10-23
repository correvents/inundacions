#!/usr/bin/env python3
"""Simulació 2D de la pluja sobre el relleu (rain-on-grid): on corre i on s'embassa l'aigua a cada moment.

Esquema inercial local de Bates et al. (2010), el de LISFLOOD-FP: a cada cara entre cel·les el cabal unitari q s'actualitza
amb el pendent de la superfície de l'aigua i la fricció de Manning, i la làmina h amb la divergència de q.
  - Relleu: la reixa de relleu.py (MDT 5 m de l'IGN) agregada a --pas metres.
  - Pluja: graella del radar OPERA de l'episodi (15 min, 1 km), corregida hora a hora perquè la mitjana sobre la conca
    de l'aforador coincideixi amb el pluviòmetre de l'ACA (el radar infravalora els nuclis convectius).
  - Pèrdues: llindar per cel·la. Terreny natural: --ia-natural mm (calibrat amb l'aforador el 3-4/10/2026).
    Zones construïdes (ESA WorldCover 2021, classe 50): --ia-urba mm.
    Després del llindar, la làmina d'aigua de cada cel·la perd --perdua-natural mm/h (infiltració) al terreny natural i
    --clavegueram mm/h a les zones construïdes (embornals), excepte on el calat supera --h-llera (llera saturada). Així l'aigua embassada també es buida quan para de ploure.
  - Fricció: n de Manning (--n-natural, --n-urba).
  - Mar i vores del domini: l'aigua surt.
Sortides (data/<conca>/sim2d/<id>/): un PNG de calat cada --cada minuts en lat/lon, calat màxim, l'hidrograma simulat a
l'aforador (per validar-lo amb l'observat), index.json i, per al model d'inundació, dues graelles numèriques en UTM:
calat màxim (cm) i minut en què cada cel·la supera 15 cm d'aigua.

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
def pas_temps(z, h, qx, qy, n, dx, dt, actiu, psx, psy, phi, limx, limy):
    """Un pas de l'esquema inercial local. qx a la cara est de cada cel·la, qy a la cara sud.
    Porositat dels edificis: q és el cabal per metre d'amplada lliure; psx/psy, la fracció lliure de cada cara, i phi, la
    fracció lliure de cada cel·la (h és el calat a la part lliure, el carrer). Sense edificis tot val 1.
    limx/limy: cabal màxim per metre a cada cara (obres de drenatge sota terraplens; infinit a la resta)."""
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
            qmax = min(hf * np.sqrt(G * hf), limx[i, j])   # limitador de Froude 1 per estabilitat
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
            qmax = min(hf * np.sqrt(G * hf), limy[i, j])
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


def soterranis(ed, tr, nf, nc, dx, hmax, t15, durada, z, z_mdt, tr_mdt, toca=5.0):
    """Edificis amb plantes sota rasant que tenen aigua a tocar: a les cel·les a menys de `toca` m de la planta de les parts amb
    soterrani, l'alçada màxima de l'aigua per sobre del terra de l'edifici (mediana del MDT sota la planta), sense passar del
    calat. Així l'aigua d'una riera encaixada al costat, més baixa que l'edifici, no compta. I minut en què hi passa de 15 cm."""
    from shapely.geometry import Polygon, box
    from shapely.ops import unary_union
    per_ref = {}
    for ref, sobre, sota, an in ed["parts"]:
        if sota > 0 and len(an[0]) >= 4:
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
    ap.add_argument("--h-llera", type=float, default=0.1, help="calat a partir del qual la cel·la és llera saturada i no infiltra (m)")
    ap.add_argument("--clavegueram", type=float, default=30, help="capacitat dels embornals a les zones construïdes (mm/h)")
    ap.add_argument("--n-natural", type=float, default=0.06, help="n de Manning del terreny natural")
    ap.add_argument("--n-urba", type=float, default=0.02, help="n de Manning de les zones construïdes")
    ap.add_argument("--extra-hores", type=float, default=0, help="hores addicionals sense pluja després de l'última dada de radar, per veure com es buida l'aigua")
    ap.add_argument("--cada", type=int, default=10, help="minuts entre fotogrames")
    ap.add_argument("--passos-cap", help="capacitat de les obres de drenatge (passos.py: data/<conca>/passos/capacitat.json)")
    ap.add_argument("--obstruccio", type=float, default=0, help="fracció obstruïda de les obres de drenatge (0 = netes)")
    ap.add_argument("--edificis", help="edificis del Cadastre (cadastre.py: data/<conca>/cadastre/edificis.json.gz): obstacles i soterranis")
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

    # cobertura del sòl (ESA WorldCover): construït → pèrdues i fricció baixes
    urba = np.zeros((nf, nc), dtype=bool)
    for lon0 in (0, 3):
        try:
            with rasterio.open(WORLDCOVER.format(lon0)) as f:
                dst = np.zeros((nf, nc), dtype=np.uint8)
                reproject(rasterio.band(f, 1), dst, dst_transform=tr, dst_crs=crs, resampling=Resampling.mode)
                urba |= dst == 50
        except Exception as e:
            print(f"WorldCover E{lon0:03d}: {type(e).__name__} {e}")
    print(f"Construït: {urba.mean() * 100:.1f} % del domini")
    ia = np.where(urba, a.ia_urba, a.ia_natural)
    sorti = np.where(urba, a.clavegueram, a.perdua_natural) / 1000 / 3600     # m/s que perd la làmina un cop superat el llindar
    n = np.where(urba, a.n_urba, a.n_natural)
    phi, psx, psy, desti, edif = obstacles(a.edificis, tr, nf, nc, dx) if a.edificis else (np.ones((nf, nc)), np.ones((nf, nc)), np.ones((nf, nc)), None, None)

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
    # obres de drenatge: per cada cara d'una cel·la de l'obra no hi passa més del cabal de disseny (× la part no obstruïda)
    limx = np.full((nf, nc), np.inf); limy = np.full((nf, nc), np.inf); n_obres = 0
    if a.passos_cap and Path(a.passos_cap).exists():
        for o in json.loads(Path(a.passos_cap).read_text())["obres"]:
            q = o["q_m3s"] * (1 - a.obstruccio) / dx                  # m²/s per metre d'amplada
            cels = {(int((tr.f - y) // dx), int((x - tr.c) // dx)) for x, y in o["xy"]}
            for i, j in cels:
                if not (0 <= i < nf and 0 <= j < nc):
                    continue
                for ii, jj, lim, ps in ((i, j, limx, psx), (i, j - 1, limx, psx), (i, j, limy, psy), (i - 1, j, limy, psy)):
                    if 0 <= ii < nf and 0 <= jj < nc:
                        lim[ii, jj] = min(lim[ii, jj], q / max(ps[ii, jj], 1e-3))
            n_obres += 1
        print(f"Obres de drenatge limitades al cabal de disseny: {n_obres} (obstrucció {a.obstruccio * 100:.0f} %)")
    actiu[0, :] = actiu[-1, :] = actiu[:, 0] = actiu[:, -1] = False   # vores del domini obertes: l'aigua surt
    pluja_acum = np.zeros((nf, nc)); neta_ant = np.zeros((nf, nc))
    hmax = np.zeros((nf, nc))
    t15 = np.full((nf, nc), np.inf)                               # primer moment amb més de 15 cm d'aigua (min des de l'inici)
    t_mullat = np.zeros((nf, nc))                                 # darrer moment amb més de 5 cm d'aigua (min des de l'inici)
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
        h += (neta - neta_ant) / 1000 / phi                       # la pluja de tota la cel·la (teulades incloses) va al carrer
        neta_ant = neta
        # al terreny natural, les lleres on ja corre aigua (> --h-llera m) tenen el fons saturat i no hi infiltren; el clavegueram sí que drena
        np.maximum(h - np.where((pluja_acum >= ia) & (urba | (h < a.h_llera)), sorti * dt / np.where(urba, phi, 1.0), 0.0), 0.0, out=h)
        pas_temps(z, h, qx, qy, n, dx, dt, actiu, psx, psy, phi, limx, limy)
        np.maximum(hmax, h, out=hmax)
        t15[(h > 0.15) & ~np.isfinite(t15)] = t_seg / 60
        mullat = h > 0.05
        t_mullat[mullat] = t_seg / 60                              # darrer minut amb més de 5 cm d'aigua
        t_seg += dt; t = ini + timedelta(seconds=t_seg)
        if t >= proxim:
            nom = f"{proxim:%d%H%M}.png"; png(h, nom)
            if cfg.get("aforador_mode") == "caixa":       # cabal que surt d'una caixa al voltant del punt, sigui quina sigui la direcció del curs
                k = max(int(150 // dx), 1); r0, r1, c0, c1 = max(gr - k, 1), min(gr + k, nf - 2), max(gc - k, 1), min(gc + k, nc - 2)
                Qx, Qy = qx * psx, qy * psy
                q = float((np.maximum(Qx[r0:r1 + 1, c1], 0).sum() + np.maximum(-Qx[r0:r1 + 1, c0 - 1], 0).sum() +
                           np.maximum(Qy[r1, c0:c1 + 1], 0).sum() + np.maximum(-Qy[r0 - 1, c0:c1 + 1], 0).sum()) * dx)
            else:
                q = float(sum(qx[r, gc] * psx[r, gc] for r in tf) * dx)
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
    np.minimum(np.rint(hmax * 100), 65534).astype("<u2").tofile(dest / "calat_max_cm.u2")
    np.where(np.isfinite(t15), np.minimum(t15, 65534), 65535).astype("<u2").tofile(dest / "minut_15cm.u2")
    (dest / "index.json").write_text(json.dumps({
        "limits": [[s, w], [nn, e]], "pas_m": dx, "inici_utc": ini.strftime("%Y-%m-%dT%H:%M"),
        "graella_utm": {"crs": str(crs), "x0": tr.c, "y0": tr.f, "nf": nf, "nc": nc, "fitxers": {"calat_max_cm.u2": "uint16 cm",
                        "minut_15cm.u2": "uint16 minuts des de inici_utc fins a superar 15 cm; 65535 = mai", "durada_min.u2": "uint16 minuts entre superar 15 cm i l'últim moment amb més de 5 cm; 65535 = mai superat"}}, "fotogrames": fotos, "q_aforador_sim": hidro, "maxim": "maxim.png", "assecat": "assecat.png", "assecat_final_mullat_km2": round(float(mullat_final.sum()) * dx * dx / 1e6, 3), "llegenda_assecat_h": [0, 1, 3, 6, 12],
        "hipotesis": {"ia_natural_mm": a.ia_natural, "ia_urba_mm": a.ia_urba, "perdua_natural_mm_h": a.perdua_natural, "clavegueram_mm_h": a.clavegueram, "obres_drenatge": n_obres, "obstruccio_obres": a.obstruccio, "h_llera_m": a.h_llera, "manning": {"natural": a.n_natural, "construit": a.n_urba},
                      "pluja": "radar OPERA corregit cada hora amb el pluviòmetre de l'ACA" if ep else "radar OPERA sense corregir", "esquema": "inercial local (Bates et al. 2010)"},
        "llegenda_calat_m": [0.05, 0.15, 0.3, 0.6, 1.2],
        "edificis": None if edif is None else {"font": "Cadastre (INSPIRE Buildings)", "metode": "porositat: els edificis treuen espai i amplada de pas a cada cel·la; la pluja de les teulades va al carrer",
                                               "ocupacio_pct": edif["ocupacio_pct"], "soterranis": "soterranis.json",
                                               "soterranis_risc": sum(1 for v in sot.values() if v[0] >= 15)}}, ensure_ascii=False))
    print(f"Fet en {time.time() - rellotge:.0f} s: {len(fotos)} fotogrames")
    return 0


if __name__ == "__main__":
    sys.exit(main())
