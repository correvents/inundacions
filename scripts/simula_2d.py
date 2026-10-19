#!/usr/bin/env python3
"""Simulació 2D de la pluja sobre el relleu (rain-on-grid): on corre i on s'embassa l'aigua a cada moment.

Esquema inercial local de Bates et al. (2010), el de LISFLOOD-FP: a cada cara entre cel·les el cabal unitari q s'actualitza
amb el pendent de la superfície de l'aigua i la fricció de Manning, i la làmina h amb la divergència de q.
  - Relleu: la reixa de relleu.py (MDT 5 m de l'IGN) agregada a --pas metres.
  - Pluja: graella del radar OPERA de l'episodi (15 min, 1 km), corregida hora a hora perquè la mitjana sobre la conca
    de l'aforador coincideixi amb el pluviòmetre de l'ACA (el radar infravalora els nuclis convectius).
  - Pèrdues: llindar per cel·la. Terreny natural: --ia-natural mm (calibrat amb l'aforador el 3-4/10/2026).
    Zones construïdes (ESA WorldCover 2021, classe 50): --ia-urba mm.
  - Fricció: n de Manning (--n-natural, --n-urba).
  - Mar i vores del domini: l'aigua surt.
Sortides (data/<conca>/sim2d/<data>/): un PNG de calat cada --cada minuts en lat/lon, calat màxim, l'hidrograma simulat a
l'aforador (per validar-lo amb l'observat) i index.json.

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
def pas_temps(z, h, qx, qy, n, dx, dt, actiu):
    """Un pas de l'esquema inercial local. qx a la cara est de cada cel·la, qy a la cara sud."""
    nf, nc = z.shape
    for i in range(nf):
        for j in range(nc - 1):
            if not (actiu[i, j] or actiu[i, j + 1]):
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
            if not (actiu[i, j] or actiu[i + 1, j]):
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
            if j > 0: ent += qx[i, j - 1]
            if j < nc - 1: ent -= qx[i, j]
            if i > 0: ent += qy[i - 1, j]
            if i < nf - 1: ent -= qy[i, j]
            h[i, j] = max(h[i, j] + dt * ent / dx, 0.0)


@njit(cache=True)
def hmax_actiu(h, actiu):
    m = 0.0
    for i in range(h.shape[0]):
        for j in range(h.shape[1]):
            if actiu[i, j] and h[i, j] > m:
                m = h[i, j]
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--mdt", required=True)
    ap.add_argument("--radar", required=True)
    ap.add_argument("--episodi", required=True)
    ap.add_argument("--pas", type=float, default=20)
    ap.add_argument("--ia-natural", type=float, default=150)
    ap.add_argument("--ia-urba", type=float, default=5)
    ap.add_argument("--n-natural", type=float, default=0.06, help="n de Manning del terreny natural")
    ap.add_argument("--n-urba", type=float, default=0.02, help="n de Manning de les zones construïdes")
    ap.add_argument("--cada", type=int, default=10, help="minuts entre fotogrames")
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    rad = json.loads(Path(a.radar).read_text())
    ep = json.loads(Path(a.episodi).read_text())

    # relleu agregat a --pas m
    with rasterio.open(a.mdt) as f:
        fac = max(int(round(a.pas / f.res[0])), 1)
        z0 = f.read(1)
        nf, nc = f.height // fac, f.width // fac
        blocs = z0[: nf * fac, : nc * fac].reshape(nf, fac, nc, fac)
        with np.errstate(all="ignore"):
            z = np.nanmin(np.nanmin(blocs, axis=3), axis=1)            # el mínim de cada bloc conserva els llits
        tr = f.transform * f.transform.scale(fac, fac)
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
    n = np.where(urba, a.n_urba, a.n_natural)

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
    t_ep = [datetime.fromisoformat(t) for t in ep["temps_utc"]]
    plv = np.array(ep["pluja_pluviometre_mm"], dtype=float)
    rdc = np.array(ep["pluja_radar_conca_mm"], dtype=float)
    def factor(t):
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
    ini = temps_rad[0] - timedelta(minutes=rad["interval_min"])
    fi = temps_rad[-1]
    dest = Path("data") / cfg["conca"] / "sim2d" / fi.strftime("%Y-%m-%d"); dest.mkdir(parents=True, exist_ok=True)
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
        if k_frame != k_taxa:                                         # la pluja només canvia amb cada fotograma de radar
            taxa = frames[k_frame][idx] * factor(temps_rad[k_frame]) / 1000 / (rad["interval_min"] * 60)   # m/s
            k_taxa = k_frame
        hm = hmax_actiu(h, actiu)
        dt = min(dt_max, 0.7 * dx / np.sqrt(G * max(hm, 0.01)))
        # pèrdues: només escorre la pluja per sobre del llindar
        pluja_acum += taxa * dt * 1000
        neta = np.maximum(pluja_acum - ia, 0.0)
        h += (neta - neta_ant) / 1000
        neta_ant = neta
        pas_temps(z, h, qx, qy, n, dx, dt, actiu)
        np.maximum(hmax, h, out=hmax)
        t_seg += dt; t = ini + timedelta(seconds=t_seg)
        if t >= proxim:
            nom = f"{proxim:%H%M}.png"; png(h, nom)
            q = float(sum(qx[r, gc] for r in tf) * dx)
            fotos.append({"temps": proxim.strftime("%Y-%m-%dT%H:%M"), "png": nom}); hidro.append(round(q, 1))
            print(f"{proxim:%d %H:%M} dt {dt:.2f} s, calat màx {hm:.2f} m, cabal aforador {q:.1f} m³/s, inundat >15 cm {(h > 0.15).sum() * dx * dx / 1e6:.2f} km² ({time.time() - rellotge:.0f} s)")
            proxim += timedelta(minutes=a.cada)
    png(hmax, "maxim.png")
    (dest / "index.json").write_text(json.dumps({
        "limits": [[s, w], [nn, e]], "pas_m": dx, "fotogrames": fotos, "q_aforador_sim": hidro, "maxim": "maxim.png",
        "hipotesis": {"ia_natural_mm": a.ia_natural, "ia_urba_mm": a.ia_urba, "manning": {"natural": a.n_natural, "construit": a.n_urba},
                      "pluja": "radar OPERA corregit cada hora amb el pluviòmetre de l'ACA", "esquema": "inercial local (Bates et al. 2010)"},
        "llegenda_calat_m": [0.05, 0.15, 0.3, 0.6, 1.2]}, ensure_ascii=False))
    print(f"Fet en {time.time() - rellotge:.0f} s: {len(fotos)} fotogrames")
    return 0


if __name__ == "__main__":
    sys.exit(main())
