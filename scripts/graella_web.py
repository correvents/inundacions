#!/usr/bin/env python3
"""Malla de relleu per al simulador de pluja del navegador (pluja.html).

A partir de la reixa de relleu.py (MDT 5 m de l'IGN) prepara, a --pas metres:
  - cota del terreny i cota amb les depressions omplertes (cm, uint16): la diferència és on l'aigua s'embassa,
  - direcció de drenatge D8 sobre el relleu condicionat (pysheds), codi 0–7 (N, NE, E, SE, S, SO, O, NO), 255 = mar o vora.
    Abans s'hi «graven» els rius i rieres d'OpenStreetMap (data/<conca>/passos/osm.json, de passos.py): a cada cel·la de la
    llera la cota passa a ser la més baixa del voltant (3×3) i no puja mai riu avall. Amb cel·les de 20 m el LiDAR no
    marca prou la llera a les planes i, sense això, el cabal fort es desvia per cunetes i camins. Es comprova quina part
    de les cel·les amb més d'1 km² de conca queda a menys de 40 m d'un riu dibuixat, abans i després;
  - zones construïdes (ESA WorldCover 2021, classe 50),
  - índex de la cel·la del radar de l'episodi que cobreix cada cel·la.
Ho desa tot en un binari (data/<conca>/web/graella.bin) amb la descripció a graella.json, i hi afegeix el factor horari
de correcció del radar amb el pluviòmetre per a cada fotograma.

Ús:  python3 scripts/graella_web.py --config config/ridaura.json --mdt /tmp/relleu/mdt.tif \\
         --radar data/ridaura/episodis/2026-10-03_radar.json --episodi data/ridaura/episodis/2026-10-04.json"""
import argparse, json, sys
from datetime import datetime, timedelta
from pathlib import Path
import numpy as np
np.in1d = np.isin
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject, transform_bounds
from rasterio.windows import from_bounds
from pyproj import Transformer

CAIXA = [2.90, 41.755, 3.08, 41.88]     # conca del Ridaura i plana de Platja d'Aro
WORLDCOVER = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_N39E{:03d}_Map.tif"
DIRMAP = (1, 2, 4, 8, 16, 32, 64, 128)  # pysheds en l'ordre N, NE, E, SE, S, SO, O, NO → codi 0–7


FONS_LLERA, FONS_FRANJA, FRANJA_M = 5.0, 2.0, 80.0   # només per a les direccions de l'aigua; l'embassament usa el relleu real


def crema_rius(z, mar, tr, crs, conca):
    """Grava al relleu (només per a les direccions) els rius i rieres d'OSM. Torna el relleu gravat i la màscara de llera."""
    f = Path("data") / conca / "passos" / "osm.json"
    if not f.exists():
        print("Sense osm.json: no es graven els rius"); return z, None
    from shapely.geometry import LineString
    a_utm = Transformer.from_crs(4326, crs, always_xy=True).transform
    nf, nc = z.shape; zb = z.copy(); llera = np.zeros(z.shape, bool); n = 0
    for el in json.loads(f.read_text()).get("elements", []):
        t = el.get("tags", {})
        if t.get("waterway") not in ("river", "stream", "canal") or "geometry" not in el or len(el["geometry"]) < 2:
            continue
        g = LineString([a_utm(p["lon"], p["lat"]) for p in el["geometry"]]); prev = np.inf; n += 1
        for d in np.arange(0, g.length + 0.1, tr.a / 2):                # en el sentit del dibuix d'OSM, que és el de l'aigua
            q = g.interpolate(d); c, r = ~tr * (q.x, q.y); r, c = int(r), int(c)
            if not (0 <= r < nf and 0 <= c < nc) or mar[r, c]:
                continue
            with np.errstate(all="ignore"):
                v = np.nanmin(z[max(r - 1, 0): r + 2, max(c - 1, 0): c + 2])
            prev = min(prev, v); zb[r, c] = min(zb[r, c], prev - FONS_LLERA); llera[r, c] = True
    # franja als costats que baixa cap a la llera (com AGREE), perquè l'aigua del vessant hi entri i no corri en paral·lel
    from scipy.ndimage import distance_transform_edt
    dist = distance_transform_edt(~llera) * tr.a
    franja = ~llera & (dist < FRANJA_M)
    zb[franja] -= FONS_FRANJA * (1 - dist[franja] / FRANJA_M)
    print(f"Rius i rieres d'OSM gravats: {n} trams, {llera.sum()} cel·les ({FONS_LLERA:.0f} m de llera, "
          f"franja de {FRANJA_M:.0f} m que hi baixa {FONS_FRANJA:.0f} m)")
    return zb, llera


def dins_rius(acc, llera, cel_km2, dist_m=40, pas=20):
    """Part de les cel·les amb més d'1 km² de conca que queden a menys de dist_m d'una llera d'OSM."""
    from scipy.ndimage import distance_transform_edt
    fort = acc * cel_km2 > 1
    dist = distance_transform_edt(~llera) * pas
    return float((dist[fort] <= dist_m).mean()) if fort.any() else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--mdt", required=True)
    ap.add_argument("--radar", required=True)
    ap.add_argument("--episodi", required=True)
    ap.add_argument("--pas", type=float, default=20)
    a = ap.parse_args()
    cfg = json.loads(Path(a.config).read_text())
    global CAIXA
    CAIXA = cfg.get("caixa_web", cfg.get("caixa", CAIXA))
    dest = Path("data") / cfg["conca"] / "web"; dest.mkdir(parents=True, exist_ok=True)

    with rasterio.open(a.mdt) as f:
        w, s, e, n = transform_bounds("EPSG:4326", f.crs, *CAIXA)
        fin = from_bounds(w, s, e, n, f.transform).round_offsets().round_lengths()
        z0 = f.read(1, window=fin)
        tr0 = f.window_transform(fin)
        fac = max(int(round(a.pas / f.res[0])), 1)
        crs = f.crs
    nf, nc = z0.shape[0] // fac, z0.shape[1] // fac
    with np.errstate(all="ignore"):
        z = np.nanmin(np.nanmin(z0[: nf * fac, : nc * fac].reshape(nf, fac, nc, fac), axis=3), axis=1)
    tr = tr0 * tr0.scale(fac, fac)
    mar = ~np.isfinite(z)
    print(f"Malla: {nc}×{nf} cel·les de {tr.a:.0f} m")

    from pysheds.grid import Grid
    def direccions(zz, nom):
        tmp = Path(f"/tmp/relleu/{nom}.tif")
        with rasterio.open(tmp, "w", driver="GTiff", height=nf, width=nc, count=1, dtype="float32", crs=crs, transform=tr, nodata=-9999) as f:
            f.write(np.where(mar, -9999, zz).astype("float32"), 1)
        g = Grid.from_raster(str(tmp))
        om = g.fill_depressions(g.fill_pits(g.read_raster(str(tmp))))
        fdir = g.flowdir(g.resolve_flats(om), dirmap=DIRMAP)
        return om, np.asarray(fdir), np.asarray(g.accumulation(fdir, dirmap=DIRMAP))
    zb, llera = crema_rius(z, mar, tr, crs, cfg["conca"])
    omplert, fd, acc = direccions(z, "web")
    if llera is not None:
        abans = dins_rius(acc, llera, tr.a ** 2 / 1e6, pas=tr.a)
        _, fd, acc = direccions(zb, "web_rius")                        # direccions amb els rius gravats; l'embassament, amb el relleu real
        print(f"Cabal fort (> 1 km² de conca) a menys de 40 m d'un riu d'OSM: {abans * 100:.0f} % sense gravar, "
              f"{dins_rius(acc, llera, tr.a ** 2 / 1e6, pas=tr.a) * 100:.0f} % gravant-los")
    codi = np.full(fd.shape, 255, dtype=np.uint8)
    for k, v in enumerate(DIRMAP):
        codi[fd == v] = k
    codi[mar] = 255
    omplert = np.asarray(omplert, dtype="float64")
    zc = np.where(mar, 65535, np.clip(np.rint(z * 100), 0, 65534)).astype("<u2")
    zfc = np.where(mar, 65535, np.clip(np.rint(np.maximum(omplert, z) * 100), 0, 65534)).astype("<u2")
    prof = (zfc.astype(int) - zc.astype(int))
    print(f"Depressions (> 5 cm): {(prof[~mar] >= 5).sum() * tr.a ** 2 / 1e6:.2f} km²")

    urba = np.zeros((nf, nc), dtype=bool)
    for lon0 in (0, 3):
        try:
            with rasterio.open(WORLDCOVER.format(lon0)) as f:
                dst = np.zeros((nf, nc), dtype=np.uint8)
                reproject(rasterio.band(f, 1), dst, dst_transform=tr, dst_crs=crs, resampling=Resampling.mode)
                urba |= dst == 50
        except Exception as ex:
            print(f"WorldCover E{lon0:03d}: {type(ex).__name__} {ex}")
    print(f"Construït: {urba[~mar].mean() * 100:.1f} %")

    rad = json.loads(Path(a.radar).read_text())
    to_ll = Transformer.from_crs(crs, 4326, always_xy=True)
    X, Y = np.meshgrid(tr.c + (np.arange(nc) + 0.5) * tr.a, tr.f - (np.arange(nf) + 0.5) * tr.a)
    LON, LAT = to_ll.transform(X, Y)
    lons, lats = np.array(rad["lons"]), np.array(rad["lats"])
    ci = np.rint((LON - lons[0]) / rad["pas_graus"]).astype(int)
    ri = np.rint((LAT - lats[0]) / rad["pas_graus"]).astype(int)
    fora = (ci < 0) | (ci >= len(lons)) | (ri < 0) | (ri >= len(lats))
    ridx = np.where(fora, 65535, ri * len(lons) + ci).astype("<u2")

    # correcció del radar: mm del pluviòmetre / mm del radar a la conca de l'aforador, en l'hora que acaba a cada fotograma
    ep = json.loads(Path(a.episodi).read_text())
    t_ep = [datetime.fromisoformat(t) for t in ep["temps_utc"]]
    plv = np.array(ep["pluja_pluviometre_mm"], dtype=float); rdc = np.array(ep["pluja_radar_conca_mm"], dtype=float)
    factors = []
    for t in rad["temps"]:
        t = datetime.fromisoformat(t)
        sel = [k for k, te in enumerate(t_ep) if t - timedelta(minutes=60) < te <= t]
        p, r = np.nansum(plv[sel]), np.nansum(rdc[sel])
        factors.append(round(float(np.clip(p / r, 0.3, 4)), 3) if r > 0.5 else 1.0)

    capes, off, blob = [], 0, b""
    for nom, arr in (("z_cm", zc), ("zomplert_cm", zfc), ("dir", codi), ("urba", urba.astype(np.uint8)), ("radar", ridx)):
        b = np.ascontiguousarray(arr).tobytes()
        pad = (-len(blob)) % 4
        blob += b"\0" * pad; off = len(blob)
        capes.append({"nom": nom, "tipus": str(arr.dtype.str), "offset": off, "bytes": len(b)})
        blob += b
    (dest / "graella.bin").write_bytes(blob)
    w4, s4 = to_ll.transform(tr.c, tr.f - nf * tr.a)
    e4, n4 = to_ll.transform(tr.c + nc * tr.a, tr.f)
    to_utm = Transformer.from_crs(4326, crs, always_xy=True)
    punts = {}
    punts_cfg = {"aforador": (cfg["aforador"]["lon"], cfg["aforador"]["lat"]), "desembocadura": (3.0671, 41.8051),
                 "fenals": (3.0712, 41.8194)} if cfg["conca"] == "ridaura" else {k: tuple(v) for k, v in cfg["punts"].items()}
    for nom, (lon, lat) in punts_cfg.items():
        x, y = to_utm.transform(lon, lat)
        punts[nom] = [int((tr.f - y) // tr.a), int((x - tr.c) // tr.a)]
    (dest / "graella.json").write_text(json.dumps({
        "nf": nf, "nc": nc, "pas_m": tr.a, "crs": str(crs), "origen_utm": [tr.c, tr.f], "limits": [[s4, w4], [n4, e4]],
        "capes": capes, "dirs": "0=N 1=NE 2=E 3=SE 4=S 5=SO 6=O 7=NO, 255=mar o vora",
        "radar": Path(a.radar).name, "factor_radar": factors, "punts": punts}, ensure_ascii=False))
    print(f"graella.bin: {len(blob) / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
