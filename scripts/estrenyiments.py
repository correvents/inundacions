#!/usr/bin/env python3
"""Estrenyiments de llera que desborden primer: capacitat de cada tram de riera amb el MET de l'ICGC, comparada amb els
cabals de 10 a 500 anys (opció; res ajustat a mà).

Per a cada conca del tall de 100 km²:
1. Drenatge amb l'MDT 5 m (Release mdt5-cat): conca drenada, recorregut més llarg i desnivell de cada punt.
2. Cabals de T = 2…500 anys a cada punt de llera amb la 5.2-IC (com scripts/passos_cat.py de la PR #27): pluja de T anys
   en 24 h del Gumbel regional de la XEMA (estacions de la conca pilot, o a 20 km a la resta), Témez, Kt, Ka, i P0 =
   absorció de l'episodi més moll del model (pilots) o 0,2·S de les cobertes ICGC (resta).
   On el riu té ponts de l'inventari oficial de l'ACA (Punts Singulars Estructurals, PR #28), els cabals Q10/Q100/Q500
   de l'ACA (models hidràulics MAPRI) substitueixen els de la 5.2-IC al tram, escalats per (A/A_pont)^0,75 entre ponts.
3. Seccions transversals cada 20 m al llarg de les lleres (data/catalunya/relleu/cNNN/rieretes.json, ≥ --area-min km²),
   perpendiculars a la llera, amb el MET 2 m de l'ICGC (Release met21-cat: LiDAR 2021-2023, 1 m agregat a 2 m; si no
   n'hi ha, met2-cat 2008-2011; si tampoc, l'MDT 5 m, amb avís: metres enters).
   - Fons: el punt més baix a prop de la línia de l'MDT 5 m.
   - Nivell de desbordament: s'omple la secció des del fons de 5 en 5 cm (només l'aigua connectada amb el fons). L'aigua
     surt de la llera al primer nivell (≥ 0,5 m) on, pujant 50 cm més, avança ≥ 10 m per un costat (terreny de menys
     del 5 %) sobre conreus, carrers, edificis, indústria o zones verdes (cobertes ICGC), i l'amplada creix ≥ 50 %.
     Escampar-se pel llit ample on corre el canal de baix cabal (prats, matollar, bosc, ribera, sòl nu, aigua) no
     compta. Arribar al final de la secció (50 + 25·√A m a cada banda, A en km²) també compta. Si no passa, es pren el
     mínim d'amplada/fondària (Riley, 1972); si el mínim és al capdamunt (vall en V), la vall el conté.
   - Capacitat a sobreeixidor: Manning amb la n mitjana de les cobertes ICGC mullades (Chow 1959, els mateixos valors de
     config/cobertes_parametres.json que la simulació 2D; 0,035 sense cobertes) i el pendent del fons en 300 m.
   - Les seccions a menys de 30 m d'un pas sota via (PR #16) o d'un pont ACA prenen també la capacitat de l'obra quan
     se'n sap la mida real (data/catalunya/obres/cNNN.json de la PR #28: ACA, Generalitat o OSM; o la capacitat
     màxima de la fitxa ACA): el MET no veu els tubs. Les mides empíriques no es fan servir.
4. Període de retorn de desbordament de cada secció (interpolació log-log entre els cabals de T anys) i estrenyiments:
   trams de ≥ 3 seccions seguides (60 m) on el període és el mínim de la llera en 500 m i la capacitat cau a menys de la
   meitat de la de 200 m aigües amunt, o bé on la llera desborda abans de T10. Primer els que són dins o a tocar d'un
   nucli (cobertes urbanes a < 30 m): fora, una llera natural sol desbordar cada 1-2 anys (cabal de sobreeixidor).

Sortides: data/catalunya/estrenyiments/cNNN.json (resum, estrenyiments ordenats per període de desbordament, ponts ACA)
          data/catalunya/estrenyiments/index.json (resum de totes les conques)
          --seccions/cNNN.json.gz (totes les seccions; Release estrenyiments-cat)
Comprovació: zones inundables de l'ACA (T10-T500) a la caixa de la conca → valida (taula de períodes model/ACA).

Ús:  python3 scripts/estrenyiments.py conca --id 125 --mdt mdt5_c125.tif --met met21_c125.tif --cobertes cobertes5_c125.tif
     python3 scripts/estrenyiments.py zones --id 125 --sortida /tmp/zones_c125
     python3 scripts/estrenyiments.py valida --id 125 --zones /tmp/zones_c125
     python3 scripts/estrenyiments.py informe [--ids 125,188,199,124] [--sortida docs/estrenyiments.md]
     python3 scripts/estrenyiments.py index"""
import argparse, gzip, json, math, sys, tempfile
from pathlib import Path
import numpy as np
np.in1d = np.isin
import rasterio
from pyproj import Transformer
from scipy import ndimage
from shapely.geometry import shape, Point, LineString

CONQUES = Path("data/catalunya/conques/conques_100km2.geojson")
RIERETES = Path("data/catalunya/relleu/c{:03d}/rieretes.json")
PASSOS = Path("data/catalunya/passos/c{:03d}.geojson")
OBRES = Path("data/catalunya/obres/c{:03d}.json")
ACA = Path("data/catalunya/obres/aca_pse.json")
MAXIMS = Path("data/catalunya/xema/maxims_diaris.json")
SORT = Path("data/catalunya/estrenyiments")
CAD_IDX = Path("data/catalunya/cadastre/municipis.json")
CAD_URL = "https://github.com/correvents/inundacions/releases/download/cadastre-cat/edificis_{}.json.gz"
CAD_TMP = Path("/tmp/cadastre_cat")
RIU_GRAN = 1000.0                           # km²: a l'informe de Catalunya, a part (estudi ACA)
CONQUES_WEB = Path("data/catalunya/conques/web/conques_100.json")
CABAL_BASE = Path("data/catalunya/cabal_base/aforadors.json")
AIGUA = (20, 36, 37, 38, 39, 40, 41)       # cobertes ICGC: zones humides, embassaments, llacs, cursos d'aigua, basses, canals, mar
PARET = 1.0                                 # m: la façana fa de paret fins a 1 m d'aigua; per sobre entra a les plantes baixes (desborda)
SECCIONS = Path("/tmp/estrenyiments_seccions")
PILOTS = {124: "ridaura", 125: "ridaura", 188: "granollers", 199: "cabrera"}
TS = (2, 5, 10, 25, 50, 100, 200, 500)
N_MANNING = 0.035
PAS = 20.0            # m entre seccions
DH = 0.05             # m, pas d'ompliment
SALT = 0.5            # salt relatiu d'amplada que vol dir sortir de la llera
D_OBRA = 30.0         # m: secció dominada per una obra
ESTESA = 10.0         # m que l'aigua avança cap a un costat en 50 cm de pujada (pendent < 5 %) per sortir de la llera
MOTES = True          # desbordament per sobre de la mota (vegeu mota())
MOTA_L = 30.0         # m seguits de plana per sota de la carena, darrere la mota
MOTA_DZ = 0.5         # m que la plana protegida queda per sota de la carena
HMIN = 0.5            # m: per sota, els salts d'amplada són el canal de baix cabal dins del llit


# ---------- pluja de T anys (Gumbel regional de la XEMA, com passos_cat.py) ----------
def gumbel(norm, mitj):
    x = np.array(norm); beta = x.std(ddof=1) * np.sqrt(6) / np.pi; mu = x.mean() - 0.5772 * beta
    return {T: 1.13 * float(mitj) * (mu - beta * np.log(-np.log(1 - 1 / T))) for T in TS}


def p24_pilot(nom):
    import pandas as pd
    norm, mitj = [], []
    for f in sorted((Path("data") / nom / "xema").glob("*.csv")):
        if f.name.startswith("episodis"):
            continue
        s = pd.read_csv(f, parse_dates=["temps"]).set_index("temps").mm
        dia = s.resample("1D").agg(["sum", "count"])
        dia = dia[dia["count"] >= 44]
        anual = dia.groupby(dia.index.year)["sum"].agg(["max", "count"])
        anual = anual[anual["count"] >= 330]["max"]
        if len(anual) >= 5:
            mitj.append(anual.mean()); norm += list(anual / anual.mean())
    return (gumbel(norm, np.mean(mitj)), len(norm)) if len(norm) >= 20 else (None, len(norm))


def p24_punt(lon, lat, _c={}):
    if "e" not in _c:
        _c["e"] = [(r["lat"], r["lon"], list(r["anys"].values())) for r in json.loads(MAXIMS.read_text()).values()
                   if len(r["anys"]) >= 5] if MAXIMS.exists() else []
    km = sorted((math.hypot((la - lat) * 111.2, (lo - lon) * 111.2 * math.cos(math.radians(lat))), v) for la, lo, v in _c["e"])
    radi = 20
    while True:
        sel = [v for d, v in km if d <= radi]
        if sum(len(v) for v in sel) >= 20 or radi > 200:
            break
        radi += 10
    norm = [x / np.mean(v) for v in sel for x in v]
    return (gumbel(norm, np.mean([np.mean(v) for v in sel])), len(norm)) if len(norm) >= 20 else (None, len(norm))


def q52ic(area, l_km, j, p24, p0):
    """Cabal de la 5.2-IC (mètode racional modificat), com q_disseny() de passos_cat.py."""
    tc = max(0.3 * (l_km / max(j, 1e-3) ** 0.25) ** 0.76, 0.25)
    ka = 1 - np.log10(area) / 15 if area > 1 else 1.0
    pd_ = p24 * ka; id_ = pd_ / 24
    i = id_ * 11 ** ((28 ** 0.1 - tc ** 0.1) / (28 ** 0.1 - 1))
    kt = 1 + tc ** 1.25 / (tc ** 1.25 + 14)
    r = pd_ / p0
    c = (r - 1) * (r + 23) / (r + 11) ** 2 if r > 1 else 0.0
    return c * i * area * kt / 3.6


# ---------- drenatge (MDT 5 m) ----------
def entrades(cid):
    """Conques del tall de 100 km² que desguassen a la conca cid: [(lon, lat, km2_amunt)] dels seus punts de sortida.
    L'MDT de cada conca només cobreix la seva superfície: sense això, un riu que hi entra d'aigües amunt (Ter, Llobregat,
    Congost...) hi tindria només la conca local."""
    F = json.loads(CONQUES_WEB.read_text())["features"]
    return [(f["properties"]["sortida_lon"], f["properties"]["sortida_lat"], f["properties"]["km2_amunt"])
            for f in F if f["properties"].get("avall") == cid]


def drenatge(cami_mdt, cami_cob, entra=()):
    """Drenatge de l'MDT. entra: [(x, y, km2)] (coordenades de l'MDT) on entra un riu de conques d'aigües amunt; s'hi
    suma la superfície que falta fins a km2 a la cel·la de més acumulació a menys de 300 m. La llargada que hi falta
    s'estima per la llei de Hack (L = 1,4·A^0,6 km) i el desnivell amb un pendent mitjà de l'1 %."""
    import numba
    from pysheds.grid import Grid
    with rasterio.open(cami_mdt) as f:
        z = f.read(1).astype("float32"); nd = f.nodata; tr = f.transform; crs = f.crs
    z[z == nd] = np.nan; z /= 100; z[z <= 0.2] = np.nan
    tmp = Path(tempfile.mkdtemp()) / "z.tif"
    with rasterio.open(tmp, "w", driver="GTiff", height=z.shape[0], width=z.shape[1], count=1, dtype="float32", crs=crs,
                       transform=tr, nodata=np.nan) as f:
        f.write(z, 1)
    gr = Grid.from_raster(str(tmp)); dem = gr.read_raster(str(tmp))
    fd = gr.flowdir(gr.resolve_flats(gr.fill_depressions(gr.fill_pits(dem))))
    acc = np.asarray(gr.accumulation(fd), dtype="float64")
    acc_loc = acc.copy(); cel = abs(tr.a * tr.e); pas = abs(tr.a)
    extra = []
    for x, y, km2 in entra:
        c, r = ~tr * (x, y); r, c = int(r), int(c); n = int(300 / pas)
        r0, r1, c0, c1 = max(r - n, 0), min(r + n + 1, z.shape[0]), max(c - n, 0), min(c + n + 1, z.shape[1])
        if r0 >= r1 or c0 >= c1:
            continue
        sub = acc[r0:r1, c0:c1]; k = np.unravel_index(np.argmax(sub), sub.shape); r, c = r0 + k[0], c0 + k[1]
        falta = km2 * 1e6 / cel - acc[r, c]
        if falta > 0:
            extra.append((r, c, falta, km2))
    if extra:
        w = np.where(np.isfinite(z), 1.0, 0.0)
        for r, c, falta, _ in extra:
            w[r, c] += falta
        from pysheds.sview import Raster
        acc = np.asarray(gr.accumulation(fd, weights=Raster(w, viewfinder=dem.viewfinder)), dtype="float64")
    fd = np.asarray(fd).astype("int64")

    @numba.njit(cache=True)
    def passa(ordre, fd, lup, zmax, w, nf, nc):
        dr = {64: -1, 128: -1, 1: 0, 2: 1, 4: 1, 8: 1, 16: 0, 32: -1}
        dc = {64: 0, 128: 1, 1: 1, 2: 1, 4: 0, 8: -1, 16: -1, 32: -1}
        for k in ordre:
            r = k // nc; c = k % nc; d = fd[r, c]
            if d not in dr:
                continue
            r2 = r + dr[d]; c2 = c + dc[d]
            if r2 < 0 or r2 >= nf or c2 < 0 or c2 >= nc:
                continue
            p = 1.4142135 if dr[d] != 0 and dc[d] != 0 else 1.0
            if lup[r, c] + p > lup[r2, c2]:
                lup[r2, c2] = lup[r, c] + p
            if zmax[r, c] > zmax[r2, c2]:
                zmax[r2, c2] = zmax[r, c]
            w[r2, c2] += w[r, c]

    lup = np.zeros(z.shape, "float32"); zmax = np.where(np.isfinite(z), z, -1e9).astype("float32")
    for r, c, _, km2 in extra:
        L = 1.4 * km2 ** 0.6 * 1000
        lup[r, c] = max(lup[r, c], L / pas); zmax[r, c] = max(zmax[r, c], z[r, c] + 0.01 * L)
    S = np.zeros(z.shape, "float64"); urba = np.zeros(z.shape, bool)
    if cami_cob and Path(cami_cob).exists():
        sys.path.insert(0, str(Path(__file__).parent))
        import cobertes as C
        from rasterio.warp import reproject, Resampling
        with rasterio.open(cami_cob) as src:
            cob = np.zeros(z.shape, np.uint8)
            reproject(rasterio.band(src, 1), cob, dst_transform=tr, dst_crs=crs, src_nodata=0, dst_nodata=0,
                      resampling=Resampling.mode)
        _, cn, tt = C.parametres(cob)
        urba = ndimage.binary_dilation(tt == 1, iterations=6)          # cel·la urbana a menys de 30 m
        S = np.where(np.isfinite(z), 25400 / np.clip(cn, 1, 100) - 254, 0).astype("float64")
    passa(np.argsort(acc, axis=None, kind="stable"), fd, lup, zmax, S, *z.shape)
    return {"z": z, "tr": tr, "crs": crs, "acc": acc, "acc_loc": acc_loc, "lup": lup, "zmax": zmax, "S": S, "urba": urba,
            "entrades": len(extra)}


def punt_drenatge(D, x, y, radi_cel=3):
    """Cel·la de màxima acumulació a prop de (x, y) (coordenades de l'MDT): conca km², L km, J, S mitjana."""
    tr = D["tr"]; c, r = ~tr * (x, y); r, c = int(r), int(c)
    nf, nc = D["z"].shape
    r0, r1, c0, c1 = max(r - radi_cel, 0), min(r + radi_cel + 1, nf), max(c - radi_cel, 0), min(c + radi_cel + 1, nc)
    if r0 >= r1 or c0 >= c1:
        return None
    sub = D["acc"][r0:r1, c0:c1]; k = np.unravel_index(np.argmax(sub), sub.shape); r, c = r0 + k[0], c0 + k[1]
    acc = float(D["acc"][r, c]); cel = abs(tr.a * tr.e)
    if not np.isfinite(D["z"][r, c]):
        return None
    l_km = max(float(D["lup"][r, c]) * abs(tr.a) / 1000, 0.05)
    j = (float(D["zmax"][r, c]) - float(D["z"][r, c])) / (l_km * 1000)
    return {"km2": acc * cel / 1e6, "l_km": l_km, "j": j, "S": float(D["S"][r, c]) / max(float(D["acc_loc"][r, c]), 1), "urba": bool(D["urba"][r, c])}


# ---------- secció transversal (MET) ----------
class Relleu:
    def __init__(self, cami, nom):
        self.nom = nom
        with rasterio.open(cami) as f:
            self.a = f.read(1).astype("float32"); self.tr = f.transform; self.crs = f.crs; nd = f.nodata
        self.a[self.a == nd] = np.nan; self.a /= 100
        self.a[self.a <= -100] = np.nan
        self.inv = ~self.tr
        self.nan = np.isnan(self.a)
        self.a[self.nan] = -9999
        self.edif = None

    def edificis(self, cad_idx=CAD_IDX):
        """Màscara dels edificis del Cadastre (parts amb ≥ 1 planta sobre rasant; Release cadastre-cat) a la reixa del
        relleu. El MET és terreny nu: sense això l'aigua s'escamparia per dins de les illes de cases."""
        import gzip, requests
        from rasterio import features
        if not Path(cad_idx).exists() or not str(self.crs).endswith("25831"):
            return 0
        h, w = self.a.shape
        x0, y1 = self.tr.c, self.tr.f; x1, y0 = x0 + w * self.tr.a, y1 + h * self.tr.e
        formes = []
        CAD_TMP.mkdir(parents=True, exist_ok=True)
        for codi, m in json.loads(Path(cad_idx).read_text()).items():
            b = m["bbox"]
            if b[2] < x0 or b[0] > x1 or b[3] < y0 or b[1] > y1:
                continue
            f = CAD_TMP / f"edificis_{codi}.json.gz"
            if not f.exists():
                for intent in range(4):
                    try:
                        r = requests.get(CAD_URL.format(codi), timeout=300); r.raise_for_status(); f.write_bytes(r.content); break
                    except requests.RequestException as e:
                        print(f"  Cadastre {codi}: {e}")
                else:
                    raise RuntimeError(f"no puc baixar els edificis de {codi}")
            for ref, sobre, sota, an in json.loads(gzip.open(f, "rt").read())["parts"]:
                if sobre < 1:
                    continue
                xa, ya = np.min(an[0], axis=0); xb, yb = np.max(an[0], axis=0)
                if xb < x0 or xa > x1 or yb < y0 or ya > y1:
                    continue
                formes.append(({"type": "Polygon", "coordinates": an}, 1))
        if formes:
            self.edif = features.rasterize(formes, out_shape=self.a.shape, transform=self.tr, fill=0, dtype="uint8").astype(bool)
        return len(formes)

    def es_edifici(self, xs, ys):
        if self.edif is None:
            return np.zeros(len(xs), bool)
        c, r = self.inv * (np.asarray(xs), np.asarray(ys))
        r = np.clip(r.astype(int), 0, self.a.shape[0] - 1); c = np.clip(c.astype(int), 0, self.a.shape[1] - 1)
        return self.edif[r, c]

    def mostra(self, xs, ys):
        c, r = self.inv * (np.asarray(xs), np.asarray(ys))
        v = ndimage.map_coordinates(self.a, [r - 0.5, c - 0.5], order=1, mode="constant", cval=-9999)
        ri, ci = np.floor(r - 0.5).astype(int), np.floor(c - 0.5).astype(int)
        dins = (ri >= 0) & (ci >= 0) & (ri < self.a.shape[0] - 1) & (ci < self.a.shape[1] - 1)
        mal = ~dins
        rr, cc = np.clip(ri, 0, self.a.shape[0] - 2), np.clip(ci, 0, self.a.shape[1] - 2)
        mal |= self.nan[rr, cc] | self.nan[rr + 1, cc] | self.nan[rr, cc + 1] | self.nan[rr + 1, cc + 1]
        v[mal] = np.nan
        return v


# cobertes ICGC on l'aigua no "surt de la llera": boscos, matollar, prats, ribera, sòl nu forestal, roquissars, platges,
# zones humides i aigua (7-20, 36-41). Desbordar és arribar a conreus, carrers, edificis, indústria o zones verdes.
LLIT = tuple(range(7, 21)) + tuple(range(36, 42))


def omple(s, zs, i0, llit=None, hmax=20.0, km2=0.0):
    """Ompliment de la secció des del fons i0. Torna el nivell de desbordament, l'àrea i el perímetre mullats, el
    costat per on surt i si la vall el conté."""
    zt = zs[i0]; ds = s[1] - s[0]; n = len(zs)
    hist = []
    ML = np.maximum.accumulate(zs[i0::-1])            # cota màxima entre el fons i cada punt, cap a l'esquerra
    MR = np.maximum.accumulate(zs[i0:])               # i cap a la dreta
    for h in np.arange(DH, hmax, DH):
        nivell = zt + h
        l = i0 - (int(np.searchsorted(ML, nivell, side="left")) - 1)
        r = i0 + (int(np.searchsorted(MR, nivell, side="left")) - 1)
        tocat = l == 0 or r == n - 1
        prof = np.clip(nivell - zs[l:r + 1], 0, None)
        A = float(prof.sum() * ds)
        W = (r - l + 1) * ds
        P = float(np.sum(np.hypot(ds, np.diff(zs[l:r + 1])))) + 2 * 0.5 * ds
        hist.append((h, l, r, W, A, P, tocat))
        if tocat:
            break
    if len(hist) < 3:
        return None
    # sortida de la llera: el primer nivell h (≥ HMIN) on, 50 cm més amunt, l'aigua ha avançat ≥ ESTESA m per un dels
    # costats (terreny més pla que el 5 %: una plana, un camp o un carrer) i l'amplada ha crescut ≥ 50 %. Les parets de
    # la llera, encara que siguin tendides (2:1), no hi arriben. Arribar al final de la secció també compta.
    k5 = int(round(0.5 / DH))
    millor = None
    fi = hist[-1][6]
    for k in range(len(hist) - (0 if fi else k5)):
        h, l, r, W, A, P, _ = hist[k]
        if k + k5 < len(hist):
            _, l2, r2, W2, _, _, _ = hist[k + k5]
        else:
            l2, r2, W2 = (0 if hist[-1][1] == 0 else l), (len(s) - 1 if hist[-1][2] == len(s) - 1 else r), len(s) * ds
        fora = (lambda a_, b_: (b_ - a_) if llit is None else int((~llit[a_:b_]).sum()))   # metres nous fora del llit
        if h >= HMIN - 1e-9 and W2 >= (1 + SALT) * W and max(fora(l2, l), fora(r + 1, r2 + 1)) * ds >= ESTESA:
            if MOTES and terrassa(zs, i0, l2, r2, zt + h, ds, km2):
                continue                    # s'estén per una terrassa tancada dins la llera gran: encara no desborda
            millor = (W2 - W, k)
            break
    if millor and MOTES:
        # motes: si l'aigua surt entre motes (conreus, prats o camins dins de la llera gran), no desborda fins que
        # passa la mota més baixa. A l'altre costat, la mota o la vessant (cota màxima de la secció) si és més alta.
        mL = mota(zs[i0::-1], ds, HMIN); mR = mota(zs[i0:], ds, HMIN)
        if mL is not None or mR is not None:
            cL = mL if mL is not None else (ML[-1] if ML[-1] > (mR or 0) else None)
            cR = mR if mR is not None else (MR[-1] if MR[-1] > (mL or 0) else None)
            if cL is not None and cR is not None:
                hm = min(cL, cR) - zt
                km = max((i for i, t in enumerate(hist) if t[0] <= hm + 1e-9), default=None)
                if km is not None and km > millor[1]:
                    h, l, r, W, A, P, _ = hist[km]
                    nx = hist[min(km + 1, len(hist) - 1)]
                    costat = "esquerre" if cL <= cR else "dret"
                    return {"h": h, "A": A, "P": P, "W": W, "costat": costat, "conte": False, "l": l, "r": r, "mota": True}
    if millor:
        k = millor[1]
        h, l, r, W, A, P, _ = hist[k]
        # el primer pas on l'amplada salta, dins dels 50 cm
        k2 = k
        for k2 in range(k, min(k + k5, len(hist) - 1)):
            if hist[k2 + 1][3] - hist[k2][3] >= 0.5 * millor[0] / k5 or hist[k2 + 1][3] > (1 + SALT) * hist[k2][3]:
                break
        h, l, r, W, A, P, _ = hist[k2]
        l2, r2 = hist[min(k2 + 1, len(hist) - 1)][1], hist[min(k2 + 1, len(hist) - 1)][2]
        costat = "esquerre" if (l - l2) > (r2 - r) else "dret"
        return {"h": h, "A": A, "P": P, "W": W, "costat": costat, "conte": False, "l": l, "r": r}
    # sense salt: mínim de W/D (Riley). Si el mínim és al final (vall en V que s'estreny), la vall el conté.
    wd = [(W * W / A if A > 0 else 1e9) for (_, _, _, W, A, _, _) in hist]
    k = int(np.argmin([v if hist[i][0] >= 0.3 else 1e9 for i, v in enumerate(wd)]))
    h, l, r, W, A, P, _ = hist[k]
    conte = k >= len(hist) - 2
    return {"h": h, "A": A, "P": P, "W": W, "costat": None, "conte": conte, "l": l, "r": r}


def terrassa(zs, i0, l2, r2, nivell, ds, km2):
    """Rius grans: l'aigua que s'estén fins a l2..r2 és dins d'una terrassa de la llera (parc fluvial, prats, camins entre
    murs o motes) si, a cada costat, a menys de 3·√A m del fons (A en km²) hi ha un marge ≥ 1 m per sobre del nivell."""
    D = 3 * math.sqrt(max(km2, 0)) / ds
    if D * ds < 2 * ESTESA:
        return False
    a, b = max(int(i0 - D), 0), min(int(i0 + D), len(zs) - 1)
    return l2 > a and r2 < b and zs[a:l2 + 1].max() >= nivell + 1 and zs[r2:b + 1].max() >= nivell + 1


def mota(z, ds, hmin):
    """Cota de la primera mota cap a fora des del fons (z[0] = fons): carena ≥ hmin sobre el fons amb, darrere, ≥ MOTA_L m
    seguits més de 0,3 m per sota i un punt ≥ MOTA_DZ m per sota (plana protegida), dins de 150 m. Si no n'hi ha, None."""
    M = np.maximum.accumulate(z); nL = max(int(MOTA_L / ds), 1); n150 = int(150 / ds)
    for j in range(1, len(z) - nL):
        if z[j] < M[j] or z[j] - z[0] < hmin:
            continue
        rere = z[j + 1:j + 1 + n150]
        sota = rere < z[j] - 0.3
        if not sota.any() or rere.min() > z[j] - MOTA_DZ:
            continue
        run = np.diff(np.flatnonzero(np.diff(np.r_[0, sota.astype(int), 0])))[::2]   # trams seguits per sota
        if run.size and run.max() >= nL:
            return float(z[j])
    return None


def t_desborda(q, qs):
    """Període de retorn (anys) al qual el cabal arriba a q (log-log entre TS)."""
    v = [(T, qs[T]) for T in TS if qs.get(T, 0) > 0]
    if not v:
        return None
    if q <= v[0][1]:
        return float(v[0][0]) if v[0][0] > 2 else 2.0
    if q >= v[-1][1]:
        return float(TS[-1]) + 1
    for (t0, q0), (t1, q1) in zip(v, v[1:]):
        if q0 <= q <= q1:
            f = (math.log(q) - math.log(q0)) / max(math.log(q1) - math.log(q0), 1e-9)
            return float(math.exp(math.log(t0) + f * (math.log(t1) - math.log(t0))))
    return None


def densifica(xy, pas):
    ln = LineString(xy)
    d = np.arange(0, ln.length, pas)
    return ln, [ln.interpolate(v) for v in d], d


# ---------- una conca ----------
def conca(cid, a):
    global MOTES
    MOTES = getattr(a, "motes", True)
    Cq = {f["properties"]["id"]: f for f in json.loads(CONQUES.read_text())["features"]}
    geom = shape(Cq[cid]["geometry"])
    t_mdt = Transformer.from_crs(4326, rasterio.open(a.mdt).crs, always_xy=True).transform
    D = drenatge(a.mdt, a.cobertes, [(*t_mdt(lo, la), k) for lo, la, k in entrades(cid)])
    if D["entrades"]:
        print(f"Rius que entren d'aigües amunt: {D['entrades']}")
    met = Relleu(a.met, a.font_met) if a.met and Path(a.met).exists() else None
    n_edif = met.edificis() if (met is not None and a.edificis) else 0
    if n_edif:
        print(f"Edificis del Cadastre: {n_edif} parts, {met.edif.sum() * abs(met.tr.a * met.tr.e) / 1e6:.2f} km²")
    base = json.loads(CABAL_BASE.read_text())["aforadors"] if CABAL_BASE.exists() else []
    # fora sortides d'embassament i valors que no poden ser (aforador mal situat o conca mal trobada)
    base = [b for b in base if b.get("l_s_km2") is not None and 0 <= b["l_s_km2"] <= 50 and b["km2"] >= 5
            and "embassament" not in b["nom"].lower()]
    qb_mediana = float(np.median([b["l_s_km2"] for b in base])) if base else 0.0
    if base:
        print(f"Cabal base: {len(base)} aforadors ACA")
    mdt = Relleu(a.mdt, "MDT 5 m IGN (metres enters)")
    cobs = None
    if a.cobertes and Path(a.cobertes).exists():
        with rasterio.open(a.cobertes) as f:
            cobs = (f.read(1), f.transform)
            sys.path.insert(0, str(Path(__file__).parent))
            import cobertes as C
            cobs = cobs + (C.parametres(np.arange(256, dtype=np.uint8))[0],)      # n de Manning de cada classe
            met_a_cob = Transformer.from_crs((met or mdt).crs, f.crs, always_xy=True).transform
    W84 = 4326
    t_met = Transformer.from_crs(W84, (met or mdt).crs, always_xy=True).transform
    met_a_mdt = Transformer.from_crs((met or mdt).crs, D["crs"], always_xy=True).transform
    met_a_ll = Transformer.from_crs((met or mdt).crs, W84, always_xy=True).transform
    a_met25831 = Transformer.from_crs(25831, (met or mdt).crs, always_xy=True).transform

    # pluja
    nom = PILOTS.get(cid)
    p24, n_ae, p0 = None, 0, None
    if nom:
        p24, n_ae = p24_pilot(nom)
        idx = Path("data") / nom / "episodis" / "index.json"
        ias = [e["ia_mm"] for e in json.loads(idx.read_text())["episodis"] if e.get("ia_mm")] if idx.exists() else []
        p0 = float(min(ias)) if ias else None
        font_pluja = f"XEMA de la conca pilot {nom} ({n_ae} anys-estació)"
    if not p24:
        c = geom.centroid
        p24, n_ae = p24_punt(c.x, c.y)
        font_pluja = f"XEMA a 20 km ({n_ae} anys-estació)"
    font_pluja += f", P0 = {p0:.0f} mm (episodi més moll del model)" if p0 else ", P0 = 0,2·S de les cobertes"
    print(f"== Conca {cid}: P24 T10/T100/T500 = {p24[10]:.0f}/{p24[100]:.0f}/{p24[500]:.0f} mm; {font_pluja}")

    # obres: passos sota vies (PR #16/#27) amb capacitat (PR #28) i ponts ACA
    obres = []
    cap = json.loads(OBRES.with_name(OBRES.name.format(cid)).read_text()) if OBRES.with_name(OBRES.name.format(cid)).exists() else {}
    pf = Path(str(PASSOS).format(cid))
    if pf.exists():
        for f in json.loads(pf.read_text())["features"]:
            p = f["properties"]
            if f["geometry"]["type"] != "Point" or p.get("tipus") != "terraplè":
                continue
            o = cap.get(p["id"], {})
            # només capacitats d'obres amb mida real (ACA, Generalitat, OSM): l'ajust empíric i el cabal de disseny
            # 5.2-IC suposen la resposta (desbordar a T100) i no diuen res de l'obra concreta
            if not o.get("font") or o["font"] == "empíric" or not o.get("q_obert_m3s"):
                continue
            obres.append({"tipus": "pas", "id": p["id"], "nom": p.get("nom") or p.get("via"), "q": o["q_obert_m3s"],
                          "km2": p.get("conca_km2"), "q_font": o["font"], "xy": t_met(*f["geometry"]["coordinates"])})
    aca = []
    if ACA.exists():
        for p in json.loads(ACA.read_text())["pse"]:
            x, y = p["utm"]
            if geom.contains(Point(Transformer.from_crs(25831, W84, always_xy=True).transform(x, y))):
                qmax = p.get("q_max_m3s") or p.get("q_max_min_m3s")
                aca.append({"tipus": "pont ACA", "id": p["codi"], "nom": f"{p['tipus']} {p['via']} sobre {p['riu']}",
                            "riu": p["riu"], "q": qmax, "q_supera_q500": bool(p.get("q_max_m3s_supera")),
                            "q_font": "fitxa ACA PSE (capacitat màxima de desguàs)", "xy": a_met25831(x, y),
                            "qT": {10: p["t10"].get("q_estructura"), 100: p["t100"].get("q_estructura"),
                                   500: p["t500"].get("q_estructura")},
                            "capacitat_aca": {k: p[k].get("capacitat") for k in ("t10", "t100", "t500")}})
    obres_tot = obres + aca
    print(f"{len(obres)} passos sota via, {len(aca)} ponts ACA")

    # lleres
    j = json.loads(Path(str(RIERETES).format(cid)).read_text())
    esc = j["escala"]
    seccions = []
    n_coberta = 0
    for it, tr in enumerate(j["trams"]):
        if tr[1] / 100 < a.area_min:
            continue
        ll = np.cumsum(np.array(tr[2:], dtype=np.int64).reshape(-1, 2), axis=0) / esc
        if len(ll) < 2:
            continue
        xy = np.column_stack(t_met(ll[:, 0], ll[:, 1]))
        ln, pts, dist = densifica(xy, PAS)
        for k, (pt, d) in enumerate(zip(pts, dist)):
            p1 = ln.interpolate(max(d - 15, 0)); p2 = ln.interpolate(min(d + 15, ln.length))
            dx, dy = p2.x - p1.x, p2.y - p1.y; nrm = math.hypot(dx, dy) or 1
            ux, uy = -dy / nrm, dx / nrm                                  # perpendicular (esquerra, mirant aigües avall)
            dr = punt_drenatge(D, *met_a_mdt(pt.x, pt.y))
            if not dr or dr["km2"] < a.area_min * 0.5:
                continue
            A = dr["km2"]
            L = min(500.0, 50 + 25 * math.sqrt(A)); ws = 4 + 1.5 * math.sqrt(A)
            s = np.arange(-L, L + 0.01, 1.0)
            font = None
            for R in ([met] if met else []) + [mdt]:
                zs = R.mostra(pt.x + ux * s, pt.y + uy * s)
                centre = np.abs(s) <= L / 2
                if np.isfinite(zs[centre]).all():
                    font = R
                    break
            if font is None:
                continue
            # retalla les puntes sense dades
            ok = np.isfinite(zs)
            if not ok.all():
                i_c = len(s) // 2
                l = i_c
                while l > 0 and ok[l - 1]:
                    l -= 1
                r = i_c
                while r < len(s) - 1 and ok[r + 1]:
                    r += 1
                s, zs = s[l:r + 1], zs[l:r + 1]
            llit = None; nv = None
            if cobs is not None:
                cx, cy = met_a_cob(pt.x + ux * s, pt.y + uy * s)
                cc, rr = ~cobs[1] * (np.asarray(cx), np.asarray(cy))
                rr = np.clip(rr.astype(int), 0, cobs[0].shape[0] - 1); cc = np.clip(cc.astype(int), 0, cobs[0].shape[1] - 1)
                cl = cobs[0][rr, cc]
                llit = np.isin(cl, LLIT)
                nv = cobs[2][cl]
            cen = np.where(np.abs(s) <= ws)[0]
            eb = font.es_edifici(pt.x + ux * s, pt.y + uy * s) if font is met and met.edif is not None else np.zeros(len(s), bool)
            zc = np.where(eb, np.inf, zs)
            i0 = int(cen[np.argmin(zc[cen])])
            if not np.isfinite(zc[i0]):
                n_coberta += 1                                      # llera coberta per edificis: el tall no la veu
                continue
            # làmina d'aigua: franja plana (± 5 cm) al fons, sobre cobertes d'aigua, de ≥ 4 m. El LiDAR no hi entra.
            ll_, rr_ = i0, i0
            while ll_ > 0 and abs(zs[ll_ - 1] - zs[i0]) <= 0.05 and not eb[ll_ - 1]:
                ll_ -= 1
            while rr_ < len(zs) - 1 and abs(zs[rr_ + 1] - zs[i0]) <= 0.05 and not eb[rr_ + 1]:
                rr_ += 1
            Ww = (rr_ - ll_ + 1) * 1.0
            aigua = Ww >= 4 and cobs is not None and bool(np.isin(cl[ll_:rr_ + 1], AIGUA).any())
            res = omple(s, np.where(eb, zs + PARET, zs), i0, llit, km2=A)
            if not res:
                continue
            res["Ww"] = Ww if aigua else 0.0
            res["parets"] = int(eb[res["l"]:res["r"] + 1].sum() == 0 and eb.any())
            xf, yf = pt.x + ux * s[i0], pt.y + uy * s[i0]
            # n de Manning: mitjana de les cobertes mullades (Chow, config/cobertes_parametres.json, com la simulació 2D)
            res["n"] = float(np.mean(nv[res["l"]:res["r"] + 1])) if nv is not None else N_MANNING
            seccions.append({"tram": it, "k": k, "d": float(d), "x": xf, "y": yf, "zf": float(zs[i0]), "km2": A,
                             "l_km": dr["l_km"], "j": dr["j"], "S": dr["S"], "urba": dr["urba"], "font": font.nom, **res})
    print(f"{len(seccions)} seccions ({n_coberta} sota edificis, {sum(1 for sc in seccions if sc['Ww'] > 0)} amb làmina d'aigua)")
    if not seccions:
        return desa(cid, {"km2": Cq[cid]["properties"]["km2"], "seccions": 0}, [], [])

    # pendent del fons en 300 m (regressió robusta) i cabals
    from collections import defaultdict
    per_tram = defaultdict(list)
    for i, sc in enumerate(seccions):
        per_tram[sc["tram"]].append(i)
    for it, ids in per_tram.items():
        d = np.array([seccions[i]["d"] for i in ids]); zf = np.array([seccions[i]["zf"] for i in ids])
        for i, di in zip(ids, d):
            m = np.abs(d - di) <= 150
            if m.sum() >= 3:
                b = np.polyfit(d[m], zf[m], 1)[0]
            else:
                b = -seccions[i]["j"]
            seccions[i]["pend"] = float(np.clip(-b, 5e-4, 0.3))
    for sc in seccions:
        R = sc["A"] / max(sc["P"], 1e-6)
        if sc["Ww"] > 0 and base:
            # sota la làmina: fondària que fa córrer el cabal base (mediana diària dels aforadors ACA, per km² de
            # conca, del més proper a menys de 40 km) amb Manning en un canal ample de l'amplada de la làmina (n de curs d'aigua 0,03)
            lon, lat = met_a_ll(sc["x"], sc["y"])
            dk, g = min(((111.2 * math.hypot((b["lon"] - lon) * math.cos(math.radians(lat)), b["lat"] - lat), b)
                         for b in base), key=lambda t: t[0])
            qb = (g["l_s_km2"] if dk <= 40 else qb_mediana) / 1000 * sc["km2"]   # a més de 40 km, la mediana
            dsub = (qb * 0.03 / (sc["Ww"] * math.sqrt(sc["pend"]))) ** 0.6 if qb > 0 else 0.0
            sc["A"] += sc["Ww"] * dsub; sc["P"] += 2 * dsub
            sc["d_sota_lamina"] = dsub; sc["q_base"] = qb
            R = sc["A"] / max(sc["P"], 1e-6)
        sc["q_cap"] = float(sc["A"] * R ** (2 / 3) * math.sqrt(sc["pend"]) / sc["n"])
        pp0 = p0 if p0 else 0.2 * sc["S"]
        sc["p0"] = pp0
        sc["qT"] = {T: q52ic(sc["km2"], sc["l_km"], sc["j"], p24[T], pp0) for T in TS}
        sc["qT_font"] = "5.2-IC"

    # cabals ACA als rius ARPSI: el pont ACA s'associa a la secció més propera (< 80 m), i el factor ACA/5.2-IC
    # s'estén al llarg del tram i dels trams connectats mentre la conca no canviï més d'un factor 2
    if aca:
        xy = np.array([(sc["x"], sc["y"]) for sc in seccions])
        ancores = []
        for o in aca:
            if not all(o["qT"].get(T) for T in (10, 100, 500)):
                continue
            dd = np.hypot(xy[:, 0] - o["xy"][0], xy[:, 1] - o["xy"][1])
            mateix = np.array([1 / 3 <= o["qT"][100] / max(sc["qT"][100], 1e-3) <= 3 for sc in seccions])
            if (mateix & (dd <= 80)).any():
                i = int(np.argmin(np.where(mateix, dd, 1e12)))
            elif dd.min() <= 30:                    # conca drenada mal tallada (riu que ve d'una altra conca del tall)
                i = int(np.argmin(dd))
            else:
                continue
            sc = seccions[i]
            f = {T: o["qT"][T] / max(sc["qT"][T], 1e-3) for T in (10, 100, 500)}
            ancores.append((i, sc["km2"], f, o, sc["tram"]))
            o["seccio"] = i
            o["q52ic"] = {T: round(sc["qT"][T], 1) for T in (10, 100, 500)}
        if ancores:
            ak = np.array([an[1] for an in ancores]); axy = xy[[an[0] for an in ancores]]
            atr = np.array([an[4] for an in ancores])
            for i, sc in enumerate(seccions):
                dd = np.hypot(axy[:, 0] - sc["x"], axy[:, 1] - sc["y"])
                rat = np.maximum(ak / sc["km2"], sc["km2"] / ak)
                ok = (rat <= 2) & (dd <= 3000) & (atr == sc["tram"])     # només al mateix tram de riu
                if not ok.any():
                    continue
                jn = int(np.argmin(np.where(ok, dd, 1e12)))
                _, k2, f, o, _ = ancores[jn]
                qa = {T: o["qT"][T] * (sc["km2"] / k2) ** 0.75 for T in (10, 100, 500)}
                # T intermedis: log-log sobre T amb els tres de l'ACA
                lt = np.log([10, 100, 500]); lq = np.log([max(qa[T], 1e-3) for T in (10, 100, 500)])
                sc["qT"] = {T: float(np.exp(np.interp(np.log(T), lt, lq, left=None, right=None)) if 10 <= T <= 500 else
                                     np.exp(lq[0] + (np.log(T) - lt[0]) * (lq[1] - lq[0]) / (lt[1] - lt[0]))) for T in TS}
                sc["qT_font"] = f"ACA {o['id']}"

    # obres: capacitat de l'obra a les seccions properes
    for sc in seccions:
        sc["obra"] = None
    if obres_tot:
        for o in obres_tot:
            if not o.get("q"):
                continue
            dd = [(math.hypot(sc["x"] - o["xy"][0], sc["y"] - o["xy"][1]), i) for i, sc in enumerate(seccions)]
            # la secció més propera de la mateixa llera (conca semblant: un pas d'un afluent no és el del riu)
            dd = [(d, i) for d, i in dd if not o.get("km2") or 1 / 3 <= seccions[i]["km2"] / max(o["km2"], 1e-3) <= 3] or [(1e9, 0)]
            d, i = min(dd)
            if d <= D_OBRA:
                sc = seccions[i]
                if sc["obra"] is None or o["q"] < sc["obra"]["q"]:
                    sc["obra"] = {k: o[k] for k in ("tipus", "id", "nom", "q", "q_font")}

    # suavitzat: mediana de 3 seccions seguides (60 m) de la capacitat de llera; període de desbordament
    for it, ids in per_tram.items():
        q = np.array([seccions[i]["q_cap"] for i in ids])
        qm = ndimage.median_filter(q, size=3, mode="nearest") if len(q) >= 3 else q
        for i, v in zip(ids, qm):
            seccions[i]["q_llera"] = float(v)
    for sc in seccions:
        sc["T_llera"] = t_desborda(sc["q_llera"], sc["qT"]) if not sc["conte"] else float(TS[-1]) + 1
        sc["T_obra"] = t_desborda(sc["obra"]["q"], sc["qT"]) if sc["obra"] else None

    # estrenyiments
    estr = []
    for it, ids in per_tram.items():
        d = np.array([seccions[i]["d"] for i in ids])
        T = np.array([seccions[i]["T_llera"] or 501 for i in ids])
        q = np.array([seccions[i]["q_llera"] for i in ids])
        k = 0
        while k < len(ids):
            # mínim local en 500 m
            m = np.abs(d - d[k]) <= 250
            amunt = (d < d[k] - 20) & (d >= d[k] - 220)
            qa = float(np.median(q[amunt])) if amunt.any() else None
            es_min = T[k] <= T[m].min() + 1e-9
            caiguda = qa is not None and q[k] < 0.5 * qa
            if es_min and T[k] <= 500 and (caiguda or T[k] < 10):
                # el tram de seccions seguides amb període ≤ 1,5·T[k]
                a0 = k
                while a0 > 0 and T[a0 - 1] <= 1.5 * T[k]:
                    a0 -= 1
                b0 = k
                while b0 < len(ids) - 1 and T[b0 + 1] <= 1.5 * T[k]:
                    b0 += 1
                if b0 - a0 + 1 >= 3:
                    sc = seccions[ids[k]]
                    estr.append({"i": ids[k], "tram": it, "llarg_m": float(d[b0] - d[a0] + PAS),
                                 "q_amunt_m3s": round(qa, 1) if qa else None})
                k = b0 + 1
                continue
            k += 1
    for o in obres_tot:
        pass
    # obres que desborden (capacitat < cabal) com a estrenyiments puntuals
    for i, sc in enumerate(seccions):
        if sc["obra"] and sc["T_obra"] is not None and sc["T_obra"] <= 500 and sc["T_obra"] < (sc["T_llera"] or 501):
            estr.append({"i": i, "tram": sc["tram"], "llarg_m": 0, "obra": True})

    def fitxa(e):
        sc = seccions[e["i"]]
        lon, lat = met_a_ll(sc["x"], sc["y"])
        ob = e.get("obra")
        T = sc["T_obra"] if ob else sc["T_llera"]
        q = sc["obra"]["q"] if ob else sc["q_llera"]
        return {"lonlat": [round(lon, 5), round(lat, 5)], "tipus": ("obra: " + sc["obra"]["tipus"]) if ob else "llera",
                "obra": sc["obra"] if ob else None, "conca_km2": round(sc["km2"], 2), "capacitat_m3s": round(q, 1),
                "T_desborda": round(T) if T <= 500 else ">500", "Q10": round(sc["qT"][10], 1), "Q100": round(sc["qT"][100], 1),
                "Q500": round(sc["qT"][500], 1), "qT_font": sc["qT_font"], "nivell_m": round(sc["h"], 2),
                "amplada_m": round(sc["W"], 1), "pendent": round(sc["pend"], 4), "surt_per": sc["costat"],
                "llarg_m": e["llarg_m"], "q_amunt_m3s": e.get("q_amunt_m3s"), "font_relleu": sc["font"], "urba": sc["urba"]}
    # primer els de dins o a tocar de nuclis (cobertes urbanes a < 30 m): fora, desbordar cada 2 anys és el normal d'una llera natural
    llista = sorted((fitxa(e) for e in estr), key=lambda f: (not f["urba"], f["T_desborda"] if f["T_desborda"] != ">500" else 999, f["capacitat_m3s"]))

    # comprovació: cabals 5.2-IC contra els de l'ACA als ponts
    comp = [{"pont": o["id"], "riu": o["riu"], "aca": o["qT"], "q52ic": o["q52ic"]} for o in aca if "q52ic" in o]
    ponts = [{"id": o["id"], "nom": o["nom"], "capacitat_m3s": o["q"], "supera_q500": o["q_supera_q500"],
              "Q10_100_500": [o["qT"][10], o["qT"][100], o["qT"][500]], "capacitat_aca": o["capacitat_aca"],
              "T_desborda": ">500" if o["q_supera_q500"] or not o["q"] else "<10" if o["q"] < o["qT"][10] else
              (lambda t: round(t) if t and t <= 500 else ">500")(t_desborda(o["q"], {10: o["qT"][10], 100: o["qT"][100], 500: o["qT"][500]})),
              "lonlat": [round(v, 5) for v in met_a_ll(*o["xy"])]} for o in aca]
    resum = {"km2": Cq[cid]["properties"]["km2"], "seccions": len(seccions),
             "seccions_met": sum(sc["font"] != mdt.nom for sc in seccions), "pluja": font_pluja,
             "p24_mm": {T: round(p24[T]) for T in (10, 100, 500)},
             "n_manning_mitja": round(float(np.mean([sc["n"] for sc in seccions])), 3),
             "llera_desborda_abans_T10": sum(1 for sc in seccions if (sc["T_llera"] or 501) < 10),
             "estrenyiments": len(llista), "estrenyiments_urbans": sum(f["urba"] for f in llista),
             "edificis_cadastre": n_edif, "seccions_sota_edificis": n_coberta,
             "seccions_lamina_aigua": sum(1 for sc in seccions if sc["Ww"] > 0),
             "seccions_mota": sum(1 for sc in seccions if sc.get("mota")),
             "seccions_corregides_lamina": sum(1 for sc in seccions if sc.get("d_sota_lamina")),
             "fondaria_sota_lamina_mediana_m": round(float(np.median([sc["d_sota_lamina"] for sc in seccions
                                                                      if sc.get("d_sota_lamina")] or [0])), 3)}
    print(json.dumps(resum, ensure_ascii=False))
    for f in llista[:15]:
        print(f"  {'U' if f['urba'] else ' '} T{str(f['T_desborda']):>4} {f['tipus'][:22]:22s} cap {f['capacitat_m3s']:7.1f}  Q10 {f['Q10']:6.1f} Q100 {f['Q100']:6.1f}"
              f"  A {f['conca_km2']:6.2f} km²  h {f['nivell_m']:.2f} W {f['amplada_m']:.0f}  {f['lonlat']}  {f['font_relleu'][:10]}")
    secs = [[round(sc["x"]), round(sc["y"]), round(sc["km2"], 2), round(sc["q_llera"], 1),
             (round(sc["T_llera"]) if sc["T_llera"] and sc["T_llera"] <= 500 else 501), round(sc["h"], 2), int(sc["conte"]),
             round(sc["W"], 1), round(sc["qT"][10], 1), round(sc["qT"][100], 1), int(sc["urba"])]
            for sc in seccions]
    return desa(cid, resum, llista, secs, comp, ponts, (met or mdt).crs.to_epsg())


def desa(cid, resum, llista, secs, comp=(), ponts=(), epsg=None):
    """cNNN.json (resum, estrenyiments i ponts: va al repositori i a la web) i les seccions a part, comprimides
    (SECCIONS/cNNN.json.gz: van a la Release estrenyiments-cat, que tot Catalunya fa ~100 MB)."""
    SORT.mkdir(parents=True, exist_ok=True); SECCIONS.mkdir(parents=True, exist_ok=True)
    (SORT / f"c{cid:03d}.json").write_text(json.dumps(
        {"resum": resum, "estrenyiments": llista, "ponts_aca": list(ponts), "comprovacio_cabals_aca": list(comp)},
        ensure_ascii=False, separators=(",", ":")))
    with gzip.open(SECCIONS / f"c{cid:03d}.json.gz", "wt") as f:
        json.dump({"epsg": epsg, "camps": ["x", "y", "conca_km2", "capacitat_m3s", "T_desborda", "nivell_m", "vall_conte",
                                           "amplada_m", "Q10", "Q100", "urba"], "files": secs}, f, separators=(",", ":"))
    return 0


def zones(cid, sortida):
    """Zones inundables de l'ACA (T10, T50, T100, T500) a la caixa de la conca (WFS, com prepara_zones.py)."""
    sys.path.insert(0, str(Path(__file__).parent))
    import prepara_zones as PZ
    Cq = {f["properties"]["id"]: f for f in json.loads(CONQUES.read_text())["features"]}
    bbox = list(shape(Cq[cid]["geometry"]).bounds)
    d = Path(sortida); d.mkdir(parents=True, exist_ok=True); n = 0
    for T in (10, 50, 100, 500):
        j = PZ.wfs(PZ.ACA, f"vishid:Zones_inundables_T{T}", bbox)
        if j is None:
            return 1
        (d / f"aca_T{T}.geojson").write_text(json.dumps(j, separators=(",", ":"))); n += len(j.get("features", []))
    print(f"c{cid:03d}: {n} polígons de zones inundables ACA")
    return 0 if n else 2


def valida(cid, dir_zones):
    """Comparació amb les zones inundables de l'ACA (T10, T50, T100, T500): a cada secció, el primer T la zona del qual
    surt clarament de la llera (cobreix ≥ 15 punts percentuals més del disc de radi amplada/2 + 30 m que la llera sola)."""
    from shapely.ops import transform as stransform
    from shapely.prepared import prep
    j = json.loads((SORT / f"c{cid:03d}.json").read_text())
    with gzip.open(SECCIONS / f"c{cid:03d}.json.gz", "rt") as f:
        sec = json.load(f)
    ep = sec["epsg"]; cp = sec["camps"]
    a = Transformer.from_crs(4326, ep, always_xy=True).transform
    zs = {}
    for T in (10, 50, 100, 500):
        f = Path(dir_zones) / f"aca_T{T}.geojson"
        if f.exists():
            from shapely.ops import unary_union
            zs[T] = unary_union([stransform(a, shape(g["geometry"])) for g in json.loads(f.read_text())["features"]]).buffer(0)
    taula, n = {}, 0
    if not zs or all(z.is_empty for z in zs.values()):
        print(f"c{cid:03d}: sense zones inundables ACA"); return 0
    for fila in sec["files"]:
        r = dict(zip(cp, fila))
        R = r["amplada_m"] / 2 + 30
        disc = Point(r["x"], r["y"]).buffer(R, 16)
        llera = min(r["amplada_m"] * 2 * R / disc.area, 1)
        Taca = None
        for T in sorted(zs):
            if zs[T].intersection(disc).area / disc.area >= llera + 0.15:
                Taca = T; break
        if not any(zs[T].intersects(disc) for T in zs):
            continue                                                    # fora dels trams estudiats per l'ACA
        n += 1
        Tm = r["T_desborda"]
        cm = "<10" if Tm < 10 else "10-50" if Tm < 50 else "50-100" if Tm < 100 else "100-500" if Tm <= 500 else ">500"
        ca = {10: "<10", 50: "10-50", 100: "50-100", 500: "100-500", None: ">500"}[Taca]
        taula.setdefault(ca, {}).setdefault(cm, 0); taula[ca][cm] += 1
    ordre = ["<10", "10-50", "50-100", "100-500", ">500"]
    print(f"{n} seccions dins dels trams de l'ACA. Files: ACA (zona que surt de la llera); columnes: aquest model")
    print("ACA \\ model " + " ".join(f"{c:>8}" for c in ordre))
    for ca in ordre:
        print(f"{ca:>11} " + " ".join(f"{taula.get(ca, {}).get(c, 0):>8}" for c in ordre))
    enc = sum(taula.get(c, {}).get(c, 0) for c in ordre)
    prop = sum(taula.get(ca, {}).get(cm, 0) for i, ca in enumerate(ordre) for k, cm in enumerate(ordre) if abs(i - k) <= 1)
    out = {"seccions": n, "taula": taula, "mateixa_classe": enc, "classe_veina_o_mateixa": prop}
    print(json.dumps(out, ensure_ascii=False))
    j["validacio_zones_aca"] = out
    (SORT / f"c{cid:03d}.json").write_text(json.dumps(j, ensure_ascii=False, separators=(",", ":")))

NOMS = {124: "Platja d'Aro (nucli)", 125: "Ridaura (Santa Cristina d'Aro, Castell d'Aro)", 188: "Congost (Granollers)",
        199: "Cabrera de Mar"}


def informe(ids, sortida):
    """Informe en Markdown amb els estrenyiments urbans i els ponts ACA de cada conca."""
    aca = {p["codi"]: p for p in json.loads(ACA.read_text())["pse"]} if ACA.exists() else {}
    o = ["# Estrenyiments de llera que desborden primer (opció)", "",
         "Generat per `scripts/estrenyiments.py informe` amb les dades de `data/catalunya/estrenyiments/`. El mètode és a "
         "la capçalera de l'script: seccions cada 20 m amb el MET 2 m de l'ICGC (LiDAR 2021-2023), nivell de desbordament "
         "quan l'aigua arriba a conreus, carrers o edificis (cobertes ICGC), Manning amb la n de les cobertes, cabals T2-T500 "
         "de la 5.2-IC amb la XEMA o els oficials de l'ACA on n'hi ha, i capacitat real dels ponts i passos coneguts. "
         "Res ajustat a mà.", ""]
    for cid in ids:
        f = SORT / f"c{cid:03d}.json"
        if not f.exists():
            continue
        j = json.loads(f.read_text()); r = j["resum"]
        if not r.get("seccions"):
            continue
        o += [f"## {NOMS.get(cid, 'Conca')} (c{cid:03d})", "",
              f"{r['seccions']} seccions; {r['llera_desborda_abans_T10']} amb la llera desbordant abans de T10. Pluja: {r['pluja']}. "
              f"{r['estrenyiments']} estrenyiments, {r['estrenyiments_urbans']} dins o a tocar de nucli.", ""]
        v = j.get("validacio_zones_aca")
        if v and v["seccions"]:
            ordre = ["<10", "10-50", "50-100", "100-500", ">500"]
            o += [f"Comparació amb les zones inundables de l'ACA: {v['seccions']} seccions; {round(100 * v['mateixa_classe'] / v['seccions'])} % "
                  f"a la mateixa classe de període i {round(100 * v['classe_veina_o_mateixa'] / v['seccions'])} % a la mateixa o la veïna.", "",
                  "| ACA \\ model | " + " | ".join(ordre) + " |", "|---|" + "---|" * 5]
            o += [f"| {a} | " + " | ".join(str(v["taula"].get(a, {}).get(m, 0)) for m in ordre) + " |" for a in ordre]
            o.append("")
        L = [e for e in j["estrenyiments"] if e["urba"]][:15]
        if L:
            o += ["| T desborda | Què | Riu (ACA) | Capacitat m³/s | Q10 | Q100 | Conca km² | Lloc |", "|---|---|---|---|---|---|---|---|"]
            for e in L:
                ob = e["obra"]; ancora = e["qT_font"][4:] if e["qT_font"].startswith("ACA") else None
                riu = aca.get(ancora, {}).get("riu", "") if ancora else ""
                que = f"llera, {e['llarg_m']:.0f} m, surt pel {e['surt_per'] or '—'}" if not ob else f"{ob['tipus']}: {ob['nom'][:50]}"
                T = e["T_desborda"]
                if ancora and e["capacitat_m3s"] < e["Q10"]:
                    T = "<10"                                            # l'ACA no dona cabals per sota de T10
                elif T == 2:
                    T = "≤2"
                la, lo = e["lonlat"][1], e["lonlat"][0]
                o.append(f"| {T} | {que} | {riu} | {e['capacitat_m3s']} | {e['Q10']} | {e['Q100']} | {e['conca_km2']} | "
                         f"[{la}, {lo}](https://www.openstreetmap.org/?mlat={la}&mlon={lo}#map=18/{la}/{lo}) |")
        P = sorted((p for p in j["ponts_aca"] if p["T_desborda"] != ">500"), key=lambda p: -1 if p["T_desborda"] == "<10" else p["T_desborda"])
        if P:
            o += ["", f"Ponts de l'inventari ACA que no passen el T500 ({len(P)}; capacitat de la fitxa oficial):", "",
                  "| T desborda | Pont | Capacitat m³/s | Q10 / Q100 / Q500 |", "|---|---|---|---|"]
            o += [f"| {p['T_desborda']} | {p['nom'][:70]} | {p['capacitat_m3s']} | {' / '.join(str(x) for x in p['Q10_100_500'])} |" for p in P[:20]]
        o.append("")
    if len(ids) > 10:
        # tot Catalunya: només el resum i els estrenyiments urbans que desborden abans
        o = o[:3] + [""]
        tots = []
        for cid in ids:
            f = SORT / f"c{cid:03d}.json"
            if f.exists():
                for e in json.loads(f.read_text())["estrenyiments"]:
                    if e["urba"]:
                        tots.append((cid, e))
        clau = lambda ce: ((-1 if ce[1]["qT_font"].startswith("ACA") and ce[1]["capacitat_m3s"] < ce[1]["Q10"] else
                            ce[1]["T_desborda"] if ce[1]["T_desborda"] != ">500" else 999), -ce[1]["conca_km2"])
        tots.sort(key=clau)
        # rius grans: el mètode hi agafa el canal de baix cabal i no la llera entre motes (capacitat massa baixa); tenen
        # estudi propi de l'ACA. Van a part, no a la llista principal (RIU_GRAN = 0 els hi torna).
        grans = [ce for ce in tots if RIU_GRAN and ce[1]["conca_km2"] > RIU_GRAN]
        tots = [ce for ce in tots if not (RIU_GRAN and ce[1]["conca_km2"] > RIU_GRAN)]
        val = [json.loads((SORT / f"c{c:03d}.json").read_text()).get("validacio_zones_aca") for c in ids if (SORT / f"c{c:03d}.json").exists()]
        val = [v for v in val if v and v.get("seccions")]
        ns = sum(v["seccions"] for v in val)
        o += [f"{len(ids)} conques; {len(tots) + len(grans)} estrenyiments dins o a tocar de nucli, {len(grans)} d'ells en "
              f"rius de més de {RIU_GRAN:.0f} km² (a part, al final).", ""]
        if ns:
            o += [f"Comparació amb les zones inundables de l'ACA a {len(val)} conques ({ns} seccions dins dels trams estudiats): "
                  f"{round(100 * sum(v['mateixa_classe'] for v in val) / ns)} % a la mateixa classe de període i "
                  f"{round(100 * sum(v['classe_veina_o_mateixa'] for v in val) / ns)} % a la mateixa o la veïna.", ""]
        o += ["Els 200 primers (els que desborden abans; a igualtat, els de més conca):", "",
              "| Conca | T desborda | Què | Capacitat m³/s | Q10 | Q100 | Conca km² | Lloc |", "|---|---|---|---|---|---|---|---|"]
        for cid, e in tots[:200]:
            ob = e["obra"]
            que = f"llera, {e['llarg_m']:.0f} m" if not ob else f"{ob['tipus']}: {ob['nom'][:45]}"
            T = "<10" if e["qT_font"].startswith("ACA") and e["capacitat_m3s"] < e["Q10"] else ("≤2" if e["T_desborda"] == 2 else e["T_desborda"])
            la, lo = e["lonlat"][1], e["lonlat"][0]
            o.append(f"| c{cid:03d} | {T} | {que} | {e['capacitat_m3s']} | {e['Q10']} | {e['Q100']} | {e['conca_km2']} | "
                     f"[{la}, {lo}](https://www.openstreetmap.org/?mlat={la}&mlon={lo}#map=18/{la}/{lo}) |")
        if grans:
            o += ["", f"## Rius de més de {RIU_GRAN:.0f} km² (estudi ACA)", "",
                  "Aquí el mètode agafa el canal de baix cabal i no la llera sencera entre motes, i la capacitat surt massa "
                  "baixa. Aquests rius tenen estudi propi de l'ACA (zones inundables T10-T500): cal mirar-ho amb aquell "
                  "estudi. Trams urbans per conca:", "", "| Conca | Trams | Conca km² màx. |", "|---|---|---|"]
            per = {}
            for cid, e in grans:
                n, k = per.get(cid, (0, 0)); per[cid] = (n + 1, max(k, e["conca_km2"]))
            for cid, (n, k) in sorted(per.items(), key=lambda t: -t[1][1]):
                o.append(f"| c{cid:03d} | {n} | {k:.0f} |")
    Path(sortida).write_text("\n".join(o) + "\n")
    print(f"{sortida}: {len(o)} línies")


def index():
    idx = {int(f.stem[1:]): json.loads(f.read_text())["resum"] for f in sorted(SORT.glob("c*.json"))}
    (SORT / "index.json").write_text(json.dumps({"conques": idx}, ensure_ascii=False, separators=(",", ":")))
    print(f"{len(idx)} conques")


def main():
    global SECCIONS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seccions", default=str(SECCIONS), help="carpeta de les seccions (cNNN.json.gz)")
    sub = ap.add_subparsers(dest="ordre", required=True)
    c = sub.add_parser("conca")
    c.add_argument("--id", type=int, required=True)
    c.add_argument("--mdt", required=True)
    c.add_argument("--met")
    c.add_argument("--font-met", default="MET 2 m ICGC (LiDAR 2021-2023)")
    c.add_argument("--cobertes")
    c.add_argument("--area-min", type=float, default=0.5)
    c.add_argument("--sense-edificis", dest="edificis", action="store_false", help="sense els edificis del Cadastre")
    c.add_argument("--sense-motes", dest="motes", action="store_false", help="sense la regla de les motes")
    z = sub.add_parser("zones")
    z.add_argument("--id", type=int, required=True)
    z.add_argument("--sortida", required=True)
    v = sub.add_parser("valida")
    v.add_argument("--id", type=int, required=True)
    v.add_argument("--zones", required=True, help="carpeta amb aca_T10/T50/T100/T500.geojson (prepara_zones.py)")
    i = sub.add_parser("informe")
    i.add_argument("--ids", default="125,188,199,124", help="ids separats per comes o «totes» (resum de tot Catalunya)")
    i.add_argument("--sortida", default="docs/estrenyiments.md")
    sub.add_parser("index")
    a = ap.parse_args()
    SECCIONS = Path(a.seccions)
    if a.ordre == "zones":
        return zones(a.id, a.sortida)
    if a.ordre == "informe":
        ids = (sorted(int(f.stem[1:]) for f in SORT.glob("c*.json")) if a.ids == "totes" else [int(x) for x in a.ids.split(",")])
        return informe(ids, a.sortida)
    return conca(a.id, a) if a.ordre == "conca" else valida(a.id, a.zones) if a.ordre == "valida" else index()


if __name__ == "__main__":
    sys.exit(main())
