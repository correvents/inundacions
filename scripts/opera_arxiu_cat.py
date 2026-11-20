#!/usr/bin/env python3
"""Arxiu OPERA retallat a Catalunya: un fitxer HDF5 per mes amb tots els compostos DBZH.

Només baixa el tros de cada compost europeu que cau sobre Catalunya (peticions HTTP Range sobre el
chunk HDF5 corresponent): ~200 kB per fitxer en lloc de ~2 MB.

Format de sortida (opera_cat_YYYY-MM.h5):
  /graella_N/dbz     uint8 (temps, y, x), gzip. dBZ = v/2 - 32; 0 = sense eco; 255 = sense cobertura
  /graella_N/temps   int64, segons UTC des de 1970 (final de l'interval del compost)
  atributs de /graella_N: projdef (proj4), xscale, yscale (m), x0, y0 (cantonada superior esquerra del retall,
  en coordenades de la projecció), font, codificacio.
Normalment hi ha una sola graella per mes (2 km/15 min fins 2023, 1 km/5 min des de 2024); si OPERA canvia
de graella dins el mes, n'hi ha més d'una.

Ús:  python3 scripts/opera_arxiu_cat.py --mes 2026-10 [--sortida dir] [--fils 16]"""
import argparse, io, re, sys, time, urllib.parse
from calendar import monthrange
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
import numpy as np, h5py, requests
from requests.adapters import HTTPAdapter
from pyproj import CRS, Transformer

S3 = "https://s3.waw3-1.cloudferro.com/openradar-archive"
CAT = (0.15, 40.50, 3.35, 42.90)     # lon/lat de Catalunya amb marge
SES = requests.Session()
SES.mount("https://", HTTPAdapter(pool_connections=32, pool_maxsize=32))
_TR = {}


def get(url, **kw):
    for i in range(5):
        try:
            r = SES.get(url, timeout=60, **kw); r.raise_for_status(); return r
        except Exception:
            if i == 4: raise
            time.sleep(2 * (i + 1))


class RangFile(io.RawIOBase):
    """Fitxer remot que només baixa els blocs que h5py llegeix."""
    def __init__(self, url, mida, bloc=16384):
        self.url, self.mida, self.bloc, self.pos, self.cache, self.bytes = url, mida, bloc, 0, {}, 0
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.mida + off; return self.pos
    def _blocs(self, a, b):
        falten = [i for i in range(a, b + 1) if i not in self.cache]
        while falten:
            i0 = i1 = falten.pop(0)
            while falten and falten[0] == i1 + 1: i1 = falten.pop(0)
            s, e = i0 * self.bloc, min((i1 + 1) * self.bloc, self.mida) - 1
            d = get(self.url, headers={"Range": f"bytes={s}-{e}"}).content
            self.bytes += len(d)
            for i in range(i0, i1 + 1): self.cache[i] = d[(i - i0) * self.bloc:(i - i0 + 1) * self.bloc]
    def readinto(self, buf):
        n = min(len(buf), self.mida - self.pos)
        if n <= 0: return 0
        a, b = self.pos // self.bloc, (self.pos + n - 1) // self.bloc
        self._blocs(a, b)
        buf[:n] = b"".join(self.cache[i] for i in range(a, b + 1))[self.pos - a * self.bloc:][:n]
        self.pos += n; return n


def llista(prefix):
    out, tok = [], None
    while True:
        p = {"list-type": 2, "prefix": prefix}
        if tok: p["continuation-token"] = tok
        x = get(S3, params=p).text
        out += [(k, int(s)) for k, s in re.findall(r"<Key>([^<]+)</Key>.*?<Size>(\d+)</Size>", x)]
        m = re.search(r"<NextContinuationToken>([^<]+)<", x)
        if not m: return out
        tok = m.group(1)


def txt(v): return v.decode() if isinstance(v, bytes) else str(v)


def retalla(k, mida):
    """(instant, uint8 retallat, signatura de graella) d'un compost."""
    t = datetime.strptime(re.search(r"OPERA@(\d{8}T\d{4})@", k).group(1), "%Y%m%dT%H%M").replace(tzinfo=timezone.utc)
    f = RangFile(S3 + "/" + urllib.parse.quote(k, safe="/@"), mida)
    with h5py.File(f, "r") as h:
        nom, wt = "dataset1/data1/data", None
        for d in [x for x in h if x.startswith("dataset")]:
            for cami in (f"{d}/data1/what", f"{d}/what"):
                if cami in h and txt(h[cami].attrs.get("quantity", b"")) == "DBZH":
                    nom, wt = f"{d}/data1/data", h[cami].attrs; break
            if wt is not None: break
        if wt is None: wt = h["dataset1/data1/what"].attrs
        wh = h["where"].attrs
        pd, xs, ys = txt(wh["projdef"]), float(wh["xscale"]), float(wh["yscale"])
        if pd not in _TR: _TR[pd] = Transformer.from_crs("EPSG:4326", CRS.from_proj4(pd), always_xy=True)
        tr = _TR[pd]
        xul, yul = tr.transform(float(wh["UL_lon"]), float(wh["UL_lat"]))
        LO, LA = np.meshgrid(np.linspace(CAT[0], CAT[2], 60), np.linspace(CAT[1], CAT[3], 60))
        X, Y = tr.transform(LO.ravel(), LA.ravel())
        c = np.floor((X - xul) / xs).astype(int); r = np.floor((yul - Y) / ys).astype(int)
        r0, r1, c0, c1 = r.min(), r.max() + 1, c.min(), c.max() + 1
        v = h[nom][r0:r1, c0:c1]
        g, o, nd, ud = (float(wt[x]) for x in ("gain", "offset", "nodata", "undetect"))
    q = np.clip(np.round((v * g + o + 32) * 2), 1, 254).astype(np.uint8)
    q[v == ud] = 0; q[v == nd] = 255
    sig = (pd, xs, ys, round(xul + c0 * xs, 3), round(yul - r0 * ys, 3), q.shape)
    return t, q, sig, f.bytes


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mes", required=True, help="YYYY-MM")
    ap.add_argument("--sortida", default=".")
    ap.add_argument("--fils", type=int, default=16)
    a = ap.parse_args()
    any_, mes = map(int, a.mes.split("-"))
    t0 = time.time()
    claus = []
    for d in range(1, monthrange(any_, mes)[1] + 1):
        claus += [(k, s) for k, s in llista(f"{any_}/{mes:02d}/{d:02d}/OPERA/COMP/") if re.search(r"@DBZH(_QIND)?\.h5$", k)]
    print(f"{a.mes}: {len(claus)} compostos DBZH a l'arxiu ({sum(s for _, s in claus)/1e9:.1f} GB Europa)", flush=True)
    if not claus:
        print("Cap compost: no es crea fitxer"); return 0
    res, errors, baixat = [], 0, 0
    with ThreadPoolExecutor(a.fils) as ex:
        futs = [ex.submit(retalla, k, s) for k, s in claus]
        for i, fu in enumerate(futs):
            try:
                t, q, sig, b = fu.result(); res.append((t, q, sig)); baixat += b
            except Exception as e:
                errors += 1; print(f"  error {claus[i][0]}: {type(e).__name__} {e}", flush=True)
            if (i + 1) % 500 == 0: print(f"  {i+1}/{len(futs)} ({time.time()-t0:.0f} s)", flush=True)
    grups = {}
    for t, q, sig in sorted(res, key=lambda x: x[0]):
        if t not in [x[0] for x in grups.get(sig, [])[-1:]]:
            grups.setdefault(sig, []).append((t, q))
    dest = Path(a.sortida); dest.mkdir(parents=True, exist_ok=True)
    fitxer = dest / f"opera_cat_{a.mes}.h5"
    with h5py.File(fitxer, "w") as h:
        h.attrs["font"] = "EUMETNET OPERA Open Radar Data (CC BY 4.0), compost DBZH, openradar-archive"
        h.attrs["caixa_lonlat"] = CAT
        for n, (sig, fil) in enumerate(sorted(grups.items(), key=lambda x: x[1][0][0])):
            pd, xs, ys, x0, y0, forma = sig
            g = h.create_group(f"graella_{n}")
            g.attrs.update({"projdef": pd, "xscale": xs, "yscale": ys, "x0": x0, "y0": y0,
                            "codificacio": "dBZ = v/2 - 32; 0 = sense eco (undetect); 255 = sense cobertura (nodata)"})
            cub = np.stack([q for _, q in fil])
            g.create_dataset("dbz", data=cub, chunks=(min(12, len(cub)), *forma), compression="gzip", compression_opts=6)
            g.create_dataset("temps", data=np.array([int(t.timestamp()) for t, _ in fil], dtype=np.int64))
            print(f"  graella_{n}: {cub.shape}, {xs:.0f} m, {fil[0][0]:%d %H:%M} - {fil[-1][0]:%d %H:%M}")
    print(f"{a.mes}: {len(res)} retallats, {errors} errors, {baixat/1e6:.0f} MB baixats, fitxer {fitxer.stat().st_size/1e6:.1f} MB, {time.time()-t0:.0f} s")
    return 1 if errors > 0.02 * len(claus) else 0


if __name__ == "__main__":
    sys.exit(main())
