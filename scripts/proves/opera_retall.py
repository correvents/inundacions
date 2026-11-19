#!/usr/bin/env python3
"""Prova 2: lectura per rangs HTTP (bytes reals) i volum retallat a Catalunya en uint8."""
import io, re, time, urllib.parse
import numpy as np, h5py, requests
from pyproj import CRS, Transformer

S3 = "https://s3.waw3-1.cloudferro.com/openradar-archive"
CAT = (0.15, 40.50, 3.35, 42.90)
SES = requests.Session()


class RangFile(io.RawIOBase):
    """Fitxer remot que només baixa els blocs que h5py demana (peticions HTTP Range)."""
    def __init__(self, url, bloc=16384):
        self.url, self.bloc, self.pos, self.cache, self.bytes, self.n = url, bloc, 0, {}, 0, 0
        self.mida = int(SES.head(url, timeout=60).headers["Content-Length"])
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.mida + off; return self.pos
    def _blocs(self, a, b):
        falten = [i for i in range(a, b + 1) if i not in self.cache]
        while falten:      # agrupa blocs consecutius en una sola petició
            i0 = i1 = falten.pop(0)
            while falten and falten[0] == i1 + 1: i1 = falten.pop(0)
            s, e = i0 * self.bloc, min((i1 + 1) * self.bloc, self.mida) - 1
            d = SES.get(self.url, headers={"Range": f"bytes={s}-{e}"}, timeout=60).content
            self.bytes += len(d); self.n += 1
            for i in range(i0, i1 + 1): self.cache[i] = d[(i - i0) * self.bloc:(i - i0 + 1) * self.bloc]
    def readinto(self, buf):
        n = min(len(buf), self.mida - self.pos)
        if n <= 0: return 0
        a, b = self.pos // self.bloc, (self.pos + n - 1) // self.bloc
        self._blocs(a, b)
        d = b"".join(self.cache[i] for i in range(a, b + 1))[self.pos - a * self.bloc:][:n]
        buf[:n] = d; self.pos += n; return n


def llista(prefix, delim=None):
    out, pre, tok = [], [], None
    while True:
        p = {"list-type": 2, "prefix": prefix}
        if delim: p["delimiter"] = delim
        if tok: p["continuation-token"] = tok
        x = SES.get(S3, params=p, timeout=60).text
        out += [(k, int(s)) for k, s in re.findall(r"<Key>([^<]+)</Key>.*?<Size>(\d+)</Size>", x)]
        pre += re.findall(r"<Prefix>([^<]+)</Prefix>", x)[1:]
        m = re.search(r"<NextContinuationToken>([^<]+)<", x)
        if not m: return out, pre
        tok = m.group(1)


def nom_dbzh(h):
    for d in [k for k in h if k.startswith("dataset")]:
        for cami in (f"{d}/data1/what", f"{d}/what"):
            if cami in h:
                q = h[cami].attrs.get("quantity", b"")
                if (q.decode() if isinstance(q, bytes) else str(q)) == "DBZH":
                    return f"{d}/data1/data", h[cami].attrs
    return "dataset1/data1/data", h["dataset1/data1/what"].attrs


def finestra(h):
    wh = h["where"].attrs
    pd = wh["projdef"]; pd = pd.decode() if isinstance(pd, bytes) else str(pd)
    tr = Transformer.from_crs("EPSG:4326", CRS.from_proj4(pd), always_xy=True)
    xul, yul = tr.transform(float(wh["UL_lon"]), float(wh["UL_lat"]))
    LO, LA = np.meshgrid(np.linspace(CAT[0], CAT[2], 50), np.linspace(CAT[1], CAT[3], 50))
    X, Y = tr.transform(LO.ravel(), LA.ravel())
    c = np.floor((X - xul) / float(wh["xscale"])).astype(int); r = np.floor((yul - Y) / float(wh["yscale"])).astype(int)
    return r.min(), r.max() + 1, c.min(), c.max() + 1


def a_uint8(v, wt):
    """dBZ quantitzat a 0,5 dB: 0 = sense eco, 255 = sense cobertura, 1..254 = -31,5..95 dBZ."""
    g, o, nd, ud = (float(wt[k]) for k in ("gain", "offset", "nodata", "undetect"))
    z = v * g + o
    q = np.clip(np.round((z + 32) * 2), 1, 254).astype(np.uint8)
    q[v == ud] = 0; q[v == nd] = 255
    return q


def dia(data, max_n=None):
    ks = sorted(k for k, _ in llista(data.replace("-", "/") + "/OPERA/COMP/")[0] if re.search(r"@DBZH(_QIND)?\.h5$", k))
    if max_n: ks = ks[:max_n]
    cub, tot_b, tot_n, mida_total, t0 = [], 0, 0, 0, time.time()
    for k in ks:
        f = RangFile(S3 + "/" + urllib.parse.quote(k, safe="/@"))
        with h5py.File(f, "r") as h:
            n, wt = nom_dbzh(h); r0, r1, c0, c1 = finestra(h)
            cub.append(a_uint8(h[n][r0:r1, c0:c1], wt))
        tot_b += f.bytes; tot_n += f.n; mida_total += f.mida
    cub = np.stack(cub)
    buf = io.BytesIO()
    with h5py.File(buf, "w") as h:
        h.create_dataset("dbz", data=cub, chunks=(min(12, len(cub)), *cub.shape[1:]), compression="gzip", compression_opts=6, shuffle=False)
    eco = float(((cub > 0) & (cub < 255) & (cub > 2 * (32 + 15))).mean())   # >15 dBZ
    print(f"{data}: {len(ks)} compostos, baixats {tot_b/1e6:.1f} MB de {mida_total/1e6:.0f} MB ({tot_b/mida_total*100:.1f}%), "
          f"{tot_n/len(ks):.1f} peticions/fitxer, {time.time()-t0:.0f} s; retall {cub.shape} uint8 gzip = {len(buf.getvalue())/1e6:.2f} MB/dia; "
          f"píxels >15 dBZ {eco*100:.1f}%")
    return len(buf.getvalue()), len(ks)


print("== Inici de l'arxiu")
anys = llista("", "/")[1]; print("anys:", anys)
for a in ("2011/", "2012/", "2013/"):
    mesos = llista(a, "/")[1]; print(a, mesos[:3], "...", len(mesos))
    if mesos:
        dies = llista(mesos[0], "/")[1]; print("  primer mes, dies:", dies[:3])

print("\n== Chunk de Catalunya dins un fitxer 2026 (offset i mida exactes)")
k = "2026/07/16/OPERA/COMP/OPERA@20260716T1200@0@DBZH.h5"
f = RangFile(S3 + "/" + urllib.parse.quote(k, safe="/@"))
with h5py.File(f, "r") as h:
    n, wt = nom_dbzh(h); ds = h[n]; r0, r1, c0, c1 = finestra(h)
    print("dtype", ds.dtype, "chunks", ds.chunks, "files", r0, r1, "cols", c0, c1)
    for i in range(ds.id.get_num_chunks()):
        ci = ds.id.get_chunk_info(i)
        if ci.chunk_offset[0] <= r0 < ci.chunk_offset[0] + ds.chunks[0] and ci.chunk_offset[1] <= c0 < ci.chunk_offset[1] + ds.chunks[1]:
            print("chunk", ci.chunk_offset, "byte_offset", ci.byte_offset, "mida", ci.size, "bytes")
    v = ds[r0:r1, c0:c1]
print(f"Lectura per rangs d'un fitxer: {f.n} peticions, {f.bytes/1e3:.0f} kB de {f.mida/1e3:.0f} kB")

print("\n== Dies sencers per rangs (retall uint8)")
mostres = {}
for d in ("2026-07-16", "2026-10-04", "2026-10-03"):
    mostres[d] = dia(d)
for d in ("2019-10-22", "2023-09-03", "2013-07-16"):
    mostres[d] = dia(d)
