#!/usr/bin/env python3
"""Prova: es pot llegir només Catalunya d'un compost OPERA sense baixar tot Europa? I quant ocupa retallat?"""
import io, re, time, zlib, gzip, urllib.parse
from datetime import datetime, timedelta, timezone
import numpy as np, h5py, requests, fsspec
from pyproj import CRS, Transformer

S3 = "https://s3.waw3-1.cloudferro.com/openradar-archive"
S24 = "https://s3.waw3-1.cloudferro.com/openradar-24h"
CAT = (0.15, 40.50, 3.35, 42.90)   # lon/lat Catalunya


def llista(s3, prefix):
    out, tok = [], None
    while True:
        p = {"list-type": 2, "prefix": prefix}
        if tok: p["continuation-token"] = tok
        x = requests.get(s3, params=p, timeout=60).text
        out += [(k, int(s)) for k, s in re.findall(r"<Key>([^<]+)</Key>.*?<Size>(\d+)</Size>", x)]
        m = re.search(r"<NextContinuationToken>([^<]+)<", x)
        if not m: return out
        tok = m.group(1)


def finestra(h, nom):
    wh = h["where"].attrs
    pd = wh["projdef"]; pd = pd.decode() if isinstance(pd, bytes) else str(pd)
    tr = Transformer.from_crs("EPSG:4326", CRS.from_proj4(pd), always_xy=True)
    xul, yul = tr.transform(float(wh["UL_lon"]), float(wh["UL_lat"]))
    lo = np.linspace(CAT[0], CAT[2], 50); la = np.linspace(CAT[1], CAT[3], 50)
    LO, LA = np.meshgrid(lo, la); X, Y = tr.transform(LO.ravel(), LA.ravel())
    c = np.floor((X - xul) / float(wh["xscale"])).astype(int); r = np.floor((yul - Y) / float(wh["yscale"])).astype(int)
    return r.min(), r.max() + 1, c.min(), c.max() + 1, pd, float(wh["xscale"])


def nom_dbzh(h):
    for d in [k for k in h if k.startswith("dataset")]:
        for cami in (f"{d}/data1/what", f"{d}/what"):
            if cami in h and (lambda q: (q.decode() if isinstance(q, bytes) else str(q)))(h[cami].attrs.get("quantity", b"")) == "DBZH":
                return f"{d}/data1/data"
    return "dataset1/data1/data"


print("== 1. Mida dels fitxers a l'arxiu (un dia de mostra per any)")
for any_ in range(2011, 2027):
    dia = f"{any_}/07/16/OPERA/COMP/"
    ks = llista(S3, dia)
    dbzh = [(k, s) for k, s in ks if re.search(r"@DBZH(_QIND)?\.h5$", k)]
    tots = sum(s for _, s in ks)
    if dbzh:
        mida = np.mean([s for _, s in dbzh])
        print(f"{any_}: {len(ks)} fitxers COMP/dia ({tots/1e6:.0f} MB), DBZH {len(dbzh)}/dia, mitjana {mida/1e6:.2f} MB -> {mida*len(dbzh)*365/1e9:.0f} GB/any Europa")
    else:
        print(f"{any_}: cap DBZH ({len(ks)} fitxers)")
    if any_ in (2011, 2026):
        print("   productes:", sorted({re.sub(r'.*@0@', '', k) for k, _ in ks}))

print("\n== 2. Fitxer real 2026-07-16 12:00")
k = [k for k, _ in llista(S3, "2026/07/16/OPERA/COMP/OPERA@20260716T1200") if "@DBZH" in k][0]
url = S3 + "/" + urllib.parse.quote(k, safe="/@")
hd = requests.head(url, timeout=60); print(k, hd.headers.get("Content-Length"), "bytes; Accept-Ranges:", hd.headers.get("Accept-Ranges"))
r = requests.get(url, headers={"Range": "bytes=0-1023"}, timeout=60); print("Petició Range 0-1023 ->", r.status_code, len(r.content), "bytes")
t0 = time.time(); tot = requests.get(url, timeout=120).content; print(f"Baixada sencera: {len(tot)/1e6:.2f} MB en {time.time()-t0:.1f} s")
with h5py.File(io.BytesIO(tot), "r") as h:
    n = nom_dbzh(h); ds = h[n]
    print("dataset", n, ds.shape, ds.dtype, "chunks", ds.chunks, "compressió", ds.compression, ds.compression_opts)
    r0, r1, c0, c1, pd, xs = finestra(h, n)
    print(f"Finestra Catalunya: files {r0}:{r1}, columnes {c0}:{c1} -> {(r1-r0)}x{(c1-c0)} px de {xs:.0f} m ({(r1-r0)*(c1-c0)/ds.size*100:.2f}% d'Europa)")
    tros = ds[r0:r1, c0:c1]
    raw = tros.tobytes(); print(f"Retall cru {len(raw)/1e3:.0f} kB, zlib-6 {len(zlib.compress(raw,6))/1e3:.1f} kB")
    if ds.chunks:
        nch = ((r1-1)//ds.chunks[0] - r0//ds.chunks[0] + 1) * ((c1-1)//ds.chunks[1] - c0//ds.chunks[1] + 1)
        print("Chunks que toca la finestra:", nch, "de", int(np.ceil(ds.shape[0]/ds.chunks[0])*np.ceil(ds.shape[1]/ds.chunks[1])))

print("\n== 3. Lectura remota per rangs HTTP (només la finestra)")
fs = fsspec.filesystem("http", block_size=64 * 1024)
llegits = {"n": 0, "b": 0}
f = fs.open(url, "rb", cache_type="readahead", block_size=64*1024)
orig = f._fetch_range
def comptador(s, e):
    llegits["n"] += 1; llegits["b"] += e - s; return orig(s, e)
f._fetch_range = comptador
t0 = time.time()
with h5py.File(f, "r") as h:
    n = nom_dbzh(h); r0, r1, c0, c1, *_ = finestra(h, n)
    remot = h[n][r0:r1, c0:c1]
print(f"Remot: {llegits['n']} peticions, {llegits['b']/1e6:.2f} MB llegits de {len(tot)/1e6:.2f} MB, {time.time()-t0:.1f} s, igual que local: {np.array_equal(remot, tros)}")

print("\n== 4. Retall d'un dia sencer (288 compostos) en un sol HDF5 comprimit")
ks = sorted(k for k, _ in llista(S3, "2026/07/16/OPERA/COMP/") if re.search(r"@DBZH\.h5$", k))
cub = []
t0 = time.time()
for k in ks:
    b = requests.get(S3 + "/" + urllib.parse.quote(k, safe="/@"), timeout=120).content
    with h5py.File(io.BytesIO(b), "r") as h:
        cub.append(h[nom_dbzh(h)][r0:r1, c0:c1])
cub = np.stack(cub)
buf = io.BytesIO()
with h5py.File(buf, "w") as h:
    h.create_dataset("dbzh", data=cub, chunks=(12, cub.shape[1], cub.shape[2]), compression="gzip", compression_opts=6, shuffle=True)
print(f"{len(ks)} compostos en {time.time()-t0:.0f} s; cub {cub.shape} {cub.dtype}; HDF5 gzip: {len(buf.getvalue())/1e6:.2f} MB/dia (dia amb tempestes?)")
print("Fracció de píxels amb eco:", float(((cub != 0) & (cub != 255)).mean()) if cub.dtype == np.uint8 else "n/a")

print("\n== 5. Directe (openradar-24h): últims compostos")
ara = datetime.now(timezone.utc)
for dh in (0, 1):
    t = ara - timedelta(hours=dh)
    ks = llista(S24, t.strftime("%Y/%m/%d/OPERA/COMP/OPERA@%Y%m%dT%H"))
    print(t.strftime("%Y-%m-%dT%H"), len(ks), [k for k, _ in ks][-2:])
print("Arrel 24h:", re.findall(r"<Key>([^<]+)</Key>", requests.get(S24, params={"list-type": 2, "max-keys": 5}, timeout=60).text))
print("Prefixos 24h del dia:", re.findall(r"<Prefix>([^<]+)</Prefix>", requests.get(S24, params={"list-type": 2, "prefix": ara.strftime("%Y/%m/%d/"), "delimiter": "/"}, timeout=60).text))
