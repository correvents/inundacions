import rasterio, numpy as np, sys
from rasterio.merge import merge
d=sys.argv[1]
srcs=[rasterio.open(f'{d}/E002.tif'),rasterio.open(f'{d}/E003.tif')]
arr,tr=merge(srcs,bounds=(2.80,41.72,3.12,41.95))
p=srcs[0].profile; p.update(height=arr.shape[1],width=arr.shape[2],transform=tr,driver='GTiff',compress='deflate',tiled=False)
p.pop('blockxsize',None);p.pop('blockysize',None)
with rasterio.open(f'{d}/crop.tif','w',**p) as o: o.write(arr)
print(arr.shape, arr.min(), arr.max())
