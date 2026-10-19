import sys, json, numpy as np
np.in1d = np.isin
from pysheds.grid import Grid
from shapely.geometry import shape, mapping, LineString
from shapely.ops import unary_union
from pyproj import Geod
W=sys.argv[1]; OUT=sys.argv[2]
g=Grid.from_raster(f'{W}/dem/crop.tif'); dem=g.read_raster(f'{W}/dem/crop.tif')
dem=g.resolve_flats(g.fill_depressions(g.fill_pits(dem)))
fd=g.flowdir(dem); acc=g.accumulation(fd)
geod=Geod(ellps='WGS84')
cell_km2=(30.9*np.cos(np.radians(41.83))*30.9)/1e6  # ~1 arcsec
def basin(lat,lon,name,thr=2000):
    x,y=g.snap_to_mask(acc>thr,(lon,lat))
    c=g.catchment(x=x,y=y,fdir=fd,xytype='coordinate')
    shapes=[shape(s) for s,v in g.polygonize(c.astype(np.uint8)) if v==1]
    poly=unary_union(shapes)
    area=abs(geod.geometry_area_perimeter(poly)[0])/1e6
    # longest flow path via distance_to_outlet
    dist=g.distance_to_outlet(x=x,y=y,fdir=fd,xytype='coordinate')
    d=np.where(np.isfinite(dist)&c,dist,np.nan)
    L_cells=np.nanmax(d); L_km=L_cells*0.0278  # approx mean cell step km
    i=np.unravel_index(np.nanargmax(d),d.shape)
    zmax=float(dem[i]); r_,c_=g.nearest_cell(x,y); zmin=float(dem[c_,r_]) if False else float(np.asarray(dem)[~np.isnan(d)].min())
    J=(zmax-zmin)/(L_km*1000)
    tc=0.3*(L_km/J**0.25)**0.76
    print(f'{name}: outlet {y:.5f},{x:.5f} area {area:.1f} km2, L {L_km:.1f} km, z {zmin:.0f}-{zmax:.0f} m, J {J:.4f}, Tc(Témez) {tc:.1f} h')
    return poly,dict(nom=name,area_km2=round(area,1),L_km=round(L_km,1),zmax=round(zmax),zmin=round(zmin),pendent=round(J,4),tc_temez_h=round(tc,2)),c
pm,propm,cm=basin(41.8051,3.0671,thr=90000,name='Ridaura a la desembocadura (Platja d\'Aro)')
pg,propg,cg=basin(41.81814,2.97912,'Ridaura a l\'aforador (Santa Cristina d\'Aro)')
# streams > 0.5 km2 inside basin
br=g.extract_river_network(fd,acc>int(0.5/cell_km2)) 
feats=[]
for f in br['features']:
    ls=shape(f['geometry'])
    if ls.intersects(pm): feats.append({'type':'Feature','properties':{},'geometry':mapping(ls.intersection(pm.buffer(0.0005)))})
fc={'type':'FeatureCollection','features':[
 {'type':'Feature','properties':propm,'geometry':mapping(pm.simplify(0.0002))},
 {'type':'Feature','properties':propg,'geometry':mapping(pg.simplify(0.0002))}]+feats}
json.dump(fc,open(f'{OUT}/conca_ridaura.geojson','w'))
print('trams xarxa:',len(feats))
