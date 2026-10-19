import json, sys, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from shapely.geometry import shape
fc=json.load(open(sys.argv[1]))
fig,ax=plt.subplots(figsize=(9,5.5))
cols=['#1f77b4','#ff7f0e']
for i,f in enumerate(fc['features'][:2]):
    g=shape(f['geometry']); gs=getattr(g,'geoms',[g])
    for p in gs: x,y=p.exterior.xy; ax.fill(x,y,alpha=.18 if i==0 else .35,color=cols[i],label=f"{f['properties']['nom']} · {f['properties']['area_km2']} km²" if p is gs[0] else None)
for f in fc['features'][2:]:
    g=shape(f['geometry'])
    for l in getattr(g,'geoms',[g]):
        if l.geom_type=='LineString': x,y=l.xy; ax.plot(x,y,color='#08519c',lw=.8)
ax.plot(2.97912,41.81814,'k^',ms=8,label='Aforador ACA 171812-001'); ax.plot(3.0672,41.805,'ro',ms=7,label='Desembocadura (Platja d\'Aro)')
ax.set_aspect(1/0.745); ax.legend(loc='lower left',fontsize=8); ax.set_title('Conca del Ridaura (MDE Copernicus 30 m)'); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig(sys.argv[2],dpi=130)
