#!/usr/bin/env python3
"""Ajust provisional del model de riu (llindar inicial + pèrdua contínua) amb els episodis 4/10/2026 i gener 2026, validació amb el catàleg XEMA+premsa i taula de llindars (config/ridaura_llindars_riu.json).
Ús:  python3 scripts/llindars_riu.py"""
import sys,json,numpy as np,pandas as pd
sys.path.insert(0,'scripts'); from model_cabal import simula, DT
M=json.load(open('config/ridaura_model_pluviometre.json')); QT=json.load(open('config/ridaura.json'))['cabals_T']['aforador']
s=pd.read_csv('data/ridaura/xema/DO.csv',parse_dates=['temps']).set_index('temps').mm; s=(s.resample(f'{DT}min').ffill()/ (30//DT))
def net(r,ia,phi):
    # llindar inicial ia (mm) i després pèrdua contínua phi (mm/h)
    out=np.zeros_like(r); acc=0
    for i,x in enumerate(r):
        if acc<ia: u=min(x,ia-acc); acc+=u; x-=u
        out[i]=max(x-phi*DT/60,0)
    return out
def pic(r,ia,phi,k=M['K'],m=M['m']): return simula(net(r,ia,phi),M['area_km2'],0.5,ia=0,c=1,k=k,m=m,n=2).max()
E1=s['2026-10-03T18:00':'2026-10-04T18:00'].values; E2=s['2026-01-17T18:00':'2026-01-21T00:00'].values
best=[]
for ia in range(0,161,10):
    for phi in np.arange(0,30.1,0.5):
        a,b=pic(E1,ia,phi),pic(E2,ia,phi)
        best.append((abs(np.log(a/303.7))+abs(np.log(b/28.6)),ia,phi,round(a),round(b,1)))
best.sort(); [print(x) for x in best[:12]]
json.dump(best[0],open('/tmp/claude-0/x/phi_best.json','w'))
print('--- amb K lliure i pluja de conca = estació ×1.14')
E1b=E1*1.14; E2b=E2*1.14
best=[]
for K in (100,200,300,400,563):
  for ia in range(60,161,10):
    for phi in np.arange(4,16.1,0.5):
        a,b=pic(E1b,ia,phi,k=K),pic(E2b,ia,phi,k=K)
        best.append((abs(np.log(a/303.7))+abs(np.log(b/28.6)),K,ia,phi,round(a),round(b,1)))
best.sort(); [print(x) for x in best[:8]]
print('--- validació amb tots els episodis (K 200, Ia 120, phi 12, ×1.14)')
K,IA,PHI=200,120,12
ep=pd.read_csv('data/ridaura/xema/episodis_pluja.csv');ep=ep[ep.estacio=='DO']
P=json.load(open('data/ridaura/model/premsa_episodis.json'))['episodis']; lab={p['inici']:p for p in P}
T=lambda q: next((k for k in ('Q500','Q100','Q50','Q10','QMCO') if q>=QT[k]),'-')
res=[]
for r in ep.itertuples():
    rr=s[pd.Timestamp(r.inici_utc)-pd.Timedelta(hours=1):pd.Timestamp(r.fi_utc)+pd.Timedelta(hours=1)].values*1.14
    q=pic(rr,IA,PHI,k=K); p=lab.get(r.inici_utc[:16],{})
    res.append((round(q,1),r.inici_utc[:10],T(q),p.get('inundacio'),p.get('ridaura_desborda'),r.mm,r.max3h))
res.sort(reverse=True); [print(x) for x in res[:15]]
print('Pluja necessària (mm, mitjana de conca) per arribar a Q10/Q50/Q100/Q500 amb pluja uniforme')
for d in (1,2,3,6,12,24):
    from scipy.optimize import brentq
    f=lambda P,k: pic(np.r_[np.full(int(d*60/DT),P/int(d*60/DT)),np.zeros(288)],IA,PHI,k=K)-QT[k]
    print(f'{d:>2} h',[round(brentq(lambda P: f(P,k),IA+PHI*d+0.1,3000)) for k in ('Q10','Q50','Q100','Q500')])
