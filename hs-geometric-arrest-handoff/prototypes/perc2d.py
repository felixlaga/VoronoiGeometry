import numpy as np
def read2d(p):
    ls=open(p).readlines(); fr=[]; i=0
    while i<len(ls):
        n,L=ls[i].split(); n=int(n); L=float(L)
        b=np.array([[float(t) for t in ls[i+1+k].split()] for k in range(n)])
        fr.append((b[:,:2],b[:,2],L)); i+=n+1
    return fr
class UF:
    def __init__(s,n): s.p=list(range(n)); s.d=np.zeros((n,2)); s.w=np.zeros(2,bool)
    def find(s,a):
        dd=np.zeros(2)
        while s.p[a]!=a: dd+=s.d[a]; a=s.p[a]
        return a,dd
    def union(s,a,b,off):
        ra,da=s.find(a); rb,db=s.find(b)
        if ra==rb: s.w|=np.abs(da-db-off)>0.5
        else: s.p[ra]=rb; s.d[ra]=db-da+off
def perc(pos,rad,L,eps):
    n=len(pos); uf=UF(n); cut=2*rad.max()+2*eps
    nc=max(3,int(L/cut)); cs=L/nc
    idx=(pos/cs).astype(int)%nc; cell={}
    for i in range(n): cell.setdefault(tuple(idx[i]),[]).append(i)
    for i in range(n):
        cx,cy=idx[i]
        for a in(-1,0,1):
            for b in(-1,0,1):
                for j in cell.get(((cx+a)%nc,(cy+b)%nc),()):
                    if j<=i: continue
                    d=pos[i]-pos[j]; img=np.round(d/L); d=d-L*img
                    if d@d<(rad[i]+rad[j]+2*eps)**2: uf.union(i,j,img)
    return uf.w.all()
def eps_star(pos,rad,L):
    s=2*rad.mean(); lo,hi=1e-5,0.3
    if not perc(pos,rad,L,hi*s): return np.nan
    for _ in range(18):
        m=np.sqrt(lo*hi)
        if perc(pos,rad,L,m*s): hi=m
        else: lo=m
    return np.sqrt(lo*hi)
