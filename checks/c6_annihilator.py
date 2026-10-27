# Proposition 7.5 for m = 3 and Section 7.5: F_3^perp for C_6 (pure Python, exact).
# An independent check of the C_6 annihilators of cycles/annihilators.py.
# Usage: python3 c6_annihilator.py
from fractions import Fraction as Fr
n=6; full=(1<<n)-1
def mul(u,v):
    w={}
    for m1,c1 in u.items():
        for m2,c2 in v.items():
            if m1&m2: continue
            w[m1|m2]=w.get(m1|m2,0)+c1*c2
    return {k:c for k,c in w.items() if c}
Y=[{1<<((i+1)%n):Fr(1),1<<i:Fr(-1),(1<<i)|(1<<((i+1)%n)):Fr(-1)} for i in range(n)]
allp=[{0:Fr(1)}]; level=[{0:Fr(1)}]
for d in range(1,4):
    level=[mul(p,y) for p in level for y in Y]; allp+=level
piv={}
for v in allp:
    v=dict(v)
    for k in sorted(piv):
        if k in v:
            c=v[k]
            for kk,cc in piv[k].items(): v[kk]=v.get(kk,0)-c*cc
            v={a:b for a,b in v.items() if b}
    if v:
        k=min(v); c=v[k]; v={a:b/c for a,b in v.items()}
        for kk in list(piv):
            if k in piv[kk]:
                cc=piv[kk][k]; r={a:piv[kk].get(a,0)-cc*v.get(a,0) for a in set(piv[kk])|set(v)}; piv[kk]={a:b for a,b in r.items() if b}
        piv[k]=v
print('dim F_3 =',len(piv),' dim F_3^perp =',64-len(piv))
def var(i): return {1<<i:Fr(1)}
def lin(*ts):
    w={}
    for c,i in ts: w[1<<i]=w.get(1<<i,0)+Fr(c)
    return w
cands={'Omega_5':{full^(1<<i):Fr(1) for i in range(n)},'phi0+phi1':lin((1,0),(1,1)),'phi2-phi0':lin((1,2),(-1,0)),
       'phi4-phi0':lin((1,4),(-1,0)),'phi3-phi1':lin((1,3),(-1,1)),'phi5-phi1':lin((1,5),(-1,1))}
for name,g in cands.items():
    print(name,'annihilates F_3:',all(mul(g,v).get(full,0)==0 for v in piv.values()))
