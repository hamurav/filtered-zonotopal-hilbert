from __future__ import annotations
from itertools import combinations
from math import factorial
from fractions import Fraction
from time import monotonic
from collections import defaultdict

GRAPHS={
'G8_1':[(1,2),(2,3),(3,4),(4,1),(1,5),(5,6),(2,7),(3,8)],
'G8_2':[(1,2),(2,3),(3,4),(4,1),(1,5),(5,6),(2,7),(4,8)],
'T1_1':[(1,2),(2,3),(3,4),(4,5),(5,6),(2,7),(4,8),(2,9),(8,10)],
'T1_2':[(1,2),(2,3),(3,4),(4,5),(5,6),(2,7),(4,8),(2,9),(1,10)],
'T2_1':[(1,2),(2,3),(3,4),(4,5),(5,6),(2,7),(3,8),(8,9),(9,10)],
'T2_2':[(1,2),(2,3),(3,4),(4,5),(5,6),(2,7),(3,8),(1,9),(9,10)],
'K4':list(combinations(range(1,5),2)),
'C3double':[(1,2),(1,2),(2,3),(3,1)],
}

def generators(edges,coeff):
 vv=sorted(set(v for e in edges for v in e)); result=[]
 for v in vv:
  inc=[(1<<i,1 if v==min(e) else -1) for i,e in enumerate(edges) if v in e]
  g={}
  for k in range(1,len(inc)+1):
   a=coeff.get(k,0)*factorial(k)
   if a:
    for sub in combinations(inc,k):
     mask=0;sg=1
     for b,s in sub: mask|=b;sg*=s
     g[mask]=a*sg
  result.append(g)
 return result

def mul(v,w,p=0):
 out=defaultdict(int)
 for i,a in v.items():
  for j,b in w.items():
   if not i&j:out[i|j]+=a*b
 if p:return {i:a%p for i,a in out.items() if a%p}
 return {i:a for i,a in out.items() if a}

class Space:
 def __init__(self,p=0):self.p=p;self.rows={}
 def insert(self,v):
  v=v.copy();p=self.p
  while v:
   pivot=min(v);a=v[pivot];w=self.rows.get(pivot)
   if w is None:
    inv=pow(int(a),-1,p) if p else 1/a
    v={i:b*inv%p if p else b*inv for i,b in v.items()}
    self.rows[pivot]=v;return v
   for i,b in w.items():
    c=(v.get(i,0)-a*b)%p if p else v.get(i,0)-a*b
    if c:v[i]=c
    elif i in v:del v[i]
  return None

def hilbert(edges,coeff,p=0,verbose=False):
 coeff={i: (int(a.numerator)*pow(int(a.denominator),-1,p)%p if p else Fraction(a)) for i,a in coeff.items()}
 gg=generators(edges,coeff);sp=Space(p);sp.insert({0:1 if p else Fraction(1)})
 frontier=list(sp.rows.values());hh=[1]
 for d in range(1,len(edges)+1):
  new=[]
  for w in frontier:
   for g in gg:
    z=sp.insert(mul(w,g,p))
    if z:new.append(z)
  if not new:break
  hh.append(len(new));frontier=new
  if verbose:print(d,len(new),len(sp.rows),flush=True)
 return hh

if __name__=='__main__':
 for name in ['C3double','K4']:
  for a,b in [(0,0),(0,1),(1,0),(1,Fraction(2,3)),(1,1),(1,Fraction(4,3)),(1,2)]:
   start=monotonic();print(name,a,b,hilbert(GRAPHS[name],{1:1,2:a,3:b}),monotonic()-start,flush=True)
