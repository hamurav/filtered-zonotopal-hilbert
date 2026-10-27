from __future__ import annotations
from itertools import product, combinations
from collections import Counter, defaultdict
from math import factorial
from sympy.polys.fields import field, FracElement
from sympy.polys.domains import QQ
import sympy as sp


def mul(v,w):
 out={}
 for i,a in v.items():
  for j,b in w.items():
   if not i&j:
    k=i|j
    out[k]=out.get(k,0)+a*b
 return {k:v for k,v in out.items() if v}

class Space:
 def __init__(self, track=None):
  self.rows={}
  self.track=track
 def insert(self,v):
  v=v.copy()
  while v:
   p=min(v); a=v[p]
   if p not in self.rows:
    if self.track is not None:
     self.track(a)
    inv=1/a
    v={j:b*inv for j,b in v.items()}
    self.rows[p]=v
    return v
   w=self.rows[p]
   for j,b in w.items():
    z=v.get(j,0)-a*b
    if z: v[j]=z
    elif j in v: del v[j]
  return None
 def extend(self,seq):
  for v in seq: self.insert(v)
 def rank(self): return len(self.rows)

def root_generators(r, coeff):
 # variables 0..r-1 inputs; variable r is outgoing boundary.
 lin={1<<i:-coeff[1] for i in range(r)}; lin[1<<r]=coeff[1]
 y={}; pw={0:coeff[1]}
 for k in range(1,r+2):
  pw=mul(pw,lin)
  for m,a in pw.items(): y[m]=y.get(m,0)+coeff.get(k,0)*a
 return y

def root_states(deltas,coeff,track=None,verbose=False):
 r=len(deltas); tbit=1<<r
 y=root_generators(r,coeff)
 pows=[{0:coeff[1]}]
 for k in range(1,r+2): pows.append(mul(pows[-1],y))
 bycost=defaultdict(list)
 offset=-sum(min(d,0) for d in deltas)
 for sub in range(1<<r):
  cost=sum(deltas[i] for i in range(r) if sub>>i&1)
  for k,w in enumerate(pows):
   z=mul({sub:coeff[1]},w)
   if z: bycost[cost+k+offset].append(z)
 maxcost=max(bycost)
 ps=Space(track); Q=[]; R=[]; pdims=[]
 for i in range(maxcost+1):
  ps.extend(bycost[i]); pdims.append(ps.rank())
  qs=Space(track); rs=Space(track)
  for p,w in ps.rows.items():
   if p<tbit: qs.insert({k:v for k,v in w.items() if k<tbit})
   else: rs.insert({k^tbit:v for k,v in w.items()})
  Q.append(qs);R.append(rs)
 if ps.rank() != 2*tbit:
  print('WARNING NONEXHAUSTIVE',deltas,pdims)
 dims={(-1,j):0 for j in range(-1,maxcost+1)}
 dims.update({(i,-1):0 for i in range(maxcost+1)})
 for i in range(maxcost+1):
  for j in range(maxcost+1):
   both=Space(track);both.extend(Q[i].rows.values());both.extend(R[j].rows.values())
   dims[i,j]=Q[i].rank()+R[j].rank()-both.rank()
 out=Counter()
 for i in range(maxcost+1):
  for j in range(maxcost+1):
   c=dims[i,j]-dims[i-1,j]-dims[i,j-1]+dims[i-1,j-1]
   assert c>=0,(deltas,i,j,c)
   if c:out[i-offset,j-offset]=c
 assert sum(out.values())==tbit
 return out

def shift(S,a):return Counter({(p+a,q+a):c for (p,q),c in S.items()})
def add(*Ss):
 C=Counter()
 for S in Ss:C.update(S)
 return C

def scale(S,c):return Counter({k:v*c for k,v in S.items()})
def trans(S,quadratic=True):
 C=Counter()
 for (p,q),w in S.items():
  if q<=p: pairs=[(p,p+1),(q,q+1)]
  elif q==p+1:
   pairs=[(p,p+2),(p+1,p+1)] if quadratic else [(p,p+1),(p+1,p+2)]
  else: pairs=[(p,q),(p+1,p+2)]
  for z in pairs:C[z]+=w
 return C

def close(S,T):
 H=Counter()
 for (p,q),a in S.items():
  for (u,v),b in T.items():
   H[p+u]+=a*b
   H[min(q+u,p+v)]+=a*b
 return [H[i] for i in range(max(H)+1)]

class Tracker:
 def __init__(self):self.factors=set()
 def __call__(self,a):
  # only elements of a rational function field have pivot factors; with
  # SymPy's flint ground types, plain rationals also have a numer() method
  if isinstance(a,FracElement):
   for p in [a.numer,a.denom]:
    if not p.is_ground:
     for f,e in p.factor_list()[1]:
      self.factors.add(str(f.monic()))

if __name__=='__main__':
 import sys,time
 F,b,c=field('b,c',QQ)
 coeff={1:F.one,2:F.one,3:b,4:c}
 tracker=Tracker()
 R={}
 for dd in [(1,1),(1,2),(2,2),(0,2),(0,0),(0,1),(1,1,1),(1,1,2),(0,1,1)]:
  now=time.monotonic(); R[dd]=root_states(dd,coeff,tracker)
  print(dd,sorted(R[dd].items()),'factors=',sorted(tracker.factors),'seconds',time.monotonic()-now,flush=True)
 A=R[1,1];B=add(R[2,2],scale(shift(R[1,2],1),6),scale(shift(R[1,1],2),9));C=add(R[1,2],scale(shift(R[1,1],1),3))
 print('T2',close(A,B),close(C,C),flush=True)
 A1=R[1,1,1];B1=add(R[2,2],scale(shift(R[0,2],1),2),shift(R[0,0],2))
 A2=add(R[1,1,2],shift(R[0,1,1],1));B2=add(R[1,2],shift(R[0,1],1))
 print('T1',close(trans(A1),B1),close(trans(A2),B2),flush=True)
