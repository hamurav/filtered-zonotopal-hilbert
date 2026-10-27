from boundary_local import Tracker,mul
class Space:
 def __init__(self,track):self.rows={};self.track=track
 def rank(self):return len(self.rows)
 def insert(self,v):
  v=v.copy()
  while v:
   p=min(v,key=lambda k:(k.bit_count(),k));a=v[p]
   if p not in self.rows:
    self.track(a);inv=1/a;v={k:z*inv for k,z in v.items()};self.rows[p]=v;return v
   for k,z in self.rows[p].items():
    c=v.get(k,0)-a*z
    if c:v[k]=c
    elif k in v:del v[k]
  return None
from direct_hilbert import generators,GRAPHS
from sympy.polys.fields import field
from sympy.polys.domains import QQ
from time import monotonic
import sys,json
F,z=field('z',QQ)

def hilbert_symbolic(edges,coef,verbose=True):
 tr=Tracker();space=Space(tr);gg=generators(edges,coef)
 space.insert({0:F.one});front=list(space.rows.values());hs=[1]
 for i in range(1,len(edges)+1):
  new=[];st=monotonic()
  for v in front:
   for g in gg:
    red=space.insert(mul(v,g))
    if red is not None:new.append(red)
  if not new:break
  hs.append(len(new));front=new
  if verbose:print('degree',i,'h',len(new),'cum',space.rank(),'bad',sorted(tr.factors),'seconds',monotonic()-st,flush=True)
 return hs,sorted(tr.factors)
if __name__=='__main__':
 for name in sys.argv[1:] or ['G8_1','G8_2']:
  print('GRAPH',name,flush=True)
  res=hilbert_symbolic(GRAPHS[name],{1:F.one,2:F.one,3:z})
  print('FINAL',name,res,flush=True)
