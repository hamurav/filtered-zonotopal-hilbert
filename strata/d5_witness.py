import sys, itertools, time
sys.path.insert(0, sys.argv[1]); sys.path.insert(0,'.')
from fractions import Fraction as Fr
from boundary_local import root_states
import importlib.util
spec=importlib.util.spec_from_file_location('w','d4_witness_lib.py'); W=importlib.util.module_from_spec(spec); spec.loader.exec_module(W)
NMAX=int(sys.argv[2])
b,c=Fr(7,3),Fr(11,2)
pts={'generic':(b,c,Fr(13,7)),
 '[u3]f^-1=0':(Fr(2),c,Fr(13,7)),
 '[u4]f^-1=0':(b,5*b-5,Fr(13,7)),
 '[u5]f^-1=0':(b,c,3*b*b-21*b+6*c+14),
 'S5: d=3bc-2b^3':(b,c,3*b*c-2*b**3),
 'S8':(b,c,(3*b**3-2*b**2-4*b*c+6*b+c*c-2*c-3)/(b-2))}
locs={}
for k,(B,C,D) in pts.items():
    cf={1:Fr(1),2:Fr(1),3:B,4:C,5:D}; cache={}
    def L(dd,cf=cf,cache=cache):
        if dd not in cache: cache[dd]=root_states(dd,cf)
        return cache[dd]
    locs[k]=L
first={}; t=time.monotonic(); cnt=0
for n,e in W.trees(NMAX,D=5):
    cnt+=1
    H={k:W.hilbert(e,locs[k]) for k in pts}
    for k in pts:
        if k!='generic' and H[k]!=H['generic'] and k not in first: first[k]=(n,e,H['generic'],H[k])
print(cnt,'trees',f'{time.monotonic()-t:.0f}s')
for k in pts:
    if k=='generic': continue
    print(k, first.get(k,'NOT SEEN'))
