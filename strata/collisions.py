# Universal Hilbert collisions among trees, using complete parameter partitions:
# D=3: four cells (paper, Theorem 6.7); D=4: eight cells (strata_probe, exploratory).
import sys, time, itertools
sys.path.insert(0, sys.argv[1]); sys.path.insert(0,'.')
from fractions import Fraction as Fr
from collections import defaultdict
from boundary_local import root_states
import importlib.util
spec=importlib.util.spec_from_file_location('w','d4_witness_lib.py'); W=importlib.util.module_from_spec(spec); spec.loader.exec_module(W)
D=int(sys.argv[2]); NMAX=int(sys.argv[3])
PTS={3:[{1:1},{1:1,3:1},{1:1,2:1,3:2},{1:1,2:1}],
     4:[{1:1},{1:1,4:1},{1:1,3:1},{1:1,3:1,4:1},{1:1,2:1,3:3,4:7},{1:1,2:1,3:2,4:7},{1:1,2:1,3:3,4:10},{1:1,2:1,3:2,4:5}]}[D]
GEN=PTS[-1] if D==3 else PTS[4]
locs=[]
for co in PTS:
    cf={k:Fr(v) for k,v in co.items()}; cache={}
    def L(dd,cf=cf,cache=cache):
        if dd not in cache: cache[dd]=root_states(dd,cf)
        return cache[dd]
    locs.append(L)
gi=PTS.index(GEN)
byn=defaultdict(list); t=time.monotonic()
for n,e in W.trees(NMAX,D=D):
    byn[n].append((e,tuple(W.hilbert(e,L) for L in locs)))
print(f'D={D} enumeration+Hilbert {time.monotonic()-t:.0f}s')
print('n  trees  profile-classes  trees-in-collisions  largest-class | generic-f classes  trees-not-separated-by-generic-f')
examples={}
for n in sorted(byn):
    prof=defaultdict(list); gen=defaultdict(list)
    for e,H in byn[n]:
        prof[H].append(e); gen[H[gi]].append(e)
    coll=[c for c in prof.values() if len(c)>1]
    print(f'{n:2d} {len(byn[n]):6d} {len(prof):8d} {sum(len(c) for c in coll):10d} {max(len(c) for c in prof.values()):6d} | {len(gen):8d} {sum(len(c) for c in gen.values() if len(c)>1):8d}')
    if coll and n not in examples: examples[n]=coll[:3]
for n in sorted(examples)[:3]:
    for c in examples[n]: print('n',n,'class:',c)
