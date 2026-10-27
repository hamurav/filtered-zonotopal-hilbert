# Four-point (complete, by Theorem 6.7) profiles of the n=10 subcubic collision pairs
# T2 and T3 (Proposition 6.15, Tables 7 and 8), by exact boundary states.
# Usage: python3 t3_profile.py ../certificates ../strata/d4_witness_lib.py
import sys, importlib.util
sys.path.insert(0, sys.argv[1])
from fractions import Fraction as Fr
from boundary_local import root_states
spec=importlib.util.spec_from_file_location('w',sys.argv[2]); W=importlib.util.module_from_spec(spec); spec.loader.exec_module(W)
PTS=[('u',{1:1}),('u+u^3',{1:1,3:1}),('u+u^2+2u^3',{1:1,2:1,3:2}),('u+u^2',{1:1,2:1})]
pairs={'T2':([(0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(3,7),(6,8),(7,9)],[(0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(4,7),(6,8),(7,9)]),
       'T3':([(0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(2,7),(3,8),(8,9)],[(0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(2,7),(4,8),(8,9)])}
for name,(a,b) in pairs.items():
    for lab,co in PTS:
        cf={k:Fr(v) for k,v in co.items()}; cache={}
        def L(dd):
            if dd not in cache: cache[dd]=root_states(dd,cf)
            return cache[dd]
        ha,hb=W.hilbert(a,L),W.hilbert(b,L)
        print(name,lab,ha,hb,'EQUAL' if ha==hb else 'DIFF')
