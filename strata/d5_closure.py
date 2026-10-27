import sys, time, itertools
sys.path.insert(0, sys.argv[1])
from boundary_local import root_states, Tracker
from sympy.polys.fields import field
from sympy.polys.domains import QQ
D=5
F,b,c,d=field('b,c,d',QQ)
coeff={1:F.one,2:F.one,3:b,4:c,5:d}
gaps={1}; done={}; tr=Tracker(); t0=time.monotonic()
while True:
    new=set()
    for r in range(1,D):
        for dd in itertools.combinations_with_replacement(sorted(gaps),r):
            if dd in done: continue
            s=time.monotonic(); S=root_states(dd,coeff,tr); done[dd]=S
            print(dd, f'{time.monotonic()-s:.1f}s', 'factors so far', sorted(tr.factors), flush=True)
            new|={q-p for (p,q) in S}
    if new<=gaps: break
    gaps|=new
print('DONE gap types',sorted(gaps),'local tables',len(done),'factors',sorted(tr.factors),f'total {time.monotonic()-t0:.0f}s',flush=True)
