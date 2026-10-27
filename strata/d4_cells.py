import sys, time, itertools, json
sys.path.insert(0, sys.argv[1])
from boundary_local import root_states, Tracker
from sympy.polys.fields import field
from sympy.polys.domains import QQ
D=4
def closure(coeff):
    gaps={1}; done={}; tr=Tracker()
    while True:
        new=set()
        for r in range(1,D):
            for dd in itertools.combinations_with_replacement(sorted(gaps),r):
                if dd in done: continue
                S=root_states(dd,coeff,tr); done[dd]=S
                new|={q-p for (p,q) in S}
        if new<=gaps: break
        gaps|=new
    return gaps, done, tr.factors
Fb,b,c=field('b,c',QQ)
Fx,x=field('x',QQ)
one=QQ(1)
cases={
 'generic (alpha=1; b,c free)': {1:Fb.one,2:Fb.one,3:b,4:c},
 'b=2 (c free)':               {1:Fx.one,2:Fx.one,3:2*Fx.one,4:x},
 'c=5b-5 (b free)':            {1:Fx.one,2:Fx.one,3:x,4:5*x-5},
 'c=b^2+b-1 (b free)':         {1:Fx.one,2:Fx.one,3:x,4:x**2+x-1},
 '(b,c)=(2,5)':                {1:one,2:one,3:QQ(2),4:QQ(5)},
 'alpha=0, beta=1 (c free)':   {1:Fx.one,2:Fx.zero,3:Fx.one,4:x},
 'alpha=0, beta=0, gamma=1':   {1:one,2:QQ(0),3:QQ(0),4:one},
 'graded f=u':                 {1:one},
}
res={}
for name,co in cases.items():
    t=time.monotonic(); g,done,fac=closure(co); res[name]=(g,done)
    print(f'{name}: gaps {sorted(g)}, {len(done)} tables, pivot factors {sorted(fac)}, {time.monotonic()-t:.0f}s',flush=True)
G=res['generic (alpha=1; b,c free)'][1]
for name,(g,done) in res.items():
    if name.startswith('generic'): continue
    diff=[dd for dd in done if dd not in G or done[dd]!=G[dd]]
    print(f'{name}: {len(diff)} tables differ from generic', diff[:12],flush=True)
