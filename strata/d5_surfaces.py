# Check each candidate D=5 hypersurface (alpha=1): recompute all local tables over its
# function field Q(x,y) and compare with the tables at a generic rational point.
import sys, time, itertools, pickle
sys.path.insert(0, sys.argv[1])
from boundary_local import root_states, Tracker
from sympy.polys.fields import field
from sympy.polys.domains import QQ
from fractions import Fraction as Fr
import sympy as sp
D=5
FACTORS=['b - 2','c - 5*b + 5','d - 3*b**2 + 21*b - 6*c - 14','b - 3','b**2 + b - c - 1',
 'b**3 - 3*b*c/2 + d/2',
 'b**3 - b**2*c/5 + b**2 - 17*b*c/10 + b + c**2/5 - 4*c/5 + d/2 - 1',
 'b**3 - b**2*c/10 - 23*b*c/15 - b*d/6 + 3*b + 4*c**2/15 - 6*c/5 + 2*d/3 - 2',
 'b**3 - 2*b**2/3 - 4*b*c/3 - b*d/3 + 2*b + c**2/3 - 2*c/3 + 2*d/3 - 1']
B,C,Dd,X,Y=sp.symbols('b c d x y')
def closure(coeff,tr=None):
    gaps={1}; done={}
    while True:
        new=set()
        for r in range(1,D):
            for dd in itertools.combinations_with_replacement(sorted(gaps),r):
                if dd in done: continue
                S=root_states(dd,coeff,tr); done[dd]=S
                new|={q-p for (p,q) in S}
        if new<=gaps: break
        gaps|=new
    return gaps,done
def param(fs):
    e=sp.sympify(fs,locals={'b':B,'c':C,'d':Dd})
    for v,(p,q) in [(Dd,(B,C)),(C,(B,Dd)),(B,(C,Dd))]:
        if sp.degree(e,v)==1:
            sol=sp.solve(e,v)[0]
            sub={p:X,q:Y}
            vals={B:B,C:C,Dd:Dd}; vals.update(sub); vals[v]=sol.subs(sub)
            return [vals[B],vals[C],vals[Dd]]
    raise ValueError(fs)
if __name__=='__main__':
    mode=sys.argv[2]
    if mode=='generic':
        pt={1:Fr(1),2:Fr(1),3:Fr(7,3),4:Fr(11,2),5:Fr(13,7)}
        e=[sp.sympify(f,locals={'b':B,'c':C,'d':Dd}).subs({B:sp.Rational(7,3),C:sp.Rational(11,2),Dd:sp.Rational(13,7)}) for f in FACTORS]
        assert all(v!=0 for v in e), e
        t=time.monotonic(); g,done=closure(pt)
        pickle.dump(done,open('d5_generic_tables.pkl','wb'))
        print('generic point tables',len(done),sorted(g),f'{time.monotonic()-t:.0f}s')
    else:
        i=int(mode); fs=FACTORS[i]; vb,vc,vd=param(fs)
        F,x,y=field('x,y',QQ)
        co={1:F.one,2:F.one,3:F.from_expr(sp.sympify(vb)) if vb.free_symbols else F.one*QQ(str(vb)),
            4:F.from_expr(vc) if vc.free_symbols else F.one*QQ(str(vc)),
            5:F.from_expr(vd) if vd.free_symbols else F.one*QQ(str(vd))}
        tr=Tracker(); t=time.monotonic(); g,done=closure(co,tr)
        G=pickle.load(open('d5_generic_tables.pkl','rb'))
        diff=[dd for dd in done if dd not in G or done[dd]!=G[dd]]
        print(f'[{i}] {fs}: param b,c,d = {vb}, {vc}, {vd}; gaps {sorted(g)}; {len(done)} tables; {len(diff)} differ from generic {diff[:8]}; subfactors {sorted(tr.factors)}; {time.monotonic()-t:.0f}s',flush=True)
