"""Complete Hilbert profiles of the subcubic graphs in tables B1-B3 (Delta<=3).
Normalized parameter space: alpha=0 gives f=u or u+u^3 (after scaling); alpha!=0 gives
f=u+u^2+lam*u^3 (lam=beta/alpha^2).  Generic ranks over Q(lam) with recorded pivot
factors; every factor is recomputed over its residue field."""
import sys, re, time
sys.path.insert(0,'../certificates')
from fractions import Fraction as Fr
from math import factorial
from itertools import combinations
from collections import defaultdict
import sympy as sp
from sympy.polys.fields import field
from sympy.polys.domains import QQ
from direct_hilbert import hilbert, mul

def fmul(v,w,zero):
    out={}
    for i,a in v.items():
        for j,b in w.items():
            if not i&j:
                k=i|j; out[k]=out.get(k,zero)+a*b
    return {k:x for k,x in out.items() if x!=zero}
class TSpace:
    def __init__(s,track,one): s.rows={}; s.track=track; s.one=one; s.zero=one-one
    def insert(s,v):
        v=dict(v)
        while v:
            p=min(v); a=v[p]; w=s.rows.get(p)
            if w is None:
                s.track(a); inv=s.one/a
                v={i:b*inv for i,b in v.items()}; s.rows[p]=v; return v
            for i,b in w.items():
                c=v.get(i,s.zero)-a*b
                if c!=s.zero: v[i]=c
                elif i in v: del v[i]
        return None
def gens(edges,coeff,one):
    vv=sorted(set(v for e in edges for v in e)); res=[]
    for v in vv:
        inc=[(1<<i,1 if v==min(e) else -1) for i,e in enumerate(edges) if v in e]
        g={}
        for k in range(1,len(inc)+1):
            a=coeff.get(k,0)*factorial(k)
            if a:
                for sub in combinations(inc,k):
                    m=0; sg=1
                    for b,s_ in sub: m|=b; sg*=s_
                    g[m]=g.get(m,one-one)+one*(a*sg)
        res.append(g)
    return res
def hilb_field(edges,coeff,one,track):
    gg=gens(edges,coeff,one); s=TSpace(track,one); s.insert({0:one}); fr=list(s.rows.values()); hh=[1]
    for d in range(1,len(edges)+1):
        new=[]
        for w in fr:
            for g in gg:
                z=s.insert(fmul(w,g,one-one))
                if z: new.append(z)
        if not new: break
        hh.append(len(new)); fr=new
    return tuple(hh)
def profile(edges):
    F,lam=field('lam',QQ)
    facs=set()
    def track(a):
        for p in (a.numer,a.denom):
            if not p.is_ground:
                for f,e in p.factor_list()[1]: facs.add(f.monic().as_expr())
    gen=hilb_field(edges,{1:1,2:1,3:lam},F.one,track)
    out={'u':tuple(hilbert(edges,{1:1})),'u+u^3':tuple(hilbert(edges,{1:1,3:1})),'generic (alpha!=0)':gen}
    for fexpr in sorted(facs, key=str):  # sorted, so that the output does not depend on hashing
        L=sp.Symbol('lam'); poly=sp.Poly(fexpr,L)
        if poly.degree()==1:
            r=sp.solve(fexpr,L)[0]; r=Fr(int(r.p),int(r.q))
            out[f'lam={r}']=tuple(hilbert(edges,{1:1,2:1,3:r}))
        else:
            K=QQ.algebraic_field(sp.CRootOf(poly,0))
            th=K.from_sympy(sp.CRootOf(poly,0))
            out[f'lam root of {fexpr}']=hilb_field(edges,{1:1,2:1,3:th},K.one,lambda a:None)
    return out
def macros(src):
    out = {}
    for m in re.finditer(r'\\newcommand\{\\(G\w+)\}\[1\]\{', src):
        i = m.end(); depth = 1
        while depth:
            depth += {'{': 1, '}': -1}.get(src[i], 0); i += 1
        out[m.group(1)] = src[m.end():i - 1]
    return out

def graph(body):
    verts = re.findall(r'\\node\[vertex\]\s*\((\w+)\)', body)
    edges = []; straight = set(); dup = []
    for d in re.findall(r'\\draw(.*?);', body, re.S):
        toks = re.findall(r'\((\w+)\)|(--)|(to\s*(?:\[[^\]]*\])?)', d)
        prev = None; conn = None
        for name, dash, to in toks:
            if name:
                if prev is not None and conn:
                    e = (prev, name)
                    # a straight segment drawn twice is one edge; parallel edges are drawn as distinct to[...] curves
                    if conn == '--' and (e in straight or e[::-1] in straight):
                        dup.append(e)
                    else:
                        edges.append(e)
                        if conn == '--': straight.add(e)
                prev, conn = name, None
            else:
                conn = dash or 'to'
    if dup: print('note: straight segment drawn twice, counted once:', dup)
    idx = {v: i + 1 for i, v in enumerate(verts)}
    assert all(a in idx and b in idx for a, b in edges), (verts, edges)
    return len(verts), [(idx[a], idx[b]) for a, b in edges]


if __name__=='__main__':
    src=open(sys.argv[1]).read(); src=re.sub(r'(?<!\\)%.*','',src)
    G={k:graph(v) for k,v in macros(src).items()}
    tables=re.split(r'\\begin\{(?:sideways)?table\}',src)[1:]
    names=[]
    for tab in tables[:3]:
        names+=re.findall(r'^\s*\\(G\w+)\{[.\d]+\}',tab,re.M)
    for name in names:
        nv,E=G[name]; t=time.monotonic()
        pr=profile(E)
        gen=pr['generic (alpha!=0)']
        special={k:v for k,v in pr.items() if k not in ('u','u+u^3','generic (alpha!=0)') and v!=gen}
        spurious=[k for k,v in pr.items() if k.startswith('lam') and v==gen]
        print(f'{name} |V|={nv} |E|={len(E)} edges={E}', flush=True)
        print(f'   u: {pr["u"]}  u+u^3: {pr["u+u^3"]}  generic: {gen}  special: {special}  spurious: {spurious}  ({time.monotonic()-t:.1f}s)', flush=True)
