#!/usr/bin/env python3
"""Exact certificates for deformed zonotopal Hilbert strata.

The main calculations are over rational-function fields and exact residue
number fields, not over sampled floating-point or finite-field parameters.
Dependencies: Python 3, SymPy; C++17 compiler and GMP development libraries
for the two eight-vertex unicycles. Run from any directory.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
from pathlib import Path
import json
import subprocess
import sympy as sp
from sympy.polys.fields import field
from sympy.polys.domains import QQ
from boundary_local import root_states, Tracker, add, shift, scale, trans, close
from direct_hilbert import GRAPHS, hilbert
from direct_symbolic_degree import hilbert_symbolic, F as SF, z as sz

HERE = Path(__file__).resolve().parent
OUT = HERE / 'verified_output'
OUT.mkdir(exist_ok=True)
RESULTS = {}


def pretty(C):
    return [[p, q, int(w)] for (p,q), w in sorted(C.items())]


def assert_same_factors(tr, allowed):
    # Normalize expressions, rather than relying on formatting of field elements.
    actual = {sp.srepr(sp.sympify(x)) for x in tr.factors}
    expect = {sp.srepr(sp.sympify(x)) for x in allowed}
    assert actual <= expect, (tr.factors, allowed)


def times1t(h, k):
    for _ in range(k):
        h = [h[0]] + [h[i]+h[i-1] for i in range(1,len(h))] + [h[-1]]
    return h


def expected_binary(special=False):
    def C(*terms):
        return Counter({(p,q):w for p,q,w in terms})
    rr={
      (0,0): C((0,1,4)),
      (0,1): C((0,2,2),(1,1,2)),
      (0,2): C((0,2,2),(1,2,2)),
      (0,3): C((0,3,2),(1,2,2)),
      (1,1): C((0,2,1),(1,2,3)),
      (1,3): C((0,3,1),(1,3,1),(1,4,1),(2,2,1)),
      (2,3): C((0,3,1),(1,4,1),(2,3,2)),
      (3,3): C((0,3,1),(1,4,1),(2,3,1),(3,4,1)),
    }
    if special:
        rr[1,2]=C((0,2,1),(1,3,2),(2,2,1))
        rr[2,2]=C((0,2,1),(1,3,1),(2,3,2))
    else:
        rr[1,2]=C((0,3,1),(1,2,1),(1,3,1),(2,2,1))
        rr[2,2]=C((0,3,1),(1,3,1),(2,2,1),(2,3,1))
    return rr


def verify_subcubic():
    F,b=field('b',QQ)
    tables={}
    for special in [False,True]:
        tr=Tracker();expected=expected_binary(special)
        coeff={1:F.one,2:F.one,3:F(2) if special else b}
        for dd,answer in expected.items():
            found=root_states(dd,coeff,tr)
            assert found==answer,(special,dd,found,answer)
            assert all(0 <= q-p <= 3 for p,q in found)
        assert_same_factors(tr, [] if special else ['b - 2'])
        tables['special' if special else 'generic']={str(k):pretty(v) for k,v in expected.items()}
        print('Subcubic binary table:', 'b=2' if special else 'b!=2', 'verified; excluded factors:',sorted(tr.factors),flush=True)
    RESULTS['subcubic_binary_tables']=tables


def tree_H(R):
    # T2: branching centers of degree three joined by one edge.
    # Side arms (1,1)|(3,3) versus (1,3)|(1,3).
    A=R[1,1]
    B=add(R[2,2],scale(shift(R[1,2],1),6),scale(shift(R[1,1],2),9))
    C=add(R[1,2],scale(shift(R[1,1],1),3))
    H2a,H2b=close(A,B),close(C,C)
    assert H2a==H2b
    # T1: centers of degree four and three joined by a two-edge path.
    # Side arms (1,1,1)|(2,2) versus (2,1,1)|(2,1).
    A1=R[1,1,1]
    B1=add(R[2,2],scale(shift(R[0,2],1),2),shift(R[0,0],2))
    A2=add(R[1,1,2],shift(R[0,1,1],1))
    B2=add(R[1,2],shift(R[0,1],1))
    H1a,H1b=close(trans(A1),B1),close(trans(A2),B2)
    assert H1a==H1b
    return H1a,H2a


def verify_tree_pairs():
    F,b,c=field('b,c',QQ)
    U,z=field('z',QQ)
    cases=[
      ('generic complement',F,b,c,
       ['b - 2','b - c/5 - 1','b**2 + b - c - 1'],0),
      ('b=2, c!=5',U,U(2),z,['z - 5'],1),
      ('c=5b-5, b!=2',U,z,5*z-5,['z - 2'],0),
      ('c=b^2+b-1, b!=2',U,z,z*z+z-1,['z - 2'],0),
      ('b=2,c=5',U,U(2),U(5),[],2),
    ]
    expected_H1=[
      [1,10,48,143,229,81],
      [1,10,48,138,225,90],
      [1,10,47,134,220,100],
    ]
    tree_results=[]
    ds=[(1,1),(1,2),(2,2),(0,2),(0,0),(0,1),(1,1,1),(1,1,2),(0,1,1)]
    for label,K,bb,cc,allowed,kind in cases:
        tracker=Tracker();R={}
        for dd in ds:R[dd]=root_states(dd,{1:K.one,2:K.one,3:bb,4:cc},tracker)
        assert_same_factors(tracker,allowed)
        H1,H2=tree_H(R)
        assert H1==expected_H1[kind]
        assert H2==([1,10,49,146,225,81] if kind==0 else [1,10,48,138,234,81])
        tree_results.append({'cell':label,'T1':H1,'T2':H2,'excluded_factors':sorted(tracker.factors),'states':{str(d):pretty(S) for d,S in R.items()}})
        print(label,'T1=',H1,'T2=',H2,'excluded:',sorted(tracker.factors),flush=True)
    # Zero-quadratic locus: use scaling to set b=1 (with c arbitrary),
    # or b=0,c=1; the origin is separate. Field extensions preserve ranks.
    for label,bb,cc,H1expected,H2expected in [
      ('a=0,b!=0',U.one,z,[1,10,47,115,157,122,51,9],[1,10,44,106,150,126,60,14,1]),
      ('a=b=0,c!=0',U.zero,U.one,[1,10,41,94,135,126,75,26,4],[1,9,36,84,126,126,84,36,9,1]),
      ('a=b=c=0',U.zero,U.zero,[1,9,36,84,126,126,84,36,9,1],[1,9,36,84,126,126,84,36,9,1]),
    ]:
        tr=Tracker();co={1:U.one,2:U.zero,3:bb,4:cc}
        A=root_states((1,1,1),co,tr);B=root_states((1,1),co,tr)
        assert not tr.factors
        H1=times1t(close(trans(A,False),B),2)
        H2=times1t(close(B,B),4)
        assert (H1,H2)==(H1expected,H2expected)
        tree_results.append({'cell':label,'T1':H1,'T2':H2})
        print(label,'T1=',H1,'T2=',H2,flush=True)
    RESULTS['tree_pair_certificates']=tree_results


def verify_small_examples():
    out={}
    z=sp.Symbol('z')
    for name in ['C3double','K4']:
        generic,bad=hilbert_symbolic(GRAPHS[name],{1:SF.one,2:SF.one,3:sz},verbose=False)
        spec=[]
        for p in bad:
            pp=sp.Poly(sp.sympify(p),z)
            assert pp.degree()==1
            b=-pp.nth(0)/pp.nth(1)
            h=hilbert(GRAPHS[name],{1:1,2:1,3:Fraction(b)})
            spec.append({'b':str(b),'H':h})
        out[name]={'generic':generic,'bad_factors':bad,'specializations':spec,
                   'a=0,b=1':hilbert(GRAPHS[name],{1:1,2:0,3:1}),
                   'a=b=0':hilbert(GRAPHS[name],{1:1,2:0,3:0})}
        print('Parameter-section example',name,json.dumps(out[name]),flush=True)
    RESULTS['parameter_section_corrections']=out


def run_cpp(exe,n,modulus=None,zero_quad=False):
    args=[str(exe),str(n)]
    if modulus is not None:args.append(modulus)
    if zero_quad:args.append('zero')
    text=subprocess.check_output(args,text=True)
    h=json.loads(next(x[2:] for x in text.splitlines() if x.startswith('H ')))
    piv=json.loads(next(x[7:] for x in text.splitlines() if x.startswith('PIVOTS ')))
    return h,piv,text


def verify_unicycles():
    exe=OUT/'certify_unicycle'
    subprocess.run(['g++','-O2','-std=c++17',str(HERE/'certify_unicycle.cpp'),'-lgmpxx','-lgmp','-o',str(exe)],check=True)
    z=sp.Symbol('z');bad=set();generic=[]
    for n in [1,2]:
        h,piv,raw=run_cpp(exe,n)
        assert h==[1,8,32,84,111,4]
        generic.append(h)
        (OUT/f'G8_{n}_generic.txt').write_text(raw)
        for cs in piv:
            poly=sum(sp.Rational(a)*z**i for i,a in enumerate(cs))
            for p,e in sp.factor_list(poly)[1]:bad.add(str(sp.Poly(p,z).monic().as_expr()))
    print('Unicycle generic ranks verified; all excluded factors:',sorted(bad),flush=True)
    expected=[1,8,32,84,111,4];specials=[]
    for idx,p in enumerate(sorted(bad)):
        pp=sp.Poly(sp.sympify(p),z)
        assert pp.is_irreducible
        modulus=','.join(str(pp.nth(i)) for i in range(pp.degree()+1))
        hs=[]
        for n in [1,2]:
            h,_,raw=run_cpp(exe,n,modulus)
            hs.append(h)
            (OUT/f'G8_{n}_residue_{idx}.txt').write_text(raw)
        assert hs[0]==hs[1]
        answer=([1,8,31,77,111,12] if pp==sp.Poly(z-2,z) else
                [1,8,32,84,107,8] if pp==sp.Poly(z-sp.Rational(2,3),z) else expected)
        assert hs[0]==answer
        specials.append({'minimal_polynomial':p,'H':hs[0]})
        print('Unicycle residue field',p,'both H=',hs[0],flush=True)
    for b in [0,1]:
        hs=[]
        for n in [1,2]:
            h,_,raw=run_cpp(exe,n,f'{-b},1',True)
            hs.append(h)
            (OUT/f'G8_{n}_zero_quad_{b}.txt').write_text(raw)
        assert hs[0]==hs[1]
        assert hs[0]==([1,7,24,50,66,55,28,8,1] if b==0 else [1,8,32,72,81,40,6])
        specials.append({'a':0,'b':b,'H':hs[0]})
        print('Unicycle zero-quadratic case',b,'both H=',hs[0],flush=True)
    RESULTS['unicycle_certificate']={'generic':generic[0],'bad_factors':sorted(bad),'specializations':specials}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-unicycles',action='store_true')
    args=parser.parse_args()
    verify_subcubic()
    verify_tree_pairs()
    verify_small_examples()
    if not args.skip_unicycles:verify_unicycles()
    (OUT/'results.json').write_text(json.dumps(RESULTS,indent=2))
    print('ALL REQUESTED EXACT CHECKS PASSED.',flush=True)

if __name__=='__main__':main()
