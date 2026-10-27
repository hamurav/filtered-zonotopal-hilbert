"""Boundary states and A_j polynomials of marked graphs; neck constants (Theorem 4.16)."""
import sys
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../certificates')
from fractions import Fraction
from direct_hilbert import generators, mul, Space, hilbert
from collections import Counter

def filtration_spaces(gens, nedges):
    sp=Space(); sp.insert({0:Fraction(1)}); frontier=list(sp.rows.values()); levels=[dict(sp.rows)]
    for d in range(1,nedges+2):
        new=[]
        for w in frontier:
            for g in gens:
                z=sp.insert(mul(w,g))
                if z: new.append(z)
        levels.append(dict(sp.rows)); frontier=new
        if not new: break
    return levels

def rank(vecs):
    s=Space()
    for v in vecs: s.insert(v)
    return len(s.rows)

def boundary_states(edges, v, coeff):
    """Attach a leaf edge r at v; return Counter of states (p,q) of the open filtration."""
    leaf=max(x for e in edges for x in e)+1
    E=list(edges)+[(v,leaf)]; ridx=len(E)-1; rbit=1<<ridx
    verts=sorted(set(x for e in E for x in e))
    gens_all=generators(E,coeff)
    gens=[g for x,g in zip(verts,gens_all) if x!=leaf]
    P=filtration_spaces(gens,len(E))
    top=len(P)-1
    def Pv(i): return list(P[min(i,top)].values()) if i>=0 else []
    def Qv(i): return [{k:a for k,a in w.items() if not k&rbit} for w in Pv(i)]
    def zP(i): return [ {k|rbit:a for k,a in w.items() if not k&rbit} for w in Pv(i)]
    # dim(Q_i ∩ R_j) = rank(rP_i) + rank(P_j) - rank(P_j + rP_i)
    D={}
    for i in range(-1,top+1):
        for j in range(-1,top+1):
            if i<0 or j<0: D[i,j]=0; continue
            a=zP(i); b=Pv(j)
            D[i,j]=rank(a)+rank(b)-rank(a+b)
    st=Counter()
    for p in range(top+1):
        for q in range(top+1):
            c=D[p,q]-D[p-1,q]-D[p,q-1]+D[p-1,q-1]
            if c: st[p,q]=c
    return st

def A(st,j):
    out=Counter()
    for (p,q),m in st.items():
        if q-p>=j: out[p]+=m
    return out

def neck_constant(st1,st2):
    """(1-t) sum_{j>=3} t^j A_j A_j  as coefficient list"""
    tot=Counter()
    for j in range(3,40):
        a,b=A(st1,j),A(st2,j)
        if not a or not b: break
        for p,x in a.items():
            for q,y in b.items():
                tot[p+q+j]+=x*y; tot[p+q+j+1]-=x*y
    return [tot[k] for k in range(max(tot)+1)] if tot else [0]

if __name__=='__main__':
    co={1:1,2:1}
    G={'K1':([],0) ,'K2 end':([(0,1)],0),'K3':([(0,1),(1,2),(0,2)],0),'double edge':([(0,1),(0,1)],0),
       'theta graph':([(0,1),(0,1),(0,1)],0),'K4':([(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],0),
       'K_{1,3} centre':([(0,1),(0,2),(0,3)],0),'K_{1,4} centre':([(0,1),(0,2),(0,3),(0,4)],0),'C4':([(0,1),(1,2),(2,3),(3,0)],0)}
    for name,(e,v) in G.items():
        if not e:
            continue
        st=boundary_states(e,v,co)
        print(name, sorted(st.items()), 'A3:',dict(A(st,3)), 'A4:',dict(A(st,4)))
