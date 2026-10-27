import sys, time, itertools
sys.path.insert(0, sys.argv[1])
from boundary_local import root_states, close
from collections import Counter, defaultdict
from fractions import Fraction as Fr

# --- nonisomorphic trees with max degree <= 4, by leaf addition + canonical form
def canon(adj):
    n=len(adj)
    # find center(s)
    deg={v:len(adj[v]) for v in adj}; leaves=[v for v in adj if deg[v]<=1]; rem=n
    while rem>2:
        nl=[]
        for v in leaves:
            rem-=1
            for w in adj[v]:
                deg[w]-=1
                if deg[w]==1: nl.append(w)
            deg[v]=0
        leaves=nl
    def enc(v,p): return '('+''.join(sorted(enc(w,v) for w in adj[v] if w!=p))+')'
    return min(enc(c,None) for c in leaves) if len(leaves)==1 else min(enc(leaves[0],None),enc(leaves[1],None))
def trees(nmax,D=4):
    cur={canon({0:[1],1:[0]}):[(0,1)]}
    yield 2,[(0,1)]
    for n in range(3,nmax+1):
        nxt={}
        for edges in cur.values():
            adj=defaultdict(list)
            for u,v in edges: adj[u].append(v); adj[v].append(u)
            for v in list(adj):
                if len(adj[v])>=D: continue
                e2=edges+[(v,n-1)]
                a2=defaultdict(list)
                for u,w in e2: a2[u].append(w); a2[w].append(u)
                k=canon(a2)
                if k not in nxt: nxt[k]=e2
        cur=nxt
        for e in cur.values(): yield n,e
def hilbert(edges,local):
    adj=defaultdict(list)
    for u,v in edges: adj[u].append(v); adj[v].append(u)
    memo={}
    def msg(v,p):
        if (v,p) in memo: return memo[v,p]
        ch=[msg(w,v) for w in adj[v] if w!=p]
        ans=Counter()
        for terms in itertools.product(*(tuple(S.items()) for S in ch)):
            base=sum(pq[0] for pq,m in terms)
            dd=tuple(sorted(pq[1]-pq[0] for pq,m in terms))
            wt=1
            for pq,m in terms: wt*=m
            for (a,b),m in local(dd).items(): ans[a+base,b+base]+=wt*m
        memo[v,p]=ans; return ans
    u,v=edges[0]
    return tuple(close(msg(u,v),msg(v,u)))
