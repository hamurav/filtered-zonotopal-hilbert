import sys, time, itertools
sys.path.insert(0, sys.argv[1])
from boundary_local import root_states, close
from collections import Counter, defaultdict
from fractions import Fraction as Fr
NMAX=int(sys.argv[2])
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
P={'generic (3,7)':{1:1,2:1,3:3,4:7},'b=2 (2,7)':{1:1,2:1,3:2,4:7},'c=5b-5 (3,10)':{1:1,2:1,3:3,4:10},
   'c=b^2+b-1 (3,11)':{1:1,2:1,3:3,4:11},'(2,5)':{1:1,2:1,3:2,4:5},
   'u+u^3+u^4':{1:1,2:0,3:1,4:1},'u+u^3':{1:1,2:0,3:1},'u+u^4':{1:1,2:0,3:0,4:1},'u':{1:1}}
locs={}
for k,co in P.items():
    cf={d:Fr(a) for d,a in co.items()}; cache={}
    def mk(cf=cf,cache=cache):
        def L(dd):
            if dd not in cache: cache[dd]=root_states(dd,cf)
            return cache[dd]
        return L
    locs[k]=mk()
first={}  # pair -> first witness
names=list(P); t0=time.monotonic(); count=Counter()
for n,e in trees(NMAX):
    count[n]+=1
    H={k:hilbert(e,locs[k]) for k in names}
    for a,b in itertools.combinations(names,2):
        if H[a]!=H[b] and (a,b) not in first: first[a,b]=(n,e,H[a],H[b])
print('trees per n',dict(count),f'{time.monotonic()-t0:.0f}s')
for a,b in itertools.combinations(names,2):
    w=first.get((a,b))
    print(f'{a:18s} vs {b:18s}:', 'NOT separated' if w is None else f'n={w[0]} {w[1]} {w[2]} {w[3]}')
