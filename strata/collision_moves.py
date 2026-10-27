import sys, time, itertools
sys.path.insert(0, sys.argv[1]); sys.path.insert(0,'.')
from fractions import Fraction as Fr
from collections import defaultdict, Counter
from boundary_local import root_states
import importlib.util
spec=importlib.util.spec_from_file_location('w','d4_witness_lib.py'); W=importlib.util.module_from_spec(spec); spec.loader.exec_module(W)
D=int(sys.argv[2]); NMAX=int(sys.argv[3])
PTS={3:[{1:1},{1:1,3:1},{1:1,2:1,3:2},{1:1,2:1}],
     4:[{1:1},{1:1,4:1},{1:1,3:1},{1:1,3:1,4:1},{1:1,2:1,3:3,4:7},{1:1,2:1,3:2,4:7},{1:1,2:1,3:3,4:10},{1:1,2:1,3:2,4:5}]}[D]
locs=[]
for co in PTS:
    cf={k:Fr(v) for k,v in co.items()}; cache={}
    def L(dd,cf=cf,cache=cache):
        if dd not in cache: cache[dd]=root_states(dd,cf)
        return cache[dd]
    locs.append(L)
def adjof(e):
    a=defaultdict(list)
    for u,v in e: a[u].append(v); a[v].append(u)
    return a
def pendant_paths(a):
    out=[]
    for leaf in [v for v in a if len(a[v])==1]:
        path=[leaf]; prev=None; cur=leaf
        while True:
            nxt=[w for w in a[cur] if w!=prev][0]
            prev,cur=cur,nxt
            if len(a[cur])!=2: break
            path.append(cur)
        if len(a[cur])>=3: out.append((cur,path))   # attached at cur, length len(path)
    return out
def dist(a,s,t):
    seen={s:0}; q=[s]
    for x in q:
        for y in a[x]:
            if y not in seen: seen[y]=seen[x]+1; q.append(y)
    return seen[t]
def swaps(e):
    a=adjof(e); res={}
    pp=pendant_paths(a)
    for (u,P),(v,Q) in itertools.combinations(pp,2):
        if u==v or len(P)==len(Q): continue
        # move path P to v and Q to u
        ee=[x for x in e if not ({x[0],x[1]}=={u,P[-1]} or {x[0],x[1]}=={v,Q[-1]})]
        ee+= [(v,P[-1]),(u,Q[-1])]
        k=W.canon(adjof(ee))
        res.setdefault(k,[]).append((min(len(P),len(Q)),max(len(P),len(Q)),dist(a,u,v)))
    return res
byn=defaultdict(list)
for n,e in W.trees(NMAX,D=D):
    byn[n].append((e,tuple(W.hilbert(e,L) for L in locs)))
stats=Counter(); unexplained=[]
for n in sorted(byn):
    prof=defaultdict(list)
    for e,H in byn[n]: prof[H].append(e)
    for cls in prof.values():
        if len(cls)<2: continue
        keys=[W.canon(adjof(e)) for e in cls]
        for i,j in itertools.combinations(range(len(cls)),2):
            sw=swaps(cls[i]).get(keys[j])
            if sw: stats[(n,)+tuple(sorted(set(sw))[0])]+=1
            else: unexplained.append((n,cls[i],cls[j]))
agg=Counter()
for k,v in stats.items(): agg[k[1:]]+=v
print('collision pairs explained by one swap of pendant paths (lenP,lenQ,dist):',dict(agg))
print('unexplained pairs:',len(unexplained), 'by n:',dict(Counter(u[0] for u in unexplained)))
for u in unexplained[:4]: print(u)
