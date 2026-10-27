import sys, itertools
sys.path.insert(0, sys.argv[1]); sys.path.insert(0,'.')
from collections import defaultdict, Counter
exec(open('collision_moves.py').read().split('byn=defaultdict(list)')[0])
byn=defaultdict(dict)
for n,e in W.trees(NMAX,D=D):
    byn[n][W.canon(adjof(e))]=(e,tuple(W.hilbert(e,L) for L in locs))
tot=Counter(); good=Counter()
for n,d in byn.items():
    for k,(e,H) in d.items():
        a=adjof(e)
        for (u,P),(v,Q) in itertools.combinations(pendant_paths(a),2):
            if u==v or len(P)==len(Q): continue
            ee=[x for x in e if not ({x[0],x[1]}=={u,P[-1]} or {x[0],x[1]}=={v,Q[-1]})]+[(v,P[-1]),(u,Q[-1])]
            k2=W.canon(adjof(ee))
            if k2==k: continue
            typ=(min(len(P),len(Q)),max(len(P),len(Q)),dist(a,u,v))
            tot[typ]+=1; good[typ]+= d[k2][1]==H
for typ in sorted(tot):
    if good[typ]: print(typ, f'{good[typ]}/{tot[typ]}')
