import sys
sys.path.insert(0, sys.argv[1]); sys.path.insert(0,'.')
from collections import defaultdict
exec(open('collision_moves.py').read().split('byn=defaultdict(list)')[0])
names={3:['u','u+u^3','u+u^2+2u^3','u+u^2 (generic)'],4:['u','u+u^4','u+u^3','u+u^3+u^4','generic','[u3]g=0','[u4]g=0','[u3]g=[u4]g=0']}[D]
byn=defaultdict(list)
for n,e in W.trees(NMAX,D=D): byn[n].append(tuple(W.hilbert(e,L) for L in locs))
for n in sorted(byn):
    if n<NMAX-2: continue
    Hs=byn[n]
    row=[len({H[i] for H in Hs}) for i in range(len(PTS))]
    full=len(set(Hs))
    print(f'D={D} n={n}: trees {len(Hs)}, full-profile classes {full}; single-point classes: '+', '.join(f'{names[i]}: {row[i]}' for i in range(len(PTS))))
