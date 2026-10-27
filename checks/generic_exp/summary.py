"""Recompute, from the recorded results, every number quoted in Section 6.10 from the modular
computations: the exponential on the graphs with at most five vertices, Remark 6.23 (six
vertices), the evidence for Conjecture 6.25 (h_2 and h_3), and the unimodality count of the
Outlook.  Needs networkx.  Usage: python3 summary.py"""
import glob, json
from collections import Counter
from math import comb
import networkx as nx

P1 = json.load(open('results/phase1.json'))
simple = {k: d for k, d in P1.items() if len(set(map(tuple, d['edges']))) == len(d['edges'])
          and not k.startswith('E')}
non = sorted(k for k, d in simple.items() if d['exp'] != d['generic'])
print(f'simple graphs with at most five vertices: {len(simple)}; exponential exceptional for {len(non)}:')
print('   ', ', '.join(non))

S = json.load(open('results/six_6_15.json'))
assert len(S) == 112 and all(r['generic'] == r['generic2'] for r in S)
bad6 = [r for r in S if r['exp'] != r['generic']]
first = Counter(next(i for i in range(30) if (r['generic'] + [0] * 30)[i] != (r['exp'] + [0] * 30)[i])
                for r in bad6)
good6 = [r for r in S if r['exp'] == r['generic']]
trees = sum(1 for r in good6 if len(r['edges']) == 5)
print(f'six vertices: exponential smaller for {len(bad6)} of 112; first difference at h_k: {dict(sorted(first.items()))}')
print(f'   generic for {len(good6)}: {trees} trees, max edges {max(len(r["edges"]) for r in good6)}; '
      f'min edges among the others {min(len(r["edges"]) for r in bad6)}; '
      f'all graphs with >= 9 edges affected: {all(r["exp"] != r["generic"] for r in S if len(r["edges"]) >= 9)}')

B = json.load(open('results/bonus_low.json'))
c2 = Counter(); c3 = Counter(); fail3 = Counter()
for r in B:
    G = nx.Graph([tuple(e) for e in r['edges']]); n = r['n']
    if nx.is_tree(G):
        continue
    deg = dict(G.degree())
    p = sum(2 if deg[next(iter(G[v]))] == 2 else 1 for v in G if deg[v] == 1)
    cycle = all(d == 2 for d in deg.values())
    b2 = comb(n + 1, 2) - p - cycle
    excluded2 = (n == 4 and G.number_of_edges() == 4)       # C_4 and the paw
    for key in ('generic', 'exp'):
        c2[(key, excluded2, r[key][2] == b2)] += 1
    m = G.number_of_edges() - n + 1
    d2 = sum(1 for d in deg.values() if d == 2)
    e22 = sum(1 for a, b in G.edges() if deg[a] == 2 and deg[b] == 2)
    b3 = comb(n + 2, 3) - d2 - e22 - (max(deg.values()) == 3)
    applies = (not nx.has_bridges(G)) and m >= 3 and not (n == 4 and m == 3)
    for key in ('generic', 'exp'):
        if applies:
            c3[(key, r[key][3] == b3)] += 1
    if not applies and r['generic'][3] != b3:
        kind = 'bridge' if nx.has_bridges(G) else ('cycle' if m == 1 else ('K4' if n == 4 else f'cyclomatic {m}'))
        fail3[kind] += 1
    if not applies and nx.has_bridges(G) is False and m == 2 and r['generic'][3] == b3:
        fail3['cyclomatic 2, attained'] += 1
print('h_2 (non-trees, 4..7 vertices): (series, excluded C4/paw, bound attained) ->', dict(c2))
print('h_3 (bridgeless, cyclomatic >= 3, not K4): (series, attained) ->', dict(c3))
print('h_3 outside the hypotheses, bound not attained (generic):', dict(fail3))

series = []
for d in P1.values():
    series += [d['generic'], d['exp']]
for fn in glob.glob('results/scan_p101_*.json'):
    series += [json.loads(h) for h in json.load(open(fn))['exceptional']]
for r in S:
    series += [r['generic'], r['exp']]


def unimodal(h):
    i = 0
    while i + 1 < len(h) and h[i + 1] >= h[i]:
        i += 1
    while i + 1 < len(h) and h[i + 1] <= h[i]:
        i += 1
    return i == len(h) - 1


print(f'unimodality: {sum(map(unimodal, series))} of {len(series)} series unimodal')
