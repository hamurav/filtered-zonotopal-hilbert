"""Find the two graphs of [NSh, Example 2] among loopless multigraphs with 4 vertices and 6 edges.

NSh give: graded series 1+3t+6t^2+9t^3+8t^4+4t^5+t^6 for both graphs,
exponential series 1+4t+10t^2+14t^3+3t^4 and 1+4t+10t^2+15t^3+2t^4,
generic series 1+4t+10t^2+15t^3+2t^4 and 1+4t+10t^2+16t^3+t^4.
The generic series is computed at f = u+u^2+3u^3+7u^4+11u^5+13u^6 and confirmed at a second point.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement
from hilbert_direct import hilbert_direct

pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
EXP = {k: F(1, 1) for k in range(1, 7)}
fact = 1
for k in range(1, 7):
    fact *= k
    EXP[k] = F(1, fact)
GEN1 = {1: 1, 2: 1, 3: 3, 4: 7, 5: 11, 6: 13}
GEN2 = {1: 1, 2: 2, 3: 5, 4: -3, 5: 17, 6: 4}
target_graded = [1, 3, 6, 9, 8, 4, 1]
seen = {}
for multiset in combinations_with_replacement(pairs, 6):
    if len(set(v for e in multiset for v in e)) < 4:
        continue
    h = hilbert_direct(4, list(multiset), {1: 1})
    if h != target_graded:
        continue
    he = hilbert_direct(4, list(multiset), EXP)
    hg = hilbert_direct(4, list(multiset), GEN1)
    assert hg == hilbert_direct(4, list(multiset), GEN2)
    key = (tuple(he), tuple(hg))
    seen.setdefault(key, multiset)
for (he, hg), m in sorted(seen.items()):
    print('edges', m, ' exp', list(he), ' generic', list(hg))
