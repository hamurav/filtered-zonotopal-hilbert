"""Low filtration degrees h_0..h_dmax by brute force: all monomials of degree <= dmax in the
Y_v = f(X_v), computed in Phi_G (sparse dicts), ranks mod the prime 1048573.
Cheap for any graph with up to ~8 vertices; used for the h_2, h_3 statistics."""
from math import factorial
from itertools import combinations, combinations_with_replacement
import numpy as np
from gz import rref_rows, PBIG


def gens_mod(n, edges, f, p=PBIG):
    cf = {}
    for k, a in f.items():
        if hasattr(a, 'denominator'):
            a = a.numerator * pow(a.denominator, p - 2, p)
        cf[k] = int(a) % p
    inc = [[] for _ in range(n)]
    for i, (a, b) in enumerate(edges):
        if a != b:
            inc[a].append((1 << i, -1))
            inc[b].append((1 << i, 1))
    gens = []
    for v in range(n):
        y = {}
        for k, a in cf.items():
            if a == 0 or k > len(inc[v]):
                continue
            c = a * factorial(k) % p
            for sub in combinations(inc[v], k):
                m, s = 0, 1
                for bit, sg in sub:
                    m |= bit
                    s *= sg
                y[m] = (y.get(m, 0) + c * s) % p
        gens.append({m: c for m, c in y.items() if c})
    return gens


def mul(x, y, p=PBIG):
    out = {}
    for i, a in x.items():
        for j, b in y.items():
            if not i & j:
                k = i | j
                out[k] = (out.get(k, 0) + a * b) % p
    return {k: c for k, c in out.items() if c}


def low_hilbert(n, edges, f, dmax=3, p=PBIG):
    Y = gens_mod(n, edges, f, p)
    mons = {(): {0: 1}}
    layers = [[()]]
    for d in range(1, dmax + 1):
        layer = []
        for w in combinations_with_replacement(range(n), d):
            mons[w] = mul(mons[w[:-1]], Y[w[-1]], p)
            layer.append(w)
        layers.append(layer)
    cols = {}
    for v in mons.values():
        for m in v:
            cols.setdefault(m, len(cols))
    ranks = []
    rows = []
    for d in range(dmax + 1):
        for w in layers[d]:
            r = np.zeros(len(cols), dtype=np.int64)
            for m, c in mons[w].items():
                r[cols[m]] = c
            rows.append(r)
        R, piv = rref_rows(np.array(rows), p)
        rows = list(R)
        ranks.append(len(piv))
    return [ranks[0]] + [ranks[d] - ranks[d - 1] for d in range(1, dmax + 1)]


if __name__ == '__main__':
    import itertools
    from fractions import Fraction as Fr
    K5 = list(itertools.combinations(range(5), 2))
    print(low_hilbert(5, K5, {1: 1, 2: 3, 3: 5, 4: 7}, 4), 'expected 1,5,15,35,67')
    print(low_hilbert(5, K5, {k: Fr(1, factorial(k)) for k in range(1, 5)}, 4), 'expected 1,5,15,35,63')
