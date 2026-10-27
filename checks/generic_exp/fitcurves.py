"""Fit the exceptional point sets of the GF(p) scans (alpha != 0 chart) by curves.

For each exceptional series, take its points (lam, mu) and find the weighted-homogeneous
polynomials F(alpha, beta, gamma) of weight w (weights 1, 2, 3), dehomogenized at alpha = 1,
that vanish on all of them: F = sum c_{jk} lam^j mu^k over 2j + 3k <= w.  Reports the smallest
weight with a nonzero solution and the solution space (mod p).  Usage: python3 fitcurves.py JSON..."""
import json, sys
import numpy as np
from gz import rref_rows


def nullspace_mod(M, p):
    R, piv = rref_rows(np.array(M, dtype=np.int64), p)
    ncols = len(M[0])
    free = [j for j in range(ncols) if j not in piv]
    basis = []
    for fj in free:
        v = [0] * ncols
        v[fj] = 1
        for r, pj in zip(R, piv):
            v[pj] = (-int(r[fj])) % p
        basis.append(v)
    return basis


def fit(points, p, wmax=12):
    for w in range(1, wmax + 1):
        mons = [(j, k) for k in range(w // 3 + 1) for j in range((w - 3 * k) // 2 + 1)]
        M = [[pow(l, j, p) * pow(m, k, p) % p for (j, k) in mons] for (l, m) in points]
        ns = nullspace_mod(M, p)
        if ns and len(points) > len(mons) + 3:
            return w, mons, ns
    return None


def show(mons, v, p):
    terms = []
    for (j, k), c in zip(mons, v):
        if c:
            c = c if c <= p // 2 else c - p
            terms.append(f'{c}*lam^{j}*mu^{k}')
    return ' + '.join(terms)


if __name__ == '__main__':
    for fn in sys.argv[1:]:
        d = json.load(open(fn))
        p = d['p']
        print('==', d['name'], 'generic', d['generic'])
        for h, fl in sorted(d['exceptional'].items(), key=lambda x: -len(x[1])):
            pts = [(dict(f)[3], dict(f)[4]) for f in [[tuple(t) for t in g] for g in fl] if 2 in dict(f)]
            other = [g for g in fl if 2 not in dict([tuple(t) for t in g])]
            if len(pts) < 8:
                print(f'  {h}: points {pts} {"(+ alpha=0 points)" if other else ""}')
                continue
            res = fit(pts, p)
            if res is None:
                print(f'  {h}: {len(pts)} points, no curve of weight <= 12')
                continue
            w, mons, ns = res
            print(f'  {h}: {len(pts)} points on weight-{w} curve(s), solution space dim {len(ns)}:')
            for v in ns:
                print('      ', show(mons, v, p))
