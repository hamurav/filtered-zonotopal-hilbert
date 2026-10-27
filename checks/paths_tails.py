#!/usr/bin/env python3
"""Paths and tails (Section 4): Table 1, Corollary 4.14, Remark 4.13, Theorem 4.9.

All series are computed exactly over Q, directly from the definition
(search/hilbert_direct.py).  The script checks

1. Table 1: the series of the path with n edges, n <= 8, for f = u and for
   three normalized f with alpha != 0 (which must agree, Proposition 3.14);
2. Corollary 4.14: the direct values for n <= NMAX satisfy
   H_n = H_(n-1) + 3t H_(n-2) - 2t H_(n-3) (n >= 4) with H_0 = 1, H_1 = 1+t,
   H_2 = 1+3t, H_3 = 1+4t+3t^2, i.e. the generating function (4.27); and,
   from the recurrence for n <= 60, deg H_n = ceil(n/2) with leading
   coefficient 3^(r-1) for n = 2r-1 and (m+8) 3^(m-2) for n = 2m;
3. Remark 4.13: for G = K3, L(H_3) = -t^3 + t^4, and for G = K_(1,3)
   marked at its center, H~_2 - H_2 = -t^2 + t^3 (f = u + alpha u^2 with
   alpha = 1 and alpha = -2);
4. Theorem 4.9, for G = K3 and G = K_(1,3) (center): L(H_n) = 0 for
   4 <= n <= 6, H~_n = H_n for n = 3, 4, and H^b_n = H_(n+1) - t H_(n-1)
   for n = 2, 3, 4.

Here G_n is G with a path of length n attached at the marked vertex,
G~_n (n >= 2) is G_(n-2) with two new edges at its last vertex, and G^b_n
is G_n with the last edge doubled; H_n, H~_n, H^b_n are their series, and
L(a_n) = a_n - a_(n-1) - 3t a_(n-2) + 2t a_(n-3) is the operator (4.18).

Usage:  python3 paths_tails.py [--nmax 11]
"""

import argparse
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'search'))
from hilbert_direct import hilbert_direct  # noqa: E402

OK = []


def check(cond, msg):
    OK.append(bool(cond))
    print(f'{"ok  " if cond else "FAIL"} {msg}', flush=True)


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def padd(*ps):
    out = [0] * max(len(p) for p in ps)
    for p in ps:
        for i, c in enumerate(p):
            out[i] += c
    return trim(out)


def tmul(p, c, k=1):
    """c t^k p."""
    return [0] * k + [c * x for x in p]


def L(H, n):
    return padd(H[n], tmul(H[n - 1], -1, 0), tmul(H[n - 2], -3),
                tmul(H[n - 3], 2))


def path(n):
    return n + 1, [(i, i + 1) for i in range(n)]


def straight(G, n):
    """G with a path v = v_0, v_1, ..., v_n attached at vertex 0."""
    nv, E = G
    verts = [0] + list(range(nv, nv + n))
    return nv + n, list(E) + [(verts[i], verts[i + 1]) for i in range(n)]


def forked(G, n):
    nv, E = straight(G, n - 2)
    last = 0 if n == 2 else nv - 1
    return nv + 2, E + [(last, nv), (last, nv + 1)]


def doubled(G, n):
    nv, E = straight(G, n)
    return nv, E + [E[-1]]


K3 = (3, [(0, 1), (1, 2), (0, 2)])
K13 = (4, [(0, 1), (0, 2), (0, 3)])

TABLE1 = {1: [1, 1], 2: [1, 3], 3: [1, 4, 3], 4: [1, 5, 10],
          5: [1, 6, 16, 9], 6: [1, 7, 23, 33], 7: [1, 8, 31, 61, 27],
          8: [1, 9, 40, 98, 108]}


def binom_row(n):
    row = [1]
    for _ in range(n):
        row = [a + b for a, b in zip(row + [0], [0] + row)]
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--nmax', type=int, default=11)
    args = ap.parse_args()
    fs = {'u+u^2': {1: 1, 2: 1}, 'u+u^2+3u^3': {1: 1, 2: 1, 3: 3},
          'u-2u^2+5u^3': {1: 1, 2: -2, 3: 5}}

    print('Table 1 (paths)')
    for n in range(1, 9):
        nv, E = path(n)
        vals = {name: hilbert_direct(nv, E, f) for name, f in fs.items()}
        lin = hilbert_direct(nv, E, {1: 1})
        check(all(v == TABLE1[n] for v in vals.values()) and
              lin == binom_row(n),
              f'n={n}: {TABLE1[n]} for the three f, (1+t)^{n} for f = u')

    print('Corollary 4.14 (paths)')
    H = {0: [1]}
    t0 = time.time()
    for n in range(1, args.nmax + 1):
        H[n] = hilbert_direct(*path(n), fs['u+u^2'])
    check(H[1] == [1, 1] and H[2] == [1, 3] and H[3] == [1, 4, 3],
          'initial values H_0, ..., H_3')
    check(all(L(H, n) == [0] for n in range(4, args.nmax + 1)),
          f'direct values satisfy the recurrence for 4 <= n <= {args.nmax} '
          f'({time.time() - t0:.0f}s)')
    R = dict(H)
    for n in range(args.nmax + 1, 61):
        R[n] = padd(R[n - 1], tmul(R[n - 2], 3), tmul(R[n - 3], -2))
    good = True
    for n in range(1, 61):
        deg = len(R[n]) - 1
        lead = R[n][-1]
        want = 3 ** ((n + 1) // 2 - 1) if n % 2 else \
            Fraction(n // 2 + 8) * Fraction(3) ** (n // 2 - 2)
        good &= deg == (n + 1) // 2 and lead == want and sum(R[n]) == 2 ** n
    check(good, 'degree ceil(n/2), leading coefficients 3^(r-1) and '
          '(m+8)3^(m-2), and H_n(1) = 2^n, for n <= 60')

    print('Remark 4.13 and Theorem 4.9 (tails)')
    for alpha in (1, -2):
        f = {1: 1, 2: alpha}
        Hk = {n: hilbert_direct(*straight(K3, n), f) for n in range(0, 7)}
        check(L(Hk, 3) == [0, 0, 0, -1, 1],
              f'alpha={alpha}: K3, L(H_3) = -t^3 + t^4')
        check(all(L(Hk, n) == [0] for n in range(4, 7)),
              f'alpha={alpha}: K3, L(H_n) = 0 for n = 4, 5, 6')
        Hs = {n: hilbert_direct(*straight(K13, n), f) for n in range(0, 7)}
        Ht2 = hilbert_direct(*forked(K13, 2), f)
        check(padd(Ht2, tmul(Hs[2], -1, 0)) == [0, 0, -1, 1],
              f'alpha={alpha}: K13, H~_2 - H_2 = -t^2 + t^3')
        check(all(L(Hs, n) == [0] for n in range(4, 7)),
              f'alpha={alpha}: K13, L(H_n) = 0 for n = 4, 5, 6')
        for name, G, Hg in (('K3', K3, Hk), ('K13', K13, Hs)):
            check(all(hilbert_direct(*forked(G, n), f) == Hg[n]
                      for n in (3, 4)),
                  f'alpha={alpha}: {name}, H~_n = H_n for n = 3, 4')
            check(all(hilbert_direct(*doubled(G, n), f) ==
                      padd(Hg[n + 1], tmul(Hg[n - 1], -1))
                      for n in (2, 3, 4)),
                  f'alpha={alpha}: {name}, H^b_n = H_(n+1) - t H_(n-1) '
                  'for n = 2, 3, 4')

    ok = all(OK)
    print('all checks passed' if ok else 'SOME CHECK FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
