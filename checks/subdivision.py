#!/usr/bin/env python3
"""Subdivision families of the Outlook (Conjecture 8.2, "Subdivision stabilization").

W_n is obtained from a graph W by subdividing a marked edge into a path with
n edges.  For W = theta (three parallel edges), W_n = Theta(1,1,n); for
W = K4, one edge of K4 is subdivided.  The script computes H_n = HS(W_n)
exactly over Q (search/hilbert_direct.py: fraction-free integer
elimination in the edge algebra) and the neck operator of (4.x),

    L(H_n) = H_n - H_{n-1} - 3t H_{n-2} + 2t H_{n-3},

and reports the lowest degree of L(H_n) - L(H_{n-1}).  The Outlook states
that in the computed ranges (theta: n <= 8, K4: n <= 5) this difference is
supported in degrees >= floor(n/2), and that the stable part of L(H_n)
begins with t^3.

All vertices of W_n have degree <= 3, so f matters only through
u + alpha u^2 + beta u^3; after scaling alpha = 1, f = u + u^2 + lam u^3.
The default values lam = 3 and lam = -1/2 avoid the special values
lam = 2 (beta = 2 alpha^2) and lam = 2/3 (the exponential); the two runs
must agree.

Usage:  python3 subdivision.py [--theta NMAX] [--k4 NMAX] [--lam 3 -1/2]
"""

import argparse
import os
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'search'))
from hilbert_direct import hilbert_direct  # noqa: E402


def theta(n):
    """Theta(1,1,n): vertices 0, 1 joined by two edges and a path of n edges."""
    if n == 1:
        path = [(0, 1)]
    else:
        path = [(0, 2)] + [(i, i + 1) for i in range(2, n)] + [(n, 1)]
    return n + 1, [(0, 1), (0, 1)] + path


def k4(n):
    """K4 on 0..3 with the edge 01 subdivided into a path of n edges."""
    rest = [(0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    if n == 1:
        path = [(0, 1)]
    else:
        path = [(0, 4)] + [(i, i + 1) for i in range(4, n + 2)] + [(n + 2, 1)]
    return n + 3, rest + path


def padd(*ps):
    out = [0] * max(len(p) for p in ps)
    for p in ps:
        for i, c in enumerate(p):
            out[i] += c
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def scale(p, c, shift=0):
    return [0] * shift + [c * x for x in p]


def L(H, n):
    """L(H_n) = H_n - H_{n-1} - 3t H_{n-2} + 2t H_{n-3}."""
    return padd(H[n], scale(H[n - 1], -1), scale(H[n - 2], -3, 1),
                scale(H[n - 3], 2, 1))


def fmt(p):
    terms = []
    for i, c in enumerate(p):
        if c:
            terms.append(f'{c:+d}' + ('' if i == 0 else f't^{i}'))
    return ' '.join(terms) if terms else '0'


def lowest(p):
    for i, c in enumerate(p):
        if c:
            return i
    return None


def run(name, family, nmax, lam):
    f = {1: 1, 2: 1, 3: lam}
    H = {}
    print(f'\n{name}, f = u + u^2 + ({lam}) u^3')
    for n in range(1, nmax + 1):
        nv, E = family(n)
        t0 = time.time()
        H[n] = hilbert_direct(nv, E, f)
        print(f'  n={n}: |V|={nv} |E|={len(E)}  H = {H[n]}  '
              f'({time.time() - t0:.1f}s)', flush=True)
    ok = True
    for n in range(4, nmax + 1):
        print(f'  L(H_{n}) = {fmt(L(H, n))}')
    for n in range(5, nmax + 1):
        D = padd(L(H, n), scale(L(H, n - 1), -1))
        lo = lowest(D)
        good = lo is None or lo >= n // 2
        ok &= good
        print(f'  L(H_{n}) - L(H_{n-1}) = {fmt(D)}   lowest degree {lo}, '
              f'floor(n/2) = {n // 2}: {"ok" if good else "FAILS"}')
    return H, ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--theta', type=int, default=8)
    ap.add_argument('--k4', type=int, default=5)
    ap.add_argument('--lam', nargs='+', default=['3', '-1/2'])
    args = ap.parse_args()
    lams = [Fraction(x) for x in args.lam]
    allok = True
    series = {}
    for lam in lams:
        for name, fam, nmax in (('Theta(1,1,n)', theta, args.theta),
                                ('K4, one edge subdivided', k4, args.k4)):
            H, ok = run(name, fam, nmax, lam)
            allok &= ok
            series[(name, lam)] = H
    for name in ('Theta(1,1,n)', 'K4, one edge subdivided'):
        same = all(series[(name, lam)] == series[(name, lams[0])]
                   for lam in lams)
        print(f'{name}: the values of lam give '
              f'{"the same series" if same else "DIFFERENT series"}')
    print('all checks passed' if allok else 'SOME CHECK FAILED')
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()
