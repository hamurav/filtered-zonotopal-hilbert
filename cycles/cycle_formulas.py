#!/usr/bin/env python3
"""Exact checks of the cycle formulas of Section 7 against the cycle data.

Reads the series H_n recorded by cycle_hilbert.py (results/cycle_series.json)
and checks, in exact integer arithmetic:

1. the transfer matrix M(t) of (4.13): its characteristic polynomial
   z^3 - z^2 - 3tz + 2t, deg tr_n <= floor(n/2), tr_n(1) = 2^n + (-1)^n l_n
   (proof of Proposition 7.9), and [t^j] tr_n = 1, n, C(n+1,2), C(n+2,3) - 3n
   for j <= 3 and n >= 7 (Theorem 7.1), for all n <= NMAX;
2. Theorem 7.1 on the data: h_1 = n, h_2 = C(n+1,2) - 1, h_3 = C(n+2,3) - 3n
   for n >= 7;
3. Theorem 7.3 on the data: the top term is t^(m+2) for n = 2m+1 and
   (2^(m-1)+1) t^(m+1) for n = 2m;
4. Conjecture 7.8, H_n = tr_n - t^2 mod t^ceil(n/2); the closed formulas
   (7.29), (7.30); the recurrence (7.28) for n >= 6; the initial values
   H_3, H_4, H_5 of Proposition 7.9;
5. the total dimension H_n(1) = 2^n - 1;
6. Remark 7.10: L(b_n) = 0 for n even and L(b_n) = l_(n-4) t^((n-1)/2) (1-t)
   for n odd, where b_n = (-1)^(n+1) l_n t^ceil(n/2), for n <= NMAX;
7. Proposition 7.7: tr M_w(t)^n = sum over the independent sets S of C_n of
   w^c(S) t^|S|, as polynomials in t and w, for 3 <= n <= 18, where M_w is
   M with the entry 2 replaced by w - 1 and c(S) is the number of vertices
   outside S with both neighbors in S.

Usage:  python3 cycle_formulas.py [results/cycle_series.json] [--nmax 40]
"""

import argparse
import json
import os
import sys
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))


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


def pmul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                out[i + j] += a * b
    return trim(out)


def mono(c, k):
    return [0] * k + [c]


def L(seq, n):
    """L(a_n) = a_n - a_(n-1) - 3t a_(n-2) + 2t a_(n-3), the operator (4.18)."""
    return padd(seq[n], pmul([-1], seq[n - 1]), pmul([0, -3], seq[n - 2]),
                pmul([0, 2], seq[n - 3]))


def traces(nmax):
    """tr_n = tr M(t)^n for M = [[0,t,0],[2,0,t],[0,1,1]] (entries in Z[t])."""
    M = [[[0], [0, 1], [0]], [[2], [0], [0, 1]], [[0], [1], [1]]]
    P = [[[1] if i == j else [0] for j in range(3)] for i in range(3)]
    tr = {0: [3]}
    for n in range(1, nmax + 1):
        P = [[padd(*[pmul(P[i][k], M[k][j]) for k in range(3)])
              for j in range(3)] for i in range(3)]
        tr[n] = padd(P[0][0], P[1][1], P[2][2])
    return tr


def weighted_trace(n):
    """tr M_w(t)^n as {(i, j): c} for c t^i w^j, M_w = [[0,t,0],[w-1,0,t],[0,1,1]]."""
    def add(*ps):
        out = {}
        for p in ps:
            for k, c in p.items():
                out[k] = out.get(k, 0) + c
        return {k: c for k, c in out.items() if c}

    def mul(p, q):
        out = {}
        for (a, b), c in p.items():
            for (e, f), d in q.items():
                out[a + e, b + f] = out.get((a + e, b + f), 0) + c * d
        return {k: c for k, c in out.items() if c}

    t, one, zero = {(1, 0): 1}, {(0, 0): 1}, {}
    M = [[zero, t, zero], [{(0, 1): 1, (0, 0): -1}, zero, t], [zero, one, one]]
    P = [[one if i == j else zero for j in range(3)] for i in range(3)]
    for _ in range(n):
        P = [[add(*[mul(P[i][k], M[k][j]) for k in range(3)])
              for j in range(3)] for i in range(3)]
    return add(P[0][0], P[1][1], P[2][2])


def weighted_independence(n):
    """sum over independent sets S of C_n of w^c(S) t^|S|, as {(|S|, c): count}."""
    out = {}
    for S in range(1 << n):
        rot = ((S >> 1) | ((S & 1) << (n - 1)))      # bit i: i+1 in S
        if S & rot:
            continue
        prev = (((S << 1) & ((1 << n) - 1)) | (S >> (n - 1)))   # bit i: i-1 in S
        c = bin(~S & rot & prev & ((1 << n) - 1)).count('1')
        k = (bin(S).count('1'), c)
        out[k] = out.get(k, 0) + 1
    return out


def lucas(nmax):
    luc = {1: 1, 2: 3}
    for k in range(3, nmax + 1):
        luc[k] = luc[k - 1] + luc[k - 2]
    return luc


def coeff(p, j):
    return p[j] if j < len(p) else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('data', nargs='?',
                    default=os.path.join(HERE, 'results', 'cycle_series.json'))
    ap.add_argument('--nmax', type=int, default=40)
    args = ap.parse_args()
    N = args.nmax
    ok = True

    def check(cond, msg):
        nonlocal ok
        ok &= bool(cond)
        print(f'{"ok  " if cond else "FAIL"} {msg}')

    tr = traces(N)
    luc = lucas(N + 4)

    # 1. the transfer matrix
    # Cayley-Hamilton for z^3 - z^2 - 3tz + 2t is the recurrence L(tr_n) = 0.
    check(all(L(tr, n) == [0] for n in range(3, N + 1)),
          f'L(tr_n) = 0 for 3 <= n <= {N} (characteristic polynomial '
          'z^3 - z^2 - 3tz + 2t)')
    check(tr[1] == [1] and tr[2] == [1, 6],
          'tr_0 = 3, tr_1 = 1, tr_2 = 1 + 6t')
    check(all(len(tr[n]) - 1 <= n // 2 for n in range(N + 1)),
          f'deg tr_n <= floor(n/2) for n <= {N}')
    check(all(sum(tr[n]) == 2 ** n + (-1) ** n * luc[n]
              for n in range(1, N + 1)),
          f'tr_n(1) = 2^n + (-1)^n l_n for 1 <= n <= {N}')
    check(all(tr[n][:4] == [1, n, comb(n + 1, 2), comb(n + 2, 3) - 3 * n]
              for n in range(7, N + 1)),
          f'[t^j] tr_n = 1, n, C(n+1,2), C(n+2,3) - 3n (j <= 3) '
          f'for 7 <= n <= {N}')

    # 6. the Lucas part
    b = {n: mono((-1) ** (n + 1) * luc[n], (n + 1) // 2)
         for n in range(1, N + 1)}
    check(all(L(b, n) == [0] for n in range(6, N + 1, 2)) and
          all(L(b, n) == pmul(mono(luc[n - 4], (n - 1) // 2), [1, -1])
              for n in range(7, N + 1, 2)),
          f'Remark 7.10: L(b_n) for 6 <= n <= {N}')

    # 7. the weighted independence form of the trace
    check(all(weighted_trace(n) == weighted_independence(n)
              for n in range(3, 19)),
          'Proposition 7.7: tr M_w(t)^n = sum_S w^c(S) t^|S| '
          'for 3 <= n <= 18')

    # data
    with open(args.data) as fh:
        raw = json.load(fh)
    H = {}
    for k, v in raw.items():
        h = v['h'] if isinstance(v, dict) else v
        if isinstance(v, dict) and not v.get('agree', True):
            check(False, f'n={k}: the primes disagree in {args.data}')
        H[int(k)] = trim(h)
    ns = sorted(H)
    print(f'cycle data for n = {ns[0]}..{ns[-1]} from {args.data}')
    check(ns == list(range(ns[0], ns[-1] + 1)), 'no gaps in the data')

    check(H.get(3) == [1, 3, 2, 1] and H.get(4) == [1, 4, 7, 3] and
          H.get(5) == [1, 5, 14, 10, 1],
          'initial values H_3, H_4, H_5 of Proposition 7.9')
    for n in ns:
        h = H[n]
        res = []
        good = sum(h) == 2 ** n - 1
        res.append('total 2^n-1' if good else 'TOTAL WRONG')
        if n >= 7:
            g = h[:4] == [1, n, comb(n + 1, 2) - 1, comb(n + 2, 3) - 3 * n]
            good &= g
            res.append('Thm 7.1' if g else 'THM 7.1 FAILS')
        m = n // 2
        top = (n + 1) // 2 + 1
        if n % 2:
            g = len(h) - 1 == m + 2 and h[-1] == 1
        else:
            g = len(h) - 1 == m + 1 and h[-1] == 2 ** (m - 1) + 1
        good &= g
        res.append('Thm 7.3' if g else 'THM 7.3 FAILS')
        k = (n + 1) // 2
        cong = all(coeff(h, j) == coeff(tr[n], j) - (j == 2)
                   for j in range(k))
        good &= cong
        res.append(f'Conj 7.8 (mod t^{k})' if cong else 'CONJ 7.8 FAILS')
        if n % 2 == 0:
            A = 2 ** (m - 1) + 1
            closed = padd(tr[n], mono(-1, 2), mono(-(luc[n] + A), m),
                          mono(A, m + 1))
        else:
            closed = padd(tr[n], mono(-1, 2), mono(luc[n] - 1, m + 1),
                          mono(1, m + 2))
        c = closed == h
        good &= c
        res.append('(7.29)/(7.30)' if c else 'CLOSED FORMULA FAILS')
        if n >= 6 and all(j in H for j in range(n - 3, n)):
            if n % 2 == 0:
                rhs = padd([0, 0, 0, 1],
                           pmul(mono(2 ** (n // 2 - 2) + 1, n // 2), [1, -1]))
            else:
                rhs = padd([0, 0, 0, 1],
                           pmul(pmul([luc[n - 4] - 1, 2], [1, -1]),
                                mono(1, (n - 1) // 2)))
            r = L(H, n) == rhs
            good &= r
            res.append('(7.28)' if r else 'RECURRENCE FAILS')
        assert top == len(h) - 1 or not good
        check(good, f'n={n}: ' + ', '.join(res))

    print('all checks passed' if ok else 'SOME CHECK FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
