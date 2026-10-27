#!/usr/bin/env python3
"""Certified annihilators of the path and cycle filtrations (Sections 7.3-7.5).

Checks, for small cases, the descriptions of annihilators under the pairing
<g, h> = [x_0 ... x_(N-1)] gh of (7.8):

  Lemma 7.4      (B_(2r-1)^(r-1))^perp = V_r, dim V_r = 3^(r-1)          (7.11)
                 (P_(2r-1)^r)^perp = k Omega_(2r-2) + (1-x_1)(1+x_(2r-3)) V_(r-1)
                                                                          (7.12)
                 P_(2r-1)^(r+1) = C_(2r-1), P_(2r)^(r+1) = C_(2r), C_N^perp = k Omega_(N-1)
                                                                    (7.9), (7.13)
  odd cycles     F_(m+1)^perp = k 1 + k Omega_(2m) for C_(2m+1)             (7.17)
                 and F_(m+2) = CC_(2m+1)
  Proposition 7.5  F_m^perp = k Omega_(2m-1) + k(a_0 + b_0) + Z_odd(a) + Z_odd(b)
                 for C_(2m)                                               (7.20)
                 and F_(m+1) = CC_(2m); for C_6 the six annihilators displayed after
                 Proposition 7.5, dim F_3 = 58.

Method (a certificate, not only a consistency check).  Let S be the list of
all products of at most d generators, computed exactly in integers, and let
A be the proposed spanning set of the annihilator.
  (a) Every pairing <a, s>, a in A, s in S, is computed exactly and must be 0.
  (b) The ranks r_S, r_A of S and A are computed modulo the prime p = 1048573.
If r_S + r_A = 2^N, then, since ranks modulo p are at most ranks over Q,
the pairing is nondegenerate, and (a) gives rk_Q S + rk_Q A <= 2^N,
    2^N = r_S + r_A <= rk_Q S + rk_Q A <= 2^N,
so equality holds throughout: span S = F_d has dimension r_S and its
annihilator is exactly span A.  The script also checks r_A against the
dimension stated in the paper.

Usage:  python3 annihilators.py [--rmax 5] [--mmax 5]
"""

import argparse
import itertools
import sys
import time

import numpy as np

P = 1048573  # largest prime below 2^20


# ---------------------------------------------------------------------------
# exact arithmetic in A_N = Q[x_0..x_(N-1)]/(x_i^2), dense integer vectors

class Alg:
    def __init__(self, N):
        self.N, self.D = N, 1 << N
        self.idx = np.arange(self.D)
        self.full = self.D - 1
        self._sel = {}

    def vec(self, terms):
        v = np.zeros(self.D, dtype=object)
        for m, c in terms.items():
            v[m] += c
        return v

    def one(self):
        return self.vec({0: 1})

    def x(self, i):
        return self.vec({1 << (i % self.N): 1})

    def mul(self, a, b):
        r = np.zeros(self.D, dtype=object)
        for m in np.flatnonzero(b):
            m = int(m)
            if m not in self._sel:
                s = self.idx[(self.idx & m) == 0]
                self._sel[m] = (s, s | m)
            s, t = self._sel[m]
            r[t] += a[s] * b[m]
        return r

    def prod(self, *vs):
        out = self.one()
        for v in vs:
            out = self.mul(out, v)
        return out

    def y(self, i, j):
        """x_j - x_i - x_i x_j."""
        return self.x(j) - self.x(i) - self.mul(self.x(i), self.x(j))

    def omega(self):
        """Omega_(N-1) = sum_i prod_(j != i) x_j."""
        return self.vec({self.full ^ (1 << i): 1 for i in range(self.N)})


def monomials(alg, gens, d):
    """All products of at most d of the generators (multisets), exactly."""
    out = [alg.one()]
    layer = {(): alg.one()}
    for _ in range(d):
        new = {}
        for w, v in layer.items():
            start = w[-1] if w else 0
            for i in range(start, len(gens)):
                new[w + (i,)] = alg.mul(v, gens[i])
        out.extend(new.values())
        layer = new
    return out


def pairings_vanish(alg, A, S):
    """Exact check that <a, s> = 0 for all a in A, s in S."""
    Am = np.array(A, dtype=object)
    Sm = np.array(S, dtype=object)[:, alg.full ^ alg.idx]
    a, s = max(abs(int(c)) for c in Am.flat), max(abs(int(c)) for c in Sm.flat)
    if a * s * alg.D < 2 ** 62:  # int64 sums cannot overflow
        Am, Sm = Am.astype(np.int64), Sm.astype(np.int64)
    return not np.any(Am.dot(Sm.T))


# ---------------------------------------------------------------------------
# ranks modulo P (float64 products are exact: inner dimension * P^2 < 2^53)

def _mm(A, B):
    if A.shape[1] == 0:
        return np.zeros((A.shape[0], B.shape[1]), dtype=np.int64)
    assert A.shape[1] * (P - 1) ** 2 < 2 ** 53
    C = A.astype(np.float64) @ B.astype(np.float64)
    return np.mod(C, P).astype(np.int64)


def _rref(M):
    M = M.copy()
    R, piv = [], []
    for i in range(M.shape[0]):
        row = M[i]
        nz = np.flatnonzero(row)
        if nz.size == 0:
            continue
        j = nz[0]
        row = row * pow(int(row[j]), P - 2, P) % P
        rest = M[i + 1:]
        rr = np.flatnonzero(rest[:, j])
        if rr.size:
            rest[rr] = (rest[rr] - np.outer(rest[rr, j], row)) % P
        for k in range(len(R)):
            if R[k][j]:
                R[k] = (R[k] - R[k][j] * row) % P
        R.append(row)
        piv.append(int(j))
    return (np.array(R, dtype=np.int64).reshape(len(R), M.shape[1]), piv)


def rank_mod(vectors, chunk=128):
    rows = (np.array(vectors, dtype=object) % P).astype(np.int64)
    D = rows.shape[1]
    B, piv = np.zeros((0, D), dtype=np.int64), []
    for s in range(0, len(rows), chunk):
        C = rows[s:s + chunk]
        if piv:
            C = (C - _mm(C[:, piv], B)) % P
        N, npiv = _rref(C)
        if npiv:
            if piv:
                B = (B - _mm(B[:, npiv], N)) % P
            B = np.vstack([B, N])
            piv = piv + npiv
    return len(piv)


# ---------------------------------------------------------------------------

RESULTS = []


def certify(label, alg, gens, d, A, expected):
    t0 = time.time()
    S = monomials(alg, gens, d)
    vanish = pairings_vanish(alg, A, S)
    rS, rA = rank_mod(S), rank_mod(A)
    ok = vanish and rS + rA == alg.D and rA == expected
    RESULTS.append(ok)
    print(f'{"ok  " if ok else "FAIL"} {label}: N={alg.N}, {len(S)} products, '
          f'pairings vanish: {vanish}, dim F = {rS}, dim of annihilator = {rA} '
          f'(expected {expected}), {rS} + {rA} = {rS + rA} vs 2^N = {alg.D} '
          f'({time.time() - t0:.1f}s)', flush=True)
    return rS


def V(alg, variables):
    """Spanning set of V_r (7.10) on the path with the given 2r-1 variables."""
    r = (len(variables) + 1) // 2
    facs = []
    for j in range(r - 1):
        a, b, c = variables[2 * j], variables[2 * j + 1], variables[2 * j + 2]
        facs.append([alg.one(), alg.x(b),
                     alg.x(a) - alg.x(c) - alg.mul(alg.x(a), alg.x(c))])
    return [alg.prod(*ch) for ch in itertools.product(*facs)]


def iota(alg, J, var):
    """iota a_J = sum_nu (-1)^(nu-1) a_(J minus j_nu)   (7.18)."""
    out = np.zeros(alg.D, dtype=object)
    for nu, jn in enumerate(J):
        out = out + (-1) ** nu * alg.prod(*[alg.x(var(j)) for j in J if j != jn])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--rmax', type=int, default=5)
    ap.add_argument('--mmax', type=int, default=5)
    args = ap.parse_args()

    print('Lemma 7.4 (paths)')
    for r in range(1, args.rmax + 1):
        N = 2 * r - 1
        alg = Alg(N)
        ys = [alg.y(i, i + 1) for i in range(N - 1)]
        closed = [alg.x(0)] + ys + [alg.x(N - 1)]
        certify(f'(7.11) r={r}: (B^{r-1}_{N})^perp = V_{r}', alg, closed,
                r - 1, V(alg, list(range(N))), 3 ** (r - 1))
        if r >= 2:
            pre = alg.mul(alg.one() - alg.x(1), alg.one() + alg.x(N - 2))
            A = [alg.omega()] + [alg.mul(pre, v)
                                 for v in V(alg, list(range(1, N - 1)))]
            certify(f'(7.12) r={r}: (P^{r}_{N})^perp', alg, ys, r, A,
                    1 + 3 ** (r - 2))
        certify(f'(7.13) r={r}: P^{r+1}_{N} = C_{N}', alg, ys, r + 1,
                [alg.omega()], 1)
        alg = Alg(N + 1)
        ys = [alg.y(i, i + 1) for i in range(N)]
        certify(f'(7.13) r={r}: P^{r+1}_{N+1} = C_{N+1}', alg, ys, r + 1,
                [alg.omega()], 1)

    print('odd cycles, (7.17) and Theorem 7.3')
    for m in range(1, args.mmax):
        n = 2 * m + 1
        alg = Alg(n)
        Y = [alg.y(i, (i + 1) % n) for i in range(n)]
        certify(f'(7.17) n={n}: F_{m+1}^perp = k1 + kOmega', alg, Y, m + 1,
                [alg.one(), alg.omega()], 2)
        certify(f'n={n}: F_{m+2} = CC_{n}', alg, Y, m + 2, [alg.omega()], 1)

    print('even cycles, Proposition 7.5 and Theorem 7.3')
    for m in range(2, args.mmax + 1):
        n = 2 * m
        alg = Alg(n)
        Y = [alg.y(i, (i + 1) % n) for i in range(n)]
        A = [alg.omega(), alg.x(0) + alg.x(1)]
        for var in (lambda j: 2 * j, lambda j: 2 * j + 1):
            for k in range(2, m + 1, 2):
                for J in itertools.combinations(range(m), k):
                    A.append(iota(alg, J, var))
        certify(f'(7.20) n={n}: F_{m}^perp', alg, Y, m, A, 2 ** (m - 1) + 2)
        if m == 3:
            x = alg.x
            six = [alg.omega(), x(0) + x(1), x(2) - x(0), x(4) - x(0),
                   x(3) - x(1), x(5) - x(1)]
            dim = certify('C_6: the six displayed annihilators', alg, Y, 3,
                          six, 6)
            RESULTS.append(dim == 58)
            print(f'{"ok  " if dim == 58 else "FAIL"} C_6: dim F_3 = {dim}')
        certify(f'n={n}: F_{m+1} = CC_{n}', alg, Y, m + 1, [alg.omega()], 1)

    ok = all(RESULTS)
    print('all checks passed' if ok else 'SOME CHECK FAILED')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
