#!/usr/bin/env python3
"""Hilbert series H_n of the filtered zonotopal algebra of the cycle C_n.

Data for Conjecture 7.8, and a check of Theorems 7.1 and 7.3.  Every f
with a nonzero quadratic coefficient gives the same series (Section 7),
so one computation per n suffices.

Usage:  python3 cycle_hilbert.py NMIN NMAX [--primes K] [--out FILE]

Method
------
Number the vertices and edges of C_n by Z/n, edge i joining the vertices
i and i+1, and put X_i = phi_i - phi_{i-1}.  After scaling, Y_i = f(X_i)
becomes exp(X_i) - 1 = X_i + X_i^2/2 (Proposition 7.12).  Write
e_i = exp(X_i), so that Y_i = e_i - 1 and (e_i - 1)^3 = 0.  The span of
the monomials of degree <= d in the Y_i equals the span of the monomials
of degree <= d in the e_i, and a monomial with an exponent >= 3 lies in
the span of monomials of no larger degree with exponents <= 2.  Hence F_d
is spanned by the product vectors

    E_m = exp(sum_i m_i X_i) = prod_i (1 + (m_i - m_{i+1}) phi_i),
    m in {0,1,2}^n, |m| <= d.

Frontier.  F_d = F_{d-1} + sum_i e_i F_{d-1}.  If B is a set of m with
|m| = d-1 such that the E_m (m in B) span F_{d-1} modulo F_{d-2}, then
F_d = F_{d-1} + span{E_{m+delta_i} : m in B, i in Z/n}, because
e_i E_m = E_{m+delta_i} and e_i F_{d-2} lies in F_{d-1}.  Only these
candidates are tested.  Candidates with an entry 3 lie in F_{d-1}, and a
candidate with all entries >= 1 equals E_{m-(1,...,1)}; both are skipped.

Rotation blocks.  Over F_p with p = 1 mod n, the rotation R splits the
edge algebra into its eigenspaces, and F_d, being R-invariant, splits
accordingly: dim F_d = sum_j dim P_j F_d with P_j = sum_k w^(-jk) R^k.
The eigenspace for w^j has a basis indexed by the R-orbits of edge sets
whose stabilizer order divides j, so each block has dimension about
2^n/n.  The signed reflection phi_i -> -phi_{-1-i} preserves F_d and
exchanges the blocks j and -j, so only j = 0, ..., floor(n/2) are
computed.  Orbit representatives m of the candidates suffice, because
P_j E_{R^k m} = w^(jk) P_j E_m.  The set B for the next degree is the
union over the blocks of the candidates that gave pivots, together with
their reflections.

Each block keeps a reduced echelon basis modulo p in the form [I | M] on
(pivot columns | free columns); a candidate is reduced to its free
coordinates x_free - x_piv M.  Products are formed in floating point with
integer entries below 2^53, hence exactly.

Ranks modulo p are lower bounds for ranks over Q.  The program repeats the
computation for K primes (default 2) and reports any discrepancy; the
result is evidence, not a certificate.  The total 2^n - 1 (the number of
spanning forests of C_n, the dimension of the algebra) is checked.
"""

import argparse
import json
import platform
import sys
import time

import numpy as np


# ----------------------------------------------------------------------
# arithmetic modulo p on float64 arrays of integers

def modp(x, p, pinv):
    """Reduce an array of integers |x| < 2^53 modulo p, returning [0, p)."""
    q = x * pinv
    np.floor(q, out=q)
    q *= p
    y = x - q
    y[y < 0] += p
    y[y >= p] -= p
    return y


def is_prime(q):
    if q < 2:
        return False
    if q % 2 == 0:
        return q == 2
    i = 3
    while i * i <= q:
        if q % i == 0:
            return False
        i += 2
    return True


def prime_factors(n):
    out, q = [], 2
    while q * q <= n:
        if n % q == 0:
            out.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        out.append(n)
    return out


def primes_one_mod(n, below, count):
    """The `count` largest primes p < below with p = 1 mod n."""
    out = []
    p = below - 1 - ((below - 2) % n)
    while len(out) < count and p > n:
        if is_prime(p):
            out.append(p)
        p -= n
    return out


def root_of_unity(n, p):
    fac = prime_factors(n)
    for g in range(2, p):
        w = pow(g, (p - 1) // n, p)
        if all(pow(w, n // q, p) != 1 for q in fac):
            return w
    raise ValueError('no primitive root of unity')


# ----------------------------------------------------------------------
# rotation orbits of edge sets (bitmasks)

def rotation_data(n):
    full = (1 << n) - 1
    masks = np.arange(1 << n, dtype=np.int64)
    rots = np.empty((n, 1 << n), dtype=np.int64)
    for k in range(n):                       # rots[k][A] = A + k
        rots[k] = ((masks << k) | (masks >> (n - k))) & full
    reps = np.flatnonzero(rots.min(axis=0) == masks)
    size = np.full(len(reps), n, dtype=np.int64)
    for k in range(n - 1, 0, -1):            # the smallest period wins
        size[rots[k][reps] == reps] = k
    stab = n // size
    rotidx = np.empty((len(reps), n), dtype=np.int64)
    for k in range(n):                       # rotidx[o][k] = A_o - k
        rotidx[:, k] = rots[(n - k) % n][reps]
    return reps, stab, rotidx


# ----------------------------------------------------------------------
# candidates m in {0,1,2}^n, kept as canonical rotation representatives

class Codes:
    def __init__(self, n):
        self.n = n
        self.pow3 = 3 ** np.arange(n, dtype=np.int64)
        self.refl = [(-i) % n for i in range(n)]

    def canonical(self, m):
        best = m @ self.pow3
        bestm = m.copy()
        for k in range(1, self.n):
            r = np.roll(m, k, axis=1)
            c = r @ self.pow3
            better = c < best
            best = np.where(better, c, best)
            bestm[better] = r[better]
        return bestm, best

    def reflect(self, m):
        return m[:, self.refl]


def expand(frontier, codes):
    """Canonical representatives of the candidates of the next degree."""
    n = codes.n
    f = frontier.shape[0]
    cand = np.repeat(frontier, n, axis=0)
    cand[np.arange(f * n), np.tile(np.arange(n), f)] += 1
    keep = (cand.max(axis=1) <= 2) & (cand.min(axis=1) == 0)
    cand = cand[keep]
    if cand.shape[0] == 0:
        return cand
    cand, code = codes.canonical(cand)
    _, first = np.unique(code, return_index=True)
    return cand[first]                       # sorted by code


# ----------------------------------------------------------------------
# block vectors and elimination

class Projector:
    """Computes the block coordinates P_j E_m for batches of candidates."""

    def __init__(self, n, p, w, chunk):
        self.n, self.p, self.chunk = n, p, chunk
        self.reps, self.stab, rotidx = rotation_data(n)
        self.rotT = np.ascontiguousarray(rotidx.T)       # (n, N)
        self.J = list(range(n // 2 + 1))
        self.cols = [np.flatnonzero(j % self.stab == 0) for j in self.J]
        self.W = np.array([[pow(w, (-j * k) % n, p) for k in range(n)]
                           for j in self.J], dtype=np.float64)
        N, J = len(self.reps), len(self.J)
        self.vbuf = np.empty((chunk, 1 << n))
        self.gbuf = np.empty((chunk, n, N))
        self.pbuf = np.empty((chunk, J, N))
        self.qbuf = np.empty((chunk, J, N))

    def blocks(self, m):
        """Array P with P[:, j, o] = (P_j E_m)[orbit o], reduced mod p."""
        n, p = self.n, self.p
        b = m.shape[0]
        c = (m - np.roll(m, -1, axis=1)).astype(np.float64)
        v = self.vbuf[:b]
        v[:, 0] = 1.0
        for i in range(n):                   # |entries| <= 2^n
            h = 1 << i
            np.multiply(v[:, :h], c[:, i:i + 1], out=v[:, h:2 * h])
        g = self.gbuf[:b]
        for i in range(b):                   # g[i, k, o] = v[i][A_o - k]
            np.take(v[i], self.rotT, out=g[i])
        P, q = self.pbuf[:b], self.qbuf[:b]
        np.matmul(self.W, g, out=P)          # exact: |sum| < n 2^n p
        np.multiply(P, 1.0 / p, out=q)
        np.floor(q, out=q)
        np.multiply(q, p, out=q)
        np.subtract(P, q, out=P)
        np.add(P, p, out=P, where=P < 0)
        np.subtract(P, p, out=P, where=P >= p)
        return P


def eliminate(S, p, pinv, room):
    """Gaussian elimination of the rows of S in order (S is changed).

    Returns the indices of the rows giving new pivots, the new rows in
    reduced echelon form, and their pivot columns (at most `room`).
    """
    rows, pivs = [], []
    live = np.flatnonzero(S.any(axis=1))
    for idx in live:
        row = S[idx]
        nz = np.flatnonzero(row)
        if nz.size == 0:
            continue
        col = nz[0]
        row *= pow(int(row[col]), p - 2, p)
        row[:] = modp(row, p, pinv)
        later = S[idx + 1:]
        fac = later[:, col].copy()
        if fac.any():
            later -= np.outer(fac, row)
            later[:] = modp(later, p, pinv)
        rows.append(idx)
        pivs.append(col)
        if len(rows) == room:
            break
    if not rows:
        return rows, None, None
    new = S[rows].copy()
    pivs = np.array(pivs, dtype=np.int64)
    for t in range(len(rows) - 1, 0, -1):          # back substitution
        fac = new[:t, pivs[t]].copy()
        if fac.any():
            new[:t] -= np.outer(fac, new[t])
            new[:t] = modp(new[:t], p, pinv)
    return rows, new, pivs


class Block:
    """Reduced echelon basis [I | M] modulo p of a subspace of F_p^N."""

    def __init__(self, cols, p, sub=32):
        """`cols`: the coordinates (orbits) belonging to this block."""
        self.N, self.p, self.pinv, self.sub = len(cols), p, 1.0 / p, sub
        self.piv = np.zeros(0, dtype=np.int64)
        self.free = np.array(cols, dtype=np.int64)
        self.M = np.zeros((0, len(cols)), dtype=np.float64)

    @property
    def r(self):
        return len(self.piv)

    def add(self, B):
        """Add the rows of B; return the indices of rows giving pivots."""
        p, pinv = self.p, self.pinv
        if self.r == self.N:
            return []
        Z = B[:, self.free]
        if self.r:
            Z = modp(Z - B[:, self.piv] @ self.M, p, pinv)
        L, lc = None, np.zeros(0, dtype=np.int64)   # local rows, pivots
        chosen = []
        for s in range(0, Z.shape[0], self.sub):
            S = Z[s:s + self.sub].copy()
            if L is not None:
                S = modp(S - S[:, lc] @ L, p, pinv)
            room = self.N - self.r - len(lc)
            rows, new, pivs = eliminate(S, p, pinv, room)
            if rows:
                if L is not None:
                    L = modp(L - L[:, pivs] @ new, p, pinv)
                    L = np.vstack([L, new])
                else:
                    L = new
                lc = np.concatenate([lc, pivs])
                chosen.extend(s + i for i in rows)
            if self.r + len(lc) == self.N:
                break
        if L is not None:                            # commit
            keep = np.ones(len(self.free), dtype=bool)
            keep[lc] = False
            if self.r:
                M = modp(self.M - self.M[:, lc] @ L, p, pinv)
                self.M = np.vstack([M[:, keep], L[:, keep]])
            else:
                self.M = np.ascontiguousarray(L[:, keep])
            self.piv = np.concatenate([self.piv, self.free[lc]])
            self.free = self.free[keep]
        return chosen


# ----------------------------------------------------------------------

def cycle_hilbert(n, p, chunk=256, log=None):
    w = root_of_unity(n, p)
    proj = Projector(n, p, w, chunk)
    codes = Codes(n)
    J = proj.J
    mult = [1 if (j == 0 or 2 * j == n) else 2 for j in J]
    blocks = [Block(proj.cols[t], p) for t in range(len(J))]
    total = (1 << n) - 1

    def process(cands):
        ranks0 = [bl.r for bl in blocks]
        chosen = set()
        for start in range(0, cands.shape[0], chunk):
            if all(bl.r == bl.N for bl in blocks):
                break
            P = proj.blocks(cands[start:start + chunk])
            for t, bl in enumerate(blocks):
                for idx in bl.add(P[:, t, :]):
                    chosen.add(start + idx)
        gain = sum(mult[t] * (blocks[t].r - ranks0[t]) for t in range(len(J)))
        return gain, sorted(chosen)

    h, steps = [], []
    cands = np.zeros((1, n), dtype=np.int64)
    d = dim = 0
    while cands.shape[0] and dim < total:
        t0 = time.time()
        gain, chosen = process(cands)
        h.append(gain)
        dim += gain
        piv = cands[chosen]
        front = np.vstack([piv, codes.canonical(codes.reflect(piv))[0]])
        _, first = np.unique(front @ codes.pow3, return_index=True)
        front = front[first]
        step = {'d': d, 'candidates': int(cands.shape[0]), 'h_d': int(gain),
                'frontier': int(front.shape[0]),
                'seconds': round(time.time() - t0, 1)}
        steps.append(step)
        if log:
            log(f'  n={n} p={p} d={d}: candidates {step["candidates"]}, '
                f'h_d={gain}, dim F_d={dim}, frontier {step["frontier"]}, '
                f'{step["seconds"]}s')
        d += 1
        cands = expand(front, codes)
    while h and h[-1] == 0:
        h.pop()
    return h, dim, steps


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('nmin', type=int)
    ap.add_argument('nmax', type=int)
    ap.add_argument('--primes', type=int, default=2,
                    help='number of primes per n (default 2)')
    ap.add_argument('--out', default=None,
                    help='JSON file to update with the results')
    args = ap.parse_args()

    def log(s):
        print(s, flush=True)

    log(f'# python {platform.python_version()}, numpy {np.__version__}, '
        f'{platform.machine()}, {platform.system()}')
    results = {}
    if args.out:
        try:
            with open(args.out) as fh:
                results = json.load(fh)
        except FileNotFoundError:
            pass
    ok = True
    for n in range(args.nmin, args.nmax + 1):
        N = len(rotation_data(n)[0])
        below = int((2 ** 51 / max(N, n)) ** 0.5)
        runs = []
        for p in primes_one_mod(n, below, args.primes):
            t0 = time.time()
            h, dim, steps = cycle_hilbert(n, p, log=log)
            el = time.time() - t0
            good = dim == (1 << n) - 1
            log(f'n={n} p={p}: H = {h}  (sum {dim}'
                f'{" = 2^n - 1" if good else ", WRONG TOTAL"}, {el:.1f}s)')
            ok &= good
            runs.append({'p': p, 'h': h, 'seconds': round(el, 1),
                         'steps': steps})
        agree = all(r['h'] == runs[0]['h'] for r in runs)
        if not agree:
            log(f'n={n}: the primes give different series')
            ok = False
        results[str(n)] = {'h': runs[0]['h'], 'agree': agree, 'runs': runs}
        if args.out:
            with open(args.out, 'w') as fh:
                json.dump(results, fh, indent=1, sort_keys=True)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
