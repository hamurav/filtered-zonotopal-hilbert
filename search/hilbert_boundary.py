"""Filtered Hilbert functions of trees by the boundary-state calculus.

Standard library only; exact arithmetic.  This is an independent
implementation of the calculus of Section 4 of the paper (Lemma 4.1,
Lemma 4.3, the closing evaluation and the local vertex rule, collected in
Theorem 4.8):

* A marked leaf edge carries boundary states [p, q] (with multiplicity).
  A single edge whose far end is open has the one state [0, 1].
* At a vertex with incoming branch states [p_j, q_j] (j = 1..r) and an
  outgoing edge z, the local filtration of
      Q[x_1, ..., x_r, z] / (x_j^2, z^2),   y = f(z - x_1 - ... - x_r),
  is spanned by x_S y^k (k = 0..r+1) in relative degree sum_{j in S} d_j + k,
  with d_j = q_j - p_j; the outgoing states are read off from the image
  filtration Q_i = P_i mod z and the kernel filtration R_j = {u : z u in P_j}
  through dim(Q_i cap R_j) = rk(z P_i) + rk(P_j) - rk(P_j + z P_i) and the
  second difference (4.4), and then shifted by sum_j p_j.
* Two marked graphs closed along their common edge contribute
  t^(p+p') + t^min(q+p', p+q') for each pair of states.

Local tables depend only on the sorted gaps (d_j) and on f, and are cached.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import factorial

from hilbert_direct import EchelonSpace, integer_coefficients


def _mul(x, y):
    out = {}
    for i, a in x.items():
        for j, b in y.items():
            if not i & j:
                out[i | j] = out.get(i | j, 0) + a * b
    return {k: c for k, c in out.items() if c}


def _rank(vectors):
    space = EchelonSpace()
    for v in vectors:
        if v:
            space.insert(v)
    return len(space)


class LocalRule:
    """Cached local vertex tables for one polynomial f (integer coefficients after scaling)."""

    def __init__(self, f):
        self.f = integer_coefficients(f)
        self.cache = {}

    def y(self, r):
        zbit = 1 << r
        lin = {1 << j: -1 for j in range(r)}
        lin[zbit] = 1
        power, out = {0: 1}, {}
        for k in range(1, r + 2):
            power = _mul(power, lin)
            a = self.f.get(k, 0)
            for m, c in power.items():
                out[m] = out.get(m, 0) + a * c
        return {m: c for m, c in out.items() if c}

    def __call__(self, gaps):
        gaps = tuple(sorted(gaps))
        if gaps not in self.cache:
            self.cache[gaps] = self._compute(gaps)
        return self.cache[gaps]

    def _compute(self, gaps):
        r = len(gaps)
        zbit = 1 << r
        y = self.y(r)
        powers = [{0: 1}]
        for _ in range(r + 1):
            powers.append(_mul(powers[-1], y))
        columns = defaultdict(list)
        for S in range(1 << r):
            base = sum(gaps[j] for j in range(r) if S >> j & 1)
            for k, w in enumerate(powers):
                c = _mul({S: 1}, w)
                if c:
                    columns[base + k].append(c)
        # Echelon form of P_i with the smallest monomial as pivot, built level by
        # level.  A row with pivot < z is z-free in its leading term: these rows,
        # reduced mod z, form a basis of Q_i.  A row with pivot >= z lies in (z):
        # these rows, divided by z, form a basis of R_i.
        qvecs, svecs = [], []          # (level, vector in U)
        Pspace = EchelonSpace()
        for i in sorted(columns):
            for c in columns[i]:
                row = Pspace.insert(c)
                if row is None:
                    continue
                if min(row) < zbit:
                    qvecs.append((i, {m: a for m, a in row.items() if not m & zbit}))
                else:
                    svecs.append((i, {m ^ zbit: a for m, a in row.items()}))
        assert len(qvecs) == len(svecs) == zbit, ('filtration not exhaustive', gaps)
        # The s-vectors have distinct smallest monomials, so they form a
        # triangular basis of U.  Write each q-vector in this basis, then
        # eliminate the q-vectors in order of Q-level with the s-vector of
        # largest R-index as pivot.  A vector entering Q at level p whose pivot
        # has R-level q is a basis vector of state [p, q]: this is the
        # second difference (4.4) of dim(Q_p cap R_q).
        spos = {min(v): k for k, (_, v) in enumerate(svecs)}
        slev = [lev for lev, _ in svecs]
        coords = []
        for p, q in qvecs:
            q = {m: Fraction(a) for m, a in q.items()}
            c = {}
            while q:
                m = min(q)
                k = spos[m]
                s = svecs[k][1]
                a = q[m] / s[m]
                c[k] = a
                for mm, b in s.items():
                    d = q.get(mm, 0) - a * b
                    if d:
                        q[mm] = d
                    else:
                        q.pop(mm, None)
            coords.append((p, c))
        rows = {}
        states = Counter()
        for p, c in coords:
            while c:
                k = max(c)
                row = rows.get(k)
                if row is None:
                    rows[k] = c
                    states[p, slev[k]] += 1
                    break
                a = c[k] / row[k]
                for kk, b in row.items():
                    d = c.get(kk, 0) - a * b
                    if d:
                        c[kk] = d
                    else:
                        c.pop(kk, None)
            else:
                raise AssertionError('Q-vectors dependent')
        assert sum(states.values()) == zbit
        return states


def close(S, T):
    H = Counter()
    for (p, q), a in S.items():
        for (pp, qq), b in T.items():
            H[p + pp] += a * b
            H[min(q + pp, p + qq)] += a * b
    return [H[d] for d in range(max(H) + 1)]


def hilbert_tree(n, edges, rule):
    """Hilbert function of a tree with vertices 0..n-1, for the LocalRule `rule`."""
    if not edges:
        return [1]
    adj = defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    memo = {}

    def message(v, parent):
        key = (v, parent)
        if key not in memo:
            incoming = [message(w, v) for w in adj[v] if w != parent]
            out = Counter()
            for combo in product(*(list(S.items()) for S in incoming)):
                base = sum(pq[0] for pq, _ in combo)
                weight = 1
                for _, m in combo:
                    weight *= m
                for (p, q), m in rule(tuple(pq[1] - pq[0] for pq, _ in combo)).items():
                    out[p + base, q + base] += weight * m
            memo[key] = out
        return memo[key]
    u, v = edges[0]
    return close(message(u, v), message(v, u))


if __name__ == '__main__':
    import time
    rule = LocalRule({1: 1, 2: 1})
    print('leaf state', dict(rule(())), 'expected {(0, 1): 1}')
    print('path b of Example 4.2', dict(rule((1,))), 'expected states [0,2] and [1,1]')
    D6 = [(0, 1), (0, 2), (0, 3), (3, 4), (4, 5)]
    E6 = [(0, 1), (0, 2), (2, 3), (0, 4), (4, 5)]
    for f in [{1: 1, 3: 1}, {1: 1, 2: 1}, {1: 1, 2: 1, 3: 2}]:
        R = LocalRule(f)
        print(f, hilbert_tree(6, D6, R), hilbert_tree(6, E6, R))
    t = time.monotonic()
    T3a = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 6), (2, 7), (3, 8), (8, 9)]
    print('T3a u+u^2', hilbert_tree(10, T3a, rule), f'{time.monotonic() - t:.3f}s')
