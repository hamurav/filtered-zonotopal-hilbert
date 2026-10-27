"""Degenerate functions (Section 3.7 of the paper, [NSh, Problem 4]).

    python3 degenerate_check.py trees --nmax 12     # Theorem 3.20 for all trees with at most nmax vertices
    python3 degenerate_check.py lemma --nmax 9      # Lemma 3.23 on random families of exponent sets
    python3 degenerate_check.py cycles              # Proposition 3.29 (cycles) and Remark 3.30

For degenerate f (no linear term) the elements f(X_v) generate a proper subalgebra of C_G, filtered by
word length.  A partial score vector of G is the out-degree vector of a partial orientation
(Definition 3.19), and s_k(G; t) counts the partial score vectors with all entries divisible by k, by
weight/k.  Theorem 3.20: for a forest and f of order k >= 2 the Hilbert series is s_k.  For every graph
s_k bounds the dimensions of the filtration spaces from above (Lemma 3.22).  Proposition 3.29: for the
cycle C_n and f of order two the Hilbert series is s_2(C_n; t), minus t^(n/2) for even n.

The checks are not used in the proofs.  `trees` compares the Hilbert functions, computed directly in
the edge algebra, with s_k at eight functions, and compares the recursion of Remark 3.27 with the
enumeration of score vectors.  `lemma` tests Lemma 3.23 on random families of exponent sets.  `cycles`
checks Proposition 3.29 for n <= 12 and produces the data of Remark 3.30.  Each run also writes its
output to results/degenerate_<command>.txt.

Standard library only; exact integer arithmetic (hilbert_direct.py).
"""
import argparse
import os
import random
import sys
import time
from collections import Counter
from itertools import combinations, permutations, product
from math import comb

from small_graphs import trees, unicycles
from hilbert_direct import hilbert_direct, multiply, EchelonSpace

FUNCTIONS = {
    'u^2': (2, {2: 1}), 'u^2+u^3': (2, {2: 1, 3: 1}), 'u^2-2u^3+5u^4': (2, {2: 1, 3: -2, 4: 5}),
    'u^3': (3, {3: 1}), 'u^3+u^4': (3, {3: 1, 4: 1}), '2u^3-u^5': (3, {3: 2, 5: -1}),
    'u^4': (4, {4: 1}), 'u^4+3u^5': (4, {4: 1, 5: 3}),
}


def partial_scores(n, edges):
    vecs = {tuple([0] * n)}
    for (a, b) in edges:
        new = set()
        for v in vecs:
            new.add(v)
            if a != b:
                w = list(v); w[a] += 1; new.add(tuple(w))
                w = list(v); w[b] += 1; new.add(tuple(w))
        vecs = new
    return vecs


def divisible_series(vecs, k):
    """Coefficients of s_k(G; t) from the set of all partial score vectors."""
    c = Counter(sum(v) // k for v in vecs if all(x % k == 0 for x in v))
    return [c[m] for m in range(max(c) + 1)]


def _padd(p, q):
    r = [0] * max(len(p), len(q))
    for i, c in enumerate(p): r[i] += c
    for i, c in enumerate(q): r[i] += c
    return r


def _pmul(p, q):
    if not p or not q: return []
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): r[i + j] += a * b
    return r


def divisible_series_tree(n, edges, k, root=0):
    """s_k of a tree by the recursion of Remark 3.27 (no enumeration of score vectors).

    Below a vertex v with parent edge g: A_v counts vectors realizable without g, B_v those realizable
    only with g given to v.  If j children are of type A, then v may take any a_v <= j, and the vector
    is of type B exactly when a_v = j + 1."""
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b); adj[b].append(a)
    order, parent = [root], {root: None}
    for v in order:
        for w in adj[v]:
            if w not in parent:
                parent[w] = v; order.append(w)
    A, B = {}, {}
    for v in reversed(order):
        Q = [[1]]  # Q[j]: children configurations with exactly j children of type A
        for c in adj[v]:
            if c == parent[v]: continue
            new = [[] for _ in range(len(Q) + 1)]
            for j, p in enumerate(Q):
                new[j] = _padd(new[j], _pmul(p, B[c]))
                new[j + 1] = _padd(new[j + 1], _pmul(p, A[c]))
            Q = new
        A[v], B[v] = [], []
        for j, p in enumerate(Q):
            A[v] = _padd(A[v], _pmul(p, [1] * (j // k + 1)))
            if (j + 1) % k == 0:
                B[v] = _padd(B[v], _pmul(p, [0] * ((j + 1) // k) + [1]))
    h = A[root]
    while len(h) > 1 and h[-1] == 0: h.pop()
    return h


def P(n, edges, a):
    """prod_v e_{a_v}(signed star of v) in Phi_G, as {edge bitmask: coefficient} (signs as in hilbert_direct)."""
    res = {0: 1}
    for v in range(n):
        if a[v] == 0: continue
        inc = [(1 << i, -1 if v == x else 1) for i, (x, y) in enumerate(edges) if v in (x, y) and x != y]
        ev = {}
        for sub in combinations(inc, a[v]):
            m, s = 0, 1
            for bit, t in sub: m |= bit; s *= t
            ev[m] = ev.get(m, 0) + s
        res = multiply(res, {m: c for m, c in ev.items() if c})
    return {m: c for m, c in res.items() if c}


def cmd_trees(nmax):
    t0 = time.monotonic()
    count, bad = Counter(), Counter()
    for n, e in trees(nmax):
        if not e: continue
        vecs = partial_scores(n, e)
        pred = {k: divisible_series(vecs, k) for k in (2, 3, 4)}
        for k in pred:
            assert divisible_series_tree(n, e, k) == pred[k], (n, e, k)
        for name, (k, f) in FUNCTIONS.items():
            count[name] += 1
            if hilbert_direct(n, e, f) != pred[k]:
                bad[name] += 1
                print('MISMATCH', name, n, e)
    print(f'trees with 2 to {nmax} vertices: {count["u^2"]}')
    print('the recursion of Remark 3.27 equals the enumeration of score vectors for k = 2, 3, 4')
    for name in FUNCTIONS:
        print(f'  {name:14s} Hilbert function = s_k for {count[name] - bad[name]} of {count[name]} trees')
    print(f'{time.monotonic() - t0:.0f}s')


def no_consecutive(d, rng):
    s, x = [], 0
    while x <= d:
        if rng.random() < 0.5:
            s.append(x); x += 2
        else:
            x += 1
    return s or [0]


def cmd_lemma(nmax, trials):
    rng = random.Random(1)
    tested = dependent = 0
    for n, e in trees(nmax):
        if not e: continue
        deg = [sum(v in ed for ed in e) for v in range(n)]
        for _ in range(trials):
            A = [no_consecutive(deg[v], rng) for v in range(n)]
            space, nonzero = EchelonSpace(), 0
            for a in product(*A):
                p = P(n, e, a)
                if p:
                    nonzero += 1
                    space.insert(p)
            tested += 1
            dependent += len(space) != nonzero
    print(f'{tested} random families (A_v without two consecutive integers) on trees with at most {nmax} vertices;'
          f' families with a linear dependence: {dependent}')
    n, e = 3, [(0, 1), (1, 2)]
    space, nonzero = EchelonSpace(), 0
    for a in product(range(2), range(3), range(2)):
        p = P(n, e, a)
        if p:
            nonzero += 1; space.insert(p)
    print(f'the hypothesis is needed: path with 2 edges, all exponents: {nonzero} nonzero P_a, rank {len(space)}')


def bipartite(n, edges):
    col = [None] * n
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b); adj[b].append(a)
    for s in range(n):
        if col[s] is not None: continue
        col[s] = 0; stack = [s]
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if col[w] is None:
                    col[w] = 1 - col[v]; stack.append(w)
                elif col[w] == col[v]:
                    return False
    return True


def connected_with_cycle(nmax):
    """Connected simple graphs with 3 to nmax vertices that are not trees, up to isomorphism (brute force)."""
    for n in range(3, nmax + 1):
        pairs = list(combinations(range(n), 2))
        perms = list(permutations(range(n)))
        seen = set()
        for m in range(n, len(pairs) + 1):
            for es in combinations(pairs, m):
                comp, stack = {0}, [0]
                while stack:
                    v = stack.pop()
                    for a, b in es:
                        for x, y in ((a, b), (b, a)):
                            if x == v and y not in comp:
                                comp.add(y); stack.append(y)
                if len(comp) < n: continue
                key = min(tuple(sorted(tuple(sorted((p[a], p[b]))) for a, b in es)) for p in perms)
                if key not in seen:
                    seen.add(key)
                    yield n, list(es)


MULTIGRAPHS = {
    'double edge': (2, [(0, 1), (0, 1)]),
    'theta graph': (2, [(0, 1), (0, 1), (0, 1)]),
    'star K13, one edge doubled': (4, [(0, 1), (0, 1), (0, 2), (0, 3)]),
    'path P3, both edges doubled': (3, [(0, 1), (0, 1), (1, 2), (1, 2)]),
    'triangle, one edge doubled': (3, [(0, 1), (0, 1), (1, 2), (0, 2)]),
    'triangle, all edges doubled': (3, [(0, 1), (0, 1), (1, 2), (1, 2), (0, 2), (0, 2)]),
    'C4, opposite edges doubled': (4, [(0, 1), (0, 1), (1, 2), (2, 3), (2, 3), (3, 0)]),
    'C4, adjacent edges doubled': (4, [(0, 1), (0, 1), (1, 2), (1, 2), (2, 3), (3, 0)]),
}


def compare(n, e, stats, label=None):
    vecs = partial_scores(n, e)
    bip = 'bipartite' if bipartite(n, e) else 'not bipartite'
    for k, fs in ((2, ('u^2', 'u^2+u^3')), (3, ('u^3', 'u^3+u^4'))):
        bound = divisible_series(vecs, k)
        hs = [hilbert_direct(n, e, FUNCTIONS[name][1]) for name in fs]
        sharp = hs[0] == bound
        same = hs[0] == hs[1]
        stats[bip, k, sharp, same] += 1
        if label and (not sharp or not same):
            print(f'  {label:30s} k={k}  {fs[0]}: {hs[0]}  {fs[1]}: {hs[1]}  s_k: {bound}')


def show(stats):
    for (bip, k, sharp, same), c in sorted(stats.items()):
        print(f'  {bip:13s} k={k}  ' + ('h = s_k ' if sharp else 'h < s_k ')
              + ('' if same else '  (and u^k, u^k+u^(k+1) differ)') + f': {c}')


def cmd_cycles(nuni, nsimple):
    t0 = time.monotonic()
    print('Proposition 3.29: cycles C_n; formula sum_m n/(n-m) binom(n-m, m) t^m - [n even] t^(n/2)')
    lucas = [2, 1]
    for n in range(2, 13):
        lucas.append(lucas[-1] + lucas[-2])
    for n in range(3, 13):
        e = [(i, (i + 1) % n) for i in range(n)]
        formula = [n * comb(n - m, m) // (n - m) for m in range(n // 2 + 1)]
        if n % 2 == 0:
            formula[-1] -= 1
        h = hilbert_direct(n, e, {2: 1})
        ok = (h == formula == hilbert_direct(n, e, {2: 1, 3: 5}) and hilbert_direct(n, e, {3: 1}) == [1]
              and sum(h) == lucas[n] - (n % 2 == 0))
        b = divisible_series(partial_scores(n, e), 2)
        print(f'  C{n:<2d} u^2: {h}, sum {sum(h)};  s_2 = {b};  formula, u^2+5u^3 and u^3 ' + ('agree' if ok else 'DISAGREE'))
        assert ok
    stats = Counter()
    for n, e in unicycles(nuni):
        compare(n, e, stats)
    print(f'connected unicyclic simple graphs with at most {nuni} vertices:')
    show(stats)
    stats = Counter()
    for n, e in connected_with_cycle(nsimple):
        compare(n, e, stats, label=None)
    print(f'connected simple graphs with at most {nsimple} vertices, not trees:')
    show(stats)
    for n, e in connected_with_cycle(nsimple):
        a, b = hilbert_direct(n, e, {2: 1}), hilbert_direct(n, e, {2: 1, 3: 1})
        if a != b:
            print(f'  u^2 and u^2+u^3 differ on {e}: {a} and {b}')
    print('multigraphs (only cases with h < s_k or a difference between the two functions are listed):')
    stats = Counter()
    for label, (n, e) in MULTIGRAPHS.items():
        compare(n, e, stats, label)
    show(stats)
    print(f'{time.monotonic() - t0:.0f}s')


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for st in self.streams:
            st.write(text)

    def flush(self):
        for st in self.streams:
            st.flush()


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('trees'); p.add_argument('--nmax', type=int, default=12)
    p = sub.add_parser('lemma'); p.add_argument('--nmax', type=int, default=9); p.add_argument('--trials', type=int, default=5)
    p = sub.add_parser('cycles'); p.add_argument('--nuni', type=int, default=8); p.add_argument('--nsimple', type=int, default=5)
    args = ap.parse_args()
    os.makedirs('results', exist_ok=True)
    with open(os.path.join('results', f'degenerate_{args.cmd}.txt'), 'w') as fh:
        sys.stdout = Tee(sys.__stdout__, fh)
        if args.cmd == 'trees':
            cmd_trees(args.nmax)
        elif args.cmd == 'lemma':
            cmd_lemma(args.nmax, args.trials)
        else:
            cmd_cycles(args.nuni, args.nsimple)
        sys.stdout = sys.__stdout__
