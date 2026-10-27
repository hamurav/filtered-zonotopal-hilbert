#!/usr/bin/env python3
"""Certificate for Proposition 6.26: Conjecture 6.25 on all graphs with 4 to 7 vertices.

Cumulative upper-bound saturation, as in (6.21).  Let F_d be the span of
the monomials of degree <= d in the Y_v.  Proposition 6.24 gives, over Q,

    dim F_1 <= U_1 = n + 1,
    dim F_2 <= U_2 = U_1 + C(n+1, 2) - p(G) - eps_2,
    dim F_3 <= U_3 = U_2 + C(n+2, 3) - d_2(G) - e_22(G) - eps_3,

for every normalized nondegenerate f (eps_2 = 1 for cycles, eps_3 = 1 when
Delta(G) = 3 and beta != 2 alpha^2, else 0).  For an integer (or
p-integral rational) parameter, the rank of the same matrix modulo p is at
most its rank over Q.  If the cumulative ranks modulo p attain U_1, U_2
(part (i)) and U_3 (part (ii)), then equality holds over Q at that
parameter; since the generic rank is at least the rank at any parameter
and at most the bound, equality also holds for generic f, and the
coefficients h_2, h_3 equal their bounds in (6.20).

The script computes, for every connected simple graph with 4 to 7
vertices (networkx graph atlas), the cumulative ranks modulo
p = 1048573 at a fixed integer parameter and at the exponential
(coefficients 1/k!), and records which graphs are certified for parts (i)
and (ii) of Conjecture 6.25.  Expected: all 968 graphs of part (i) and
all 557 graphs of part (ii) for the integer parameter; at the exponential,
all 968 for part (i) and 555 of the 557 for part (ii).

A modular rank below the bound decides nothing over Q.  For every graph
not certified at the exponential, the script therefore recomputes the
Hilbert function at the exponential exactly, in rational arithmetic
(../../search/hilbert_direct.py).  Expected: the two graphs left over
(atlas G150 and G152) have dim F_3 = U_3 - 1 over Q as well, so the
exponential does not attain the cubic bound on them.

Usage:  python3 low_certificate.py [--jobs N]
Output: results/low_certificate.json, results/low_certificate.txt
"""

import argparse
import json
import random
import sys
import time
from fractions import Fraction
from math import comb, factorial
from multiprocessing import Pool

import os

import networkx as nx

from gz import PBIG
from lowdeg import low_hilbert

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'search'))
from hilbert_direct import hilbert_direct  # noqa: E402  (exact, over Q)


def invariants(G):
    deg = dict(G.degree())
    n, m = G.number_of_nodes(), G.number_of_edges()
    p = 0
    for v in G:
        if deg[v] == 1:
            p += 1
            (w,) = G[v]
            if deg[w] == 2:
                p += 1
    d2 = sum(1 for v in G if deg[v] == 2)
    e22 = sum(1 for a, b in G.edges() if deg[a] == 2 and deg[b] == 2)
    Delta = max(deg.values())
    is_cycle = all(d == 2 for d in deg.values())
    bridgeless = not any(True for _ in nx.bridges(G))
    return dict(n=n, m=m, p=p, d2=d2, e22=e22, Delta=Delta, cycle=is_cycle,
                bridgeless=bridgeless, cyclomatic=m - n + 1)


def bounds(inv, alpha, beta):
    n = inv['n']
    eps2 = 1 if (inv['cycle'] and alpha != 0) else 0
    eps3 = 1 if (inv['Delta'] == 3 and beta != 2 * alpha * alpha) else 0
    U1 = n + 1
    U2 = U1 + comb(n + 1, 2) - inv['p'] - eps2
    U3 = U2 + comb(n + 2, 3) - inv['d2'] - inv['e22'] - eps3
    return [U1, U2, U3]


def part_i(inv, G):
    return (inv['cyclomatic'] >= 1 and inv['n'] >= 4
            and not nx.is_isomorphic(G, nx.cycle_graph(4))
            and not nx.is_isomorphic(G, PAW))


def part_ii(inv, G):
    return (inv['bridgeless'] and inv['cyclomatic'] >= 3
            and not nx.is_isomorphic(G, nx.complete_graph(4)))


PAW = nx.Graph([(0, 1), (1, 2), (2, 0), (2, 3)])


def cumulative(h):
    out, s = [], 0
    for x in h:
        s += x
        out.append(s)
    return out


def job(args):
    idx, edges = args
    G = nx.Graph(edges)
    inv = invariants(G)
    n, D = inv['n'], inv['Delta']
    rnd = random.Random(1000 + idx)
    fg = {1: 1, **{k: rnd.randrange(1, PBIG) for k in range(2, D + 1)}}
    fe = {k: Fraction(1, factorial(k)) for k in range(1, D + 1)}
    out = dict(idx=idx, edges=edges, **inv, applicable_i=part_i(inv, G),
               applicable_ii=part_ii(inv, G), generic_parameter=fg)
    for tag, f in (('generic', fg), ('exp', fe)):
        alpha, beta = Fraction(f.get(2, 0)), Fraction(f.get(3, 0))
        U = bounds(inv, alpha, beta)
        r = cumulative(low_hilbert(n, edges, f, 3))[1:4]
        out[tag] = dict(bounds=U, ranks=r,
                        cert_i=r[0] == U[0] and r[1] == U[1],
                        cert_ii=r[0] == U[0] and r[2] == U[2])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--jobs', type=int, default=4)
    args = ap.parse_args()
    jobs = []
    for i, g in enumerate(nx.graph_atlas_g()):
        n = g.number_of_nodes()
        if 4 <= n <= 7 and nx.is_connected(g) and g.number_of_edges() >= n:
            jobs.append((i, sorted(tuple(sorted(e)) for e in g.edges())))
    t0 = time.time()
    with Pool(args.jobs) as pool:
        res = sorted(pool.map(job, jobs, chunksize=4), key=lambda r: r['idx'])
    lines = [f'# prime {PBIG}; {len(res)} connected non-tree graphs with '
             f'4 <= n <= 7; {time.time() - t0:.0f}s']
    ok = True
    for tag in ('generic', 'exp'):
        A = [r for r in res if r['applicable_i']]
        B = [r for r in res if r['applicable_ii']]
        ci = sum(r[tag]['cert_i'] for r in A)
        cii = sum(r[tag]['cert_ii'] for r in B)
        lines.append(f'{tag}: part (i) certified for {ci} of {len(A)} graphs; '
                     f'part (ii) certified for {cii} of {len(B)} graphs')
        for r in B:
            if not r[tag]['cert_ii']:
                lines.append(f'  part (ii) not certified ({tag}): atlas '
                             f'G{r["idx"]}, edges {r["edges"]}, ranks '
                             f'{r[tag]["ranks"]}, bounds {r[tag]["bounds"]}')
        for r in A:
            if not r[tag]['cert_i']:
                lines.append(f'  part (i) not certified ({tag}): atlas '
                             f'G{r["idx"]}, edges {r["edges"]}')
    counts_ok = (len([r for r in res if r['applicable_i']]) == 968 and
              len([r for r in res if r['applicable_ii']]) == 557)
    gen = all(r['generic']['cert_i'] for r in res if r['applicable_i']) and \
        all(r['generic']['cert_ii'] for r in res if r['applicable_ii'])
    lines.append('generic parameter: all applicable graphs certified'
                 if gen else 'generic parameter: SOME GRAPHS NOT CERTIFIED')
    # exact recomputation at the exponential where the modular test failed
    below = True
    for r in res:
        e = r['exp']
        if (r['applicable_i'] and not e['cert_i']) or \
                (r['applicable_ii'] and not e['cert_ii']):
            h = hilbert_direct(r['n'], [tuple(x) for x in r['edges']],
                               {k: Fraction(1, factorial(k))
                                for k in range(1, r['Delta'] + 1)})
            exact, U = cumulative(h)[1:4], e['bounds']
            e['exact_hilbert'], e['exact_ranks'] = h, exact
            if r['applicable_i'] and not e['cert_i']:
                below = below and exact[:2] != U[:2]
            if r['applicable_ii'] and not e['cert_ii']:
                below = below and (exact[0], exact[2]) != (U[0], U[2])
            lines.append(f'  exact over Q (exp): atlas G{r["idx"]}: h = {h}, '
                         f'cumulative ranks {exact}, bounds {e["bounds"]}')
    lines.append('exponential: every graph not certified modulo p is '
                 'below the bound over Q as well' if below else
                 'exponential: SOME EXACT VALUE ATTAINS THE BOUND')
    ok = counts_ok and gen and below
    text = '\n'.join(lines)
    print(text)
    with open('results/low_certificate.txt', 'w') as fh:
        fh.write(text + '\n')
    with open('results/low_certificate.json', 'w') as fh:
        json.dump(res, fh, default=str)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
