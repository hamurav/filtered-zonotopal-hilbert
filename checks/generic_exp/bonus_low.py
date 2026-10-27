"""Low Hilbert coefficients of generic f and of the exponential for all connected simple graphs
on 4..7 vertices (networkx graph atlas).  h_0..h_4 for n <= 6, h_0..h_3 for n = 7.
Writes results/bonus_low.json with simple graph invariants for fitting."""
import json, random, sys, time
from fractions import Fraction as Fr
from math import factorial
from multiprocessing import Pool
import networkx as nx
from lowdeg import low_hilbert
from gz import PBIG


def pendant_vertices(G):
    """Vertices on pendant paths: leaves and the degree-2 vertices of paths leading to a leaf."""
    P = set()
    for v in G:
        if G.degree(v) == 1:
            prev, cur = None, v
            while True:
                P.add(cur)
                nbrs = [w for w in G[cur] if w != prev]
                if len(nbrs) != 1:
                    break
                nxt = nbrs[0]
                if G.degree(nxt) != 2 or nxt in P:
                    break
                prev, cur = cur, nxt
    return P


def job(args):
    idx, n, E = args
    G = nx.Graph(E)
    rnd = random.Random(idx)
    D = max(d for _, d in G.degree())
    dmax = 4 if n <= 6 else 3
    f1 = {1: 1, **{k: rnd.randrange(1, PBIG) for k in range(2, D + 1)}}
    f2 = {1: 1, **{k: rnd.randrange(1, PBIG) for k in range(2, D + 1)}}
    fe = {k: Fr(1, factorial(k)) for k in range(1, D + 1)}
    g1 = low_hilbert(n, E, f1, dmax)
    g2 = low_hilbert(n, E, f2, dmax)
    he = low_hilbert(n, E, fe, dmax)
    degs = sorted((d for _, d in G.degree()), reverse=True)
    return dict(idx=idx, n=n, edges=E, N=len(E), degs=degs, Delta=D,
                leaves=sum(1 for d in degs if d == 1), pendant=len(pendant_vertices(G)),
                bridges=len(list(nx.bridges(G))), triangles=sum(nx.triangles(G).values()) // 3,
                generic=g1, generic2=g2, exp=he)


if __name__ == '__main__':
    A = nx.graph_atlas_g()
    jobs = []
    for i, g in enumerate(A):
        n = g.number_of_nodes()
        if 4 <= n <= 7 and nx.is_connected(g):
            jobs.append((i, n, sorted(tuple(sorted(e)) for e in g.edges())))
    t = time.time()
    out = []
    with Pool(4) as pool:
        for r in pool.imap_unordered(job, jobs, chunksize=4):
            out.append(r)
            if len(out) % 100 == 0:
                print(len(out), 'done', '%.0fs' % (time.time() - t), flush=True)
    out.sort(key=lambda r: r['idx'])
    json.dump(out, open('results/bonus_low.json', 'w'))
    bad = [r['idx'] for r in out if r['generic'] != r['generic2']]
    print(len(out), 'graphs; generic draws disagree for', bad)
