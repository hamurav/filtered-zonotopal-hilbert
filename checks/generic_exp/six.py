"""Generic and exponential Hilbert series for all connected simple graphs on NV vertices
with at most MAXE edges (mod the prime 1048573; generic = two random draws that must agree).
Usage: python3 six.py NV MAXE   ->  results/six_{NV}_{MAXE}.json"""
import json, random, sys, time
from fractions import Fraction as Fr
from math import factorial
from multiprocessing import Pool
from gz import GraphAlg, connected_graphs, PBIG


def job(args):
    nv, E, seed = args
    rnd = random.Random(seed)
    t = time.time()
    A = GraphAlg(nv, E)
    D = A.Delta
    gens = []
    for trial in range(2):
        f = {1: 1}
        for k in range(2, D + 1):
            f[k] = rnd.randrange(1, PBIG)
        gens.append(A.hilbert(f))
    hexp = A.hilbert({k: Fr(1, factorial(k)) for k in range(1, D + 1)})
    return dict(n=nv, edges=E, Delta=D, forests=A.D, generic=gens[0], generic2=gens[1], exp=hexp,
                time=time.time() - t)


if __name__ == '__main__':
    nv, maxe = int(sys.argv[1]), int(sys.argv[2])
    gs = [E for E in connected_graphs(nv) if len(E) <= maxe]
    gs.sort(key=len, reverse=True)
    print(len(gs), 'graphs', flush=True)
    out = []
    with Pool(4) as pool:
        # imap keeps the input order, so that the output does not depend on timing
        for r in pool.imap(job, [(nv, E, 1000 + i) for i, E in enumerate(gs)]):
            out.append(r)
            flag = '' if r['generic'] == r['exp'] else '  EXP NON-GENERIC'
            ok = '' if r['generic'] == r['generic2'] else '  (generic draws disagree!)'
            print(len(r['edges']), r['edges'], r['generic'], r['exp'], flag + ok, '%.0fs' % r['time'], flush=True)
    json.dump(out, open(f'results/six_{nv}_{maxe}.json', 'w'))
