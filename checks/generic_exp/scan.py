"""Complete scan of the normalized parameter space over GF(p).

Delta = 3:  alpha != 0:  f = u + u^2 + lam u^3,           lam in GF(p)
            alpha  = 0:  f = u + u^3, f = u
Delta = 4:  alpha != 0:  f = u + u^2 + lam u^3 + mu u^4,   (lam, mu) in GF(p)^2
            alpha  = 0:  f = u + b u^3 + nu u^4 with b in {1, a non-residue}, nu in GF(p);
                         f = u + u^4, f = u
Every point with a series different from the char-0 generic one is recorded.
Usage: python3 scan.py p name1 [name2 ...]   (names as in graphs_list; 'ALL4', 'ALL3' allowed)
Writes results/scan_p{p}_{slug}.json"""
import json, re, sys, time
from multiprocessing import Pool
from gz import GraphAlg, number_of_forests
from graphs_list import all_graphs

GR = all_graphs()
PH1 = json.load(open('results/phase1.json'))


def slug(name):
    return re.sub(r'[^A-Za-z0-9]+', '_', name).strip('_')


_A = {}


def get_alg(name, p):
    key = (name, p)
    if key not in _A:
        n, E = GR[name]
        _A[key] = GraphAlg(n, E, p)
    return _A[key]


def work(args):
    name, p, pts = args
    A = get_alg(name, p)
    out = []
    for f in pts:
        out.append((f, A.hilbert({k: v for k, v in f})))
    return out


def nonresidue(p):
    for r in range(2, p):
        if pow(r, (p - 1) // 2, p) == p - 1:
            return r


def run(name, p, pool):
    n, E = GR[name]
    A = get_alg(name, p)
    assert A.D == number_of_forests(n, E), ('bad prime for C_G', name, p)
    gen = PH1[name]['generic']
    D = A.Delta
    pts = []
    if D == 3:
        pts += [((1, 1), (2, 1), (3, lam)) for lam in range(p)]
        pts += [((1, 1), (3, 1)), ((1, 1),)]
    elif D == 4:
        pts += [((1, 1), (2, 1), (3, lam), (4, mu)) for lam in range(p) for mu in range(p)]
        r = nonresidue(p)
        pts += [((1, 1), (3, b), (4, nu)) for b in (1, r) for nu in range(p)]
        pts += [((1, 1), (4, 1)), ((1, 1),)]
    else:
        pts += [((1, 1), (2, 1)), ((1, 1),)]
    chunks = [(name, p, pts[i::16]) for i in range(16)]
    res = []
    for part in pool.map(work, chunks):
        res += part
    exc = {}
    for f, h in res:
        if h != gen:
            exc.setdefault(str(h), []).append(f)
    return gen, exc, len(pts)


if __name__ == '__main__':
    p = int(sys.argv[1])
    names = []
    for a in sys.argv[2:]:
        if a in ('ALL3', 'ALL4'):
            names += [k for k, v in PH1.items() if v['Delta'] == int(a[-1])]
        else:
            names.append(a)
    with Pool(4) as pool:
        for name in names:
            t = time.time()
            gen, exc, npts = run(name, p, pool)
            json.dump(dict(name=name, p=p, generic=gen, exceptional=exc, points=npts),
                      open(f'results/scan_p{p}_{slug(name)}.json', 'w'))
            print(f'{name}: p={p}, {npts} points, generic {gen}; exceptional series:', flush=True)
            for h, fl in sorted(exc.items(), key=lambda x: -len(x[1])):
                print(f'   {h}: {len(fl)} points, e.g. {fl[:4]}', flush=True)
            print(f'   ({time.time() - t:.0f}s)', flush=True)
