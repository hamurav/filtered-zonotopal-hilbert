"""Exact checks for Section 6.10: Proposition 6.20, Table 11 and the rational points of
Remark 6.21 (graphs with five vertices and maximum degree four, and the 4-cycle with two adjacent
edges doubled).

Standard library only; all Hilbert functions are computed by exact rational arithmetic with
hilbert_direct.py.  For every graph of Table 11 the script computes

  * the exponential e^u - 1, truncated at u^4, and f5 = u + u^2 + u^4 (a generic column of
    Table B4 in ../tables/), and checks that the exponential has smaller cumulative dimensions than f5
    exactly for the graphs where Table 11 says it is exceptional;
  * the normalized functions f = u + u^2 + lam u^3 + mu u^4 at 20 rational points: two on each
    of the curves A, B, C, Q, L, four on W, the points E and K, two on the control line
    h_3 = 0 and two random points.

The expected value at a point is the generic series unless the point lies on an exceptional
locus of the graph (Remark 6.21), in which case it is the series recorded for that locus.  The
loci were found by complete scans modulo 101 and 61 (not part of this directory); this script
only confirms them exactly at the listed points.

Usage:  python3 exp_small.py [--jobs N] [graph ...]
"""
import argparse, sys, time
from fractions import Fraction as Fr
from hilbert_direct import hilbert_direct

GRAPHS = {
    'S5': dict(n=5, edges=[(0, 1), (0, 2), (0, 3), (0, 4)],
        generic=[1, 5, 10], exp=[1, 5, 10],
        loci={}),
    'cricket': dict(n=5, edges=[(0, 1), (1, 2), (2, 0), (0, 3), (0, 4)],
        generic=[1, 5, 13, 8, 1], exp=[1, 5, 13, 8, 1],
        loci={'K': [1, 5, 12, 9, 1]}),
    'butterfly': dict(n=5, edges=[(0, 1), (0, 2), (1, 2), (0, 3), (0, 4), (3, 4)],
        generic=[1, 5, 15, 21, 7], exp=[1, 5, 15, 21, 7],
        loci={'A': [1, 5, 15, 19, 9], 'K': [1, 5, 14, 20, 9]}),
    'dart': dict(n=5, edges=[(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (0, 4)],
        generic=[1, 5, 14, 25, 3], exp=[1, 5, 14, 23, 5],
        loci={'A': [1, 5, 14, 23, 5], 'C': [1, 5, 14, 23, 5], 'E': [1, 5, 14, 23, 5], 'K': [1, 5, 13, 23, 6]}),
    'fan': dict(n=5, edges=[(0, 1), (1, 2), (2, 3), (4, 0), (4, 1), (4, 2), (4, 3)],
        generic=[1, 5, 15, 33, 28], exp=[1, 5, 15, 33, 28],
        loci={'B': [1, 5, 15, 32, 29], 'A': [1, 5, 15, 31, 30], 'K': [1, 5, 14, 28, 33, 1]}),
    'K113': dict(n=5, edges=[(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4)],
        generic=[1, 5, 15, 32, 28], exp=[1, 5, 15, 32, 28],
        loci={'B': [1, 5, 15, 31, 29], 'A': [1, 5, 15, 30, 28, 2], 'K': [1, 5, 14, 27, 32, 2]}),
    'K4+leaf': dict(n=5, edges=[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (3, 4)],
        generic=[1, 5, 14, 29, 21, 5, 1], exp=[1, 5, 14, 28, 20, 6, 2],
        loci={'B': [1, 5, 14, 28, 22, 5, 1], 'C': [1, 5, 14, 28, 22, 5, 1], 'A': [1, 5, 14, 27, 23, 5, 1], 'Q': [1, 5, 14, 29, 19, 6, 2], 'E': [1, 5, 14, 28, 20, 6, 2], 'K': [1, 5, 13, 25, 24, 6, 2]}),
    'wheel': dict(n=5, edges=[(0, 1), (1, 2), (2, 3), (3, 0), (4, 0), (4, 1), (4, 2), (4, 3)],
        generic=[1, 5, 15, 35, 56, 18, 4], exp=[1, 5, 15, 35, 51, 20, 7],
        loci={'B': [1, 5, 15, 34, 56, 19, 4], 'A': [1, 5, 15, 34, 56, 19, 4], 'Q': [1, 5, 15, 35, 52, 19, 7], 'E': [1, 5, 15, 35, 51, 20, 7], 'K': [1, 5, 14, 30, 51, 26, 7]}),
    'houseX': dict(n=5, edges=[(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3), (0, 4), (1, 4)],
        generic=[1, 5, 15, 34, 53, 19, 1], exp=[1, 5, 15, 34, 49, 20, 4],
        loci={'B': [1, 5, 15, 33, 52, 21, 1], 'Q': [1, 5, 15, 34, 50, 20, 3], 'E': [1, 5, 15, 34, 49, 20, 4], 'K': [1, 5, 14, 29, 46, 28, 5]}),
    'K5-e': dict(n=5, edges=[(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4)],
        generic=[1, 5, 15, 35, 65, 64, 13], exp=[1, 5, 15, 35, 60, 59, 23],
        loci={'B': [1, 5, 15, 34, 62, 68, 13], 'Q': [1, 5, 15, 35, 64, 60, 18], 'E': [1, 5, 15, 35, 60, 59, 23], 'K': [1, 5, 14, 30, 53, 62, 33]}),
    'K5': dict(n=5, edges=[(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)],
        generic=[1, 5, 15, 35, 67, 91, 48, 15, 9, 4, 1], exp=[1, 5, 15, 35, 63, 84, 59, 15, 9, 4, 1],
        loci={'B': [1, 5, 15, 34, 64, 90, 53, 15, 9, 4, 1], 'W': [1, 5, 15, 35, 67, 86, 53, 15, 9, 4, 1], 'E': [1, 5, 15, 35, 63, 84, 59, 15, 9, 4, 1], 'K': [1, 5, 14, 30, 55, 80, 77, 15, 9, 4, 1]}),
    'C4adj2': dict(n=4, edges=[(0, 1), (0, 1), (1, 2), (1, 2), (2, 3), (3, 0)],
        generic=[1, 4, 10, 16, 1], exp=[1, 4, 10, 15, 2],
        loci={'L': [1, 4, 10, 15, 2], 'E': [1, 4, 10, 15, 2], 'K': [1, 4, 9, 15, 3]}),
}

# (curve, lam, mu) in the chart f = u + u^2 + lam u^3 + mu u^4
POINTS = [('A', '2', '1/3'), ('A', '2', '10/11'), ('B', '-17', '-90'), ('B', '14', '65'),
          ('C', '2/3', '3/5'), ('C', '2/3', '-11/3'), ('Q', '-13', '1465/8'), ('Q', '-29/7', '6561/392'),
          ('L', '13/2', '83/4'), ('L', '-5', '-39/2'), ('W', '8/3', '28/3'), ('W', '8/3', '13/3'),
          ('W', '26/15', '329/75'), ('W', '26/15', '269/75'), ('H', '30/7', '53/7'),
          ('H', '-33/10', '-38/5'), ('E', '2/3', '1/3'), ('K', '2', '5'),
          ('generic', '-25/4', '40/11'), ('generic', '34', '33/10')]


def loci_of(lam, mu):
    """The loci of Remark 6.21 through the point; h = f'(f^{-1}(u)) = 1 + h1 u + h2 u^2 + h3 u^3."""
    h1, h2, h3 = 2, 3 * lam - 2, 4 * (1 - 2 * lam + mu)
    out = set()
    if h2 == h1 ** 2: out.add('A')
    if h3 == 2 * h1 * h2 - h1 ** 3: out.add('B')
    if h2 == 0: out.add('C')
    if h1 * h3 == h2 ** 2: out.add('Q')
    if h3 == h1 * h2: out.add('L')
    if (h1**2 * h3**2 + 2 * h1 * h2**2 * h3 - 4 * h2**4 + 8 * h1**2 * h2**3 - 6 * h1**3 * h2 * h3
            - 3 * h1**4 * h2**2 + 2 * h1**5 * h3) == 0: out.add('W')
    if h2 == 0 and h3 == 0: out.add('E')
    if h2 == h1 ** 2 and h3 == h1 ** 3: out.add('K')
    return out


def expected(g, lam, mu):
    on = loci_of(lam, mu)
    loci = g['loci']
    for p in ('E', 'K'):
        if p in on and p in loci:
            return loci[p]
    vals = {tuple(loci[c]) for c in on if c in loci and c not in ('E', 'K')}
    if len(vals) > 1:
        raise ValueError('point on two exceptional curves with different series')
    return list(vals.pop()) if vals else g['generic']


def cum(h):
    out, s = [], 0
    for x in h:
        s += x
        out.append(s)
    return out


def smaller(h, k):
    """True if h has smaller cumulative dimension than k in some degree and never larger."""
    a, b = cum(h), cum(k)
    m = max(len(a), len(b))
    a += [a[-1]] * (m - len(a))
    b += [b[-1]] * (m - len(b))
    return all(x <= y for x, y in zip(a, b)) and a != b


def job(args):
    name, lab, f = args
    g = GRAPHS[name]
    return name, lab, hilbert_direct(g['n'], g['edges'], f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=1)
    ap.add_argument('graphs', nargs='*')
    a = ap.parse_args()
    names = a.graphs or list(GRAPHS)
    jobs = []
    for name in names:
        jobs.append((name, 'exp', {k: Fr(1, [1, 1, 2, 6, 24][k]) for k in range(1, 5)}))
        jobs.append((name, 'f5', {1: 1, 2: 1, 4: 1}))
        for lab, l, m in POINTS:
            jobs.append((name, (lab, l, m), {1: 1, 2: 1, 3: Fr(l), 4: Fr(m)}))
    t = time.time()
    if a.jobs > 1:
        from multiprocessing import Pool
        with Pool(a.jobs) as pool:
            res = pool.map(job, jobs)
    else:
        res = [job(j) for j in jobs]
    bad = 0
    by = {}
    for name, lab, h in res:
        by.setdefault(name, {})[lab if isinstance(lab, str) else lab] = h
    for name in names:
        g, r = GRAPHS[name], by[name]
        print('==', name)
        ok_exp = r['exp'] == g['exp']
        exc = g['exp'] != g['generic']
        ok_cmp = smaller(r['exp'], r['f5']) if exc else r['exp'] == r['f5']
        ok_f5 = r['f5'] == g['generic']
        print(f"   e^u-1   {r['exp']}  {'ok' if ok_exp else 'MISMATCH'}")
        print(f"   f5      {r['f5']}  {'ok' if ok_f5 else 'MISMATCH'}")
        print(f"   exponential {'exceptional (smaller than f5)' if exc else 'same as f5'}: {'ok' if ok_cmp else 'FAILED'}")
        bad += (not ok_exp) + (not ok_cmp) + (not ok_f5)
        for lab, l, m in POINTS:
            h = r[(lab, l, m)]
            e = expected(g, Fr(l), Fr(m))
            flag = 'ok' if h == e else 'MISMATCH (expected %s)' % e
            bad += h != e
            print(f"   {lab:8s} lam={l:7s} mu={m:9s} {h}  {flag}")
    print('checked %d Hilbert functions in %.0f s; %d mismatches' % (len(res), time.time() - t, bad))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
