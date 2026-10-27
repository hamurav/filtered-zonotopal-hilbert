"""Attribute every exceptional point of the GF(p) scans (chart alpha = 1) to the candidate loci
A, B, C, Q, L, W and the points E, K (see verify_exact.py for definitions), and check that each
locus is either entirely exceptional with one series (apart from its special points) or not at all.
Usage: python3 classify.py p  ->  prints a per-graph summary, writes results/classify_p{p}.json"""
import glob, json, sys
from collections import defaultdict

p = int(sys.argv[1])
inv = lambda x: pow(x % p, p - 2, p)


def hcoords(l, m):
    return 2, (3 * l - 2) % p, 4 * (1 - 2 * l + m) % p


def loci(l, m):
    h1, h2, h3 = hcoords(l, m)
    out = set()
    if l % p == 2: out.add('A')
    if (m - 5 * l + 5) % p == 0: out.add('B')
    if (3 * l - 2) % p == 0: out.add('C')
    if (h1 * h3 - h2 * h2) % p == 0: out.add('Q')
    if (h3 - h1 * h2) % p == 0: out.add('L')
    if (h1**2 * h3**2 + 2 * h1 * h2**2 * h3 - 4 * h2**4 + 8 * h1**2 * h2**3 - 6 * h1**3 * h2 * h3
            - 3 * h1**4 * h2**2 + 2 * h1**5 * h3) % p == 0: out.add('W')
    if h2 == 0 and h3 == 0: out.add('E')
    if l % p == 2 and m % p == 5: out.add('K')
    return out


summary = {}
for fn in sorted(glob.glob(f'results/scan_p{p}_*.json')):
    d = json.load(open(fn))
    pts = defaultdict(list)
    alpha0 = {}
    for h, fl in d['exceptional'].items():
        for f in fl:
            f = dict(tuple(t) for t in f)
            if 2 in f and 4 in f:
                pts[h].append((f[3], f.get(4, 0)))
            elif 2 in f:
                pts[h].append((f[3], None))
            else:
                alpha0.setdefault(h, []).append(f)
    lines = []
    for h, P in pts.items():
        pats = defaultdict(int)
        for l, m in P:
            if m is None:   # Delta = 3 graph
                pats['lam=%d' % l] += 1
            else:
                pats[''.join(sorted(loci(l, m))) or 'NONE'] += 1
        lines.append((h, dict(pats)))
    # exhaustiveness per locus: for Delta=4, count points of each curve that are exceptional
    summary[d['name']] = dict(generic=d['generic'], alpha_ne_0=lines,
                              alpha_0={h: [sorted(f.items()) for f in fl] for h, fl in alpha0.items()})
    print('==', d['name'], 'generic', d['generic'])
    for h, pats in lines:
        print('   ', h, pats)
    for h, fl in alpha0.items():
        print('    alpha=0:', h, 'at', [sorted(f.items()) for f in fl][:3], '...' if len(fl) > 3 else '')
json.dump(summary, open(f'results/classify_p{p}.json', 'w'), indent=1)
