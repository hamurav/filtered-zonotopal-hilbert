"""Exhaustive search for Hilbert collisions among small trees and unicyclic graphs.

    python3 search.py trees     --nmax 10 [--direct] [--jobs 4]
    python3 search.py unicycles --nmax 8  [--jobs 4]

Every graph in the class is generated up to isomorphism (small_graphs.py) and
its filtered Hilbert function is computed exactly for each test function f.
Trees use the boundary-state calculus (hilbert_boundary.py); with --direct
every tree is also computed directly in the edge algebra (hilbert_direct.py)
and the two answers are compared.  Unicyclic graphs are computed directly.

Graphs with the same number of vertices whose Hilbert functions agree for all
test functions are reported as unseparated classes.  Two graphs separated by
some f are not universally Hilbert equivalent, so the unseparated classes
contain every universal collision.  Results are written to
results/<kind>_n<nmax><tag>.json and a summary to results/<kind>_n<nmax><tag>.txt.

The default test functions are u (graded), the four functions
u+u^2, u+u^2+u^3, u+u^2+2u^3, u+u^3 of an earlier search with Sage, and
u+u^4.  The last one is needed for trees with a vertex of degree four or
more: the first five functions do not separate, for example, the star with
four leaves from the tree with degrees (3,2,1,1,1).  To see this, run
    python3 search.py trees --nmax 10 --tag _without_u4 \
        --functions u u+u^2 u+u^2+u^3 u+u^2+2u^3 u+u^3
(recorded output: results/trees_n10_without_u4.txt).
"""
import argparse
import json
import os
import time
from collections import defaultdict
from fractions import Fraction
from multiprocessing import Pool

from small_graphs import trees, unicycles
from hilbert_direct import hilbert_direct, number_of_forests
from hilbert_boundary import LocalRule, hilbert_tree

FUNCTIONS = {
    'u': {1: 1},
    'u+u^2': {1: 1, 2: 1},
    'u+u^2+u^3': {1: 1, 2: 1, 3: 1},
    'u+u^2+2u^3': {1: 1, 2: 1, 3: 2},
    'u+u^3': {1: 1, 3: 1},
    'u+u^4': {1: 1, 4: 1},
}


def parse_function(s):
    """'u+u^2+2u^3' or '1,1,2' (coefficients of u, u^2, ...) -> {k: Fraction}."""
    if s in FUNCTIONS:
        return FUNCTIONS[s]
    coeffs = [Fraction(c) for c in s.split(',')]
    return {k + 1: c for k, c in enumerate(coeffs) if c}


def _direct_job(args):
    n, edges, f = args
    return hilbert_direct(n, edges, f)


def run(kind, nmax, names, direct, jobs, tag=''):
    fs = [FUNCTIONS.get(name) or parse_function(name) for name in names]
    graphs = list(trees(nmax)) if kind == 'trees' else list(unicycles(nmax))
    t0 = time.monotonic()
    profiles = {}
    if kind == 'trees':
        rules = [LocalRule(f) for f in fs]
        for idx, (n, edges) in enumerate(graphs):
            profiles[idx] = [hilbert_tree(n, edges, R) for R in rules]
    t_boundary = time.monotonic() - t0
    checked = 0
    if kind == 'unicycles' or direct:
        tasks = [(n, edges, f) for (n, edges) in graphs for f in fs]
        with Pool(jobs) as pool:
            out = pool.map(_direct_job, tasks, chunksize=1)
        for idx in range(len(graphs)):
            dh = out[idx * len(fs):(idx + 1) * len(fs)]
            if kind == 'trees':
                assert dh == profiles[idx], ('boundary and direct computations disagree', graphs[idx])
                checked += len(fs)
            else:
                profiles[idx] = dh
    t_total = time.monotonic() - t0
    for idx, (n, edges) in enumerate(graphs):
        forests = 2 ** len(edges) if kind == 'trees' else number_of_forests(n, edges)
        for name, h in zip(names, profiles[idx]):
            if name != 'u' and fs[names.index(name)].get(1):
                assert sum(h) == forests, ('filtration not exhaustive', edges, name)
    byn = defaultdict(lambda: defaultdict(list))
    for idx, (n, edges) in enumerate(graphs):
        byn[n][tuple(tuple(h) for h in profiles[idx])].append(idx)
    lines = [f'{kind}, n <= {nmax}; test functions: ' + ', '.join(names),
             f'Hilbert functions: {"boundary-state calculus" if kind == "trees" else "direct"}'
             + (f'; direct edge-algebra check of all {checked} values: agree' if checked else '')
             + f'; exact arithmetic; {t_total:.0f} s with {jobs} processes',
             '', ' n  graphs  classes  unseparated classes (sizes)']
    unseparated = []
    for n in sorted(byn):
        cls = list(byn[n].values())
        bad = [c for c in cls if len(c) > 1]
        unseparated += bad
        lines.append(f'{n:2d} {len(byn[n]) and sum(len(c) for c in cls):7d} {len(cls):8d}  '
                     + (', '.join(str(len(c)) for c in bad) if bad else 'none'))
    lines.append('')
    for c in unseparated:
        n = graphs[c[0]][0]
        lines.append(f'unseparated class with {n} vertices:')
        for idx in c:
            lines.append(f'  edges {graphs[idx][1]}')
        for name, h in zip(names, profiles[c[0]]):
            lines.append(f'  {name:12s} {h}')
    report = '\n'.join(lines)
    os.makedirs('results', exist_ok=True)
    base = f'results/{kind}_n{nmax}{tag}'
    with open(base + '.txt', 'w') as fh:
        fh.write(report + '\n')
    data = {'kind': kind, 'nmax': nmax, 'functions': {nm: {str(k): str(v) for k, v in f.items()} for nm, f in zip(names, fs)},
            'graphs': [{'n': n, 'edges': edges, 'hilbert': dict(zip(names, profiles[idx]))}
                       for idx, (n, edges) in enumerate(graphs)]}
    with open(base + '.json', 'w') as fh:
        json.dump(data, fh)
    print(report)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('kind', choices=['trees', 'unicycles'])
    ap.add_argument('--nmax', type=int, default=10)
    ap.add_argument('--functions', nargs='+', default=list(FUNCTIONS))
    ap.add_argument('--direct', action='store_true', help='trees: also compute every value directly and compare')
    ap.add_argument('--jobs', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--tag', default='', help='suffix for the names of the output files')
    a = ap.parse_args()
    run(a.kind, a.nmax, a.functions, a.direct, a.jobs, a.tag)
