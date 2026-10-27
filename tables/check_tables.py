"""Recompute every entry of the tables B1-B4 exactly (Fraction arithmetic).

Graphs are read from the TikZ macros in the table file itself: vertices are
the \\node[vertex] names, edges are the -- and to[...] segments of \\draw paths.
Hilbert series are computed with direct_hilbert.py from certificates/ (exact).

Usage:  python3 check_tables.py TABLES.tex [--certificates DIR] [--mod P] [--json FILE]
  TABLES.tex       e.g. tables345_clean.tex
Handles the 4-column (old), 5-column (complete cubic profile) and 9-column tables.
  --certificates   directory containing direct_hilbert.py
                   (default: the directory of this script)
  --mod P          compute modulo the prime P instead (fast sample check)
  --json FILE      also write every entry (graph, edge list, f, Hilbert
                   function) to FILE and to the .csv file of the same name
The four tables of the source are Tables B1-B4 (formerly Appendix B of the
paper, Tables 12-15 of its version of 21 October 2026).
"""
import argparse, os, re, sys, time
from fractions import Fraction as Fr
ap = argparse.ArgumentParser()
ap.add_argument('tables')
ap.add_argument('--certificates', default=os.path.dirname(os.path.abspath(__file__)))
ap.add_argument('--mod', type=int, default=0)
ap.add_argument('--json', default='')
args = ap.parse_args()
sys.path.insert(0, args.certificates)
from direct_hilbert import hilbert
P = args.mod
src = open(args.tables).read()
src = re.sub(r'(?<!\\)%.*', '', src)          # strip comments

def macros(src):
    out = {}
    for m in re.finditer(r'\\newcommand\{\\(G\w+)\}\[1\]\{', src):
        i = m.end(); depth = 1
        while depth:
            depth += {'{': 1, '}': -1}.get(src[i], 0); i += 1
        out[m.group(1)] = src[m.end():i - 1]
    return out

def graph(body):
    verts = re.findall(r'\\node\[vertex\]\s*\((\w+)\)', body)
    edges = []; straight = set(); dup = []
    for d in re.findall(r'\\draw(.*?);', body, re.S):
        toks = re.findall(r'\((\w+)\)|(--)|(to\s*(?:\[[^\]]*\])?)', d)
        prev = None; conn = None
        for name, dash, to in toks:
            if name:
                if prev is not None and conn:
                    e = (prev, name)
                    # a straight segment drawn twice is one edge; parallel edges are drawn as distinct to[...] curves
                    if conn == '--' and (e in straight or e[::-1] in straight):
                        dup.append(e)
                    else:
                        edges.append(e)
                        if conn == '--': straight.add(e)
                prev, conn = name, None
            else:
                conn = dash or 'to'
    if dup: print('note: straight segment drawn twice, counted once:', dup)
    idx = {v: i + 1 for i, v in enumerate(verts)}
    assert all(a in idx and b in idx for a, b in edges), (verts, edges)
    return len(verts), [(idx[a], idx[b]) for a, b in edges]

G = {k: graph(v) for k, v in macros(src).items()}
F4 = [{1: 1, 2: 1}, {1: 1, 3: 1}, {1: 1, 2: Fr(1, 2), 3: Fr(1, 3)}, {1: 1, 2: 1, 3: 1}]
F5 = [{1: 1}, {1: 1, 3: 1}, {1: 1, 2: 1}, {1: 1, 2: 1, 3: 2}, {1: 1, 2: Fr(1, 2), 3: Fr(1, 6)}]
F9 = [{1: 1, 4: 1}, {1: 1, 2: 1}, {1: 1, 3: 1}, {1: 1, 3: 1, 4: 1}, {1: 1, 2: 1, 4: 1},
      {1: 1, 2: Fr(1, 2), 3: Fr(1, 3)}, {1: 1, 2: 1, 3: 1}, {1: 1, 2: Fr(1, 2), 3: Fr(1, 3), 4: Fr(1, 4)},
      {1: 1, 2: 1, 3: 1, 4: 1}]
tables = re.split(r'\\begin\{(?:sideways)?table\}', src)[1:]
bad = total = 0
records = []
for tno, tab in enumerate(tables, 1):
    for row in re.findall(r'^\s*\\(G\w+)\{[.\d]+\}(.*?)\\\\\[\d+pt\]', tab, re.M | re.S):
        name, rest = row
        cells = [tuple(int(x) for x in re.sub(r'\$\\\\\$|\$|\s|\\mathbf|[{}]', '', c).strip('()').split(','))
                 for c in re.findall(r'\\begin\{tabular\}\{@\{\}c@\{\}\}(.*?)\\end\{tabular\}', rest, re.S)]
        fs = {9: F9, 5: F5, 4: F4}[len(cells)]
        nv, E = G[name]
        assert len({v for e in E for v in e}) == nv, name   # connected drawing uses every vertex
        for j, (f, want) in enumerate(zip(fs, cells), 1):
            t0 = time.perf_counter(); got = tuple(hilbert(E, f, P)); total += 1
            ok = got == want; bad += not ok
            records.append(dict(table=f'B{tno}', graph=name, vertices=nv,
                                edges=[list(e) for e in E], column=j,
                                f={str(k): str(Fr(a)) for k, a in sorted(f.items())},
                                hilbert_function=list(want), recomputed=list(got), ok=ok))
            print(f'table B{tno} {name:5s} |V|={nv} |E|={len(E):2d} col {j}: {"ok " if ok else "MISMATCH"} {got}' + ('' if ok else f' table says {want}') + f'  ({time.perf_counter()-t0:.1f}s)', flush=True)
print(f'{total} entries checked, {bad} mismatches', '(exact rational arithmetic)' if not P else f'(mod {P})')
if args.json:
    import csv, json
    with open(args.json, 'w') as fh:
        json.dump(dict(source=os.path.basename(args.tables),
                       arithmetic='exact over Q' if not P else f'mod {P}',
                       note='vertices are numbered 1..n in the order of the TikZ nodes; '
                            'f maps degree k to the coefficient of u^k',
                       entries=records), fh, indent=1)
    with open(os.path.splitext(args.json)[0] + '.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['table', 'graph', 'vertices', 'edges', 'column', 'f', 'hilbert_function'])
        for r in records:
            fstr = ' + '.join((f'{c}*' if c != '1' else '') + ('u' if k == '1' else f'u^{k}')
                              for k, c in r['f'].items())
            w.writerow([r['table'], r['graph'], r['vertices'],
                        ' '.join(f'{a}-{b}' for a, b in r['edges']), r['column'], fstr,
                        ','.join(map(str, r['hilbert_function']))])
sys.exit(1 if bad else 0)
