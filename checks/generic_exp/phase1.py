"""Generic and exponential Hilbert series for all graphs of graphs_list (mod a large prime).

Generic: f = u + a_2 u^2 + ... + a_D u^D with random a_k mod p, D = max degree; three independent
draws must agree (each value is a lower bound for the generic cumulative dimensions, attained
outside a proper Zariski closed set).  Exponential: f = e^u - 1 truncated at u^D.
Output: results/phase1.txt and results/phase1.json."""
import json, random, sys, time
from fractions import Fraction as Fr
from math import factorial
from gz import GraphAlg, number_of_forests, PBIG
from graphs_list import all_graphs

random.seed(20261020)
G = all_graphs()
out = {}
lines = []
for name, (n, E) in G.items():
    t = time.time()
    A = GraphAlg(n, E)
    D = A.Delta
    assert A.D == number_of_forests(n, E), name
    gens = []
    for trial in range(3):
        f = {1: 1}
        for k in range(2, D + 1):
            f[k] = random.randrange(1, PBIG)
        gens.append(A.hilbert(f))
    assert gens[0] == gens[1] == gens[2], (name, gens)
    hexp = A.hilbert({k: Fr(1, factorial(k)) for k in range(1, D + 1)})
    hu = A.hilbert({1: 1})
    out[name] = dict(n=n, edges=E, Delta=D, forests=A.D, generic=gens[0], exp=hexp, u=hu)
    flag = '' if hexp == gens[0] else '   <-- exp NON-GENERIC'
    line = f'{name:45s} |V|={n} |E|={len(E):2d} Delta={D}  generic {gens[0]}  exp {hexp}{flag}'
    lines.append(line)
    print(line, f'({time.time() - t:.1f}s)', flush=True)
open('results/phase1.txt', 'w').write('\n'.join(lines) + '\n')
json.dump(out, open('results/phase1.json', 'w'), indent=1)
