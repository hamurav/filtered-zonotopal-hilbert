"""The graphs of this probe: all connected simple graphs with 2-5 vertices, the multigraphs of the
tables in ../../tables/ (formerly Appendix B), the two multigraphs of [NSh, Example 2], the doubled triangle and the
two-vertex multigraphs E_2, E_3 (theta), E_4.  Disconnected graphs reduce to their components."""
import itertools
from gz import connected_graphs, canon

# names of the simple graphs (macro of tables345_clean.tex in brackets)
NAMED = {
    'K2': (2, [(0, 1)]),
    'P3': (3, [(0, 1), (1, 2)]),
    'K3': (3, [(0, 1), (1, 2), (0, 2)]),
    'P4': (4, [(0, 1), (1, 2), (2, 3)]),
    'C4': (4, [(0, 1), (1, 2), (2, 3), (3, 0)]),
    'K13 (star, GDfo)': (4, [(0, 1), (0, 2), (0, 3)]),
    'paw (GPaw)': (4, [(0, 1), (0, 2), (0, 3), (2, 3)]),
    'diamond (GDia)': (4, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)]),
    'K4 (GKfo)': (4, list(itertools.combinations(range(4), 2))),
    'P5': (5, [(0, 1), (1, 2), (2, 3), (3, 4)]),
    'C5': (5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
    'fork (GDfi)': (5, [(0, 1), (1, 2), (2, 3), (2, 4)]),
    'K14 (star, Gsta)': (5, [(0, 1), (0, 2), (0, 3), (0, 4)]),
    'banner (Gban)': (5, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4)]),
    'bull (Gbul)': (5, [(0, 1), (1, 2), (2, 0), (1, 3), (2, 4)]),
    'lollipop L32 (Glol)': (5, [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4)]),
    'cricket (Gcri)': (5, [(0, 1), (1, 2), (2, 0), (0, 3), (0, 4)]),
    'K23 (Gbip)': (5, [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4)]),
    'kite (Gkit)': (5, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 4)]),
    'house (Ghou)': (5, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4)]),
    'butterfly (Gbut)': (5, [(0, 1), (0, 2), (1, 2), (0, 3), (0, 4), (3, 4)]),
    'dart (Gdar)': (5, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (0, 4)]),
    'K4 subdivided (Gsub)': (5, [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 4), (4, 3)]),
    'gem (Ggem)': (5, [(0, 1), (1, 2), (2, 3), (4, 0), (4, 1), (4, 2), (4, 3)]),
    'K113 (GKKK)': (5, [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4)]),
    'lollipop L41 (Gllo)': (5, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (3, 4)]),
    'wheel W5 (Gwhe)': (5, [(0, 1), (1, 2), (2, 3), (3, 0), (4, 0), (4, 1), (4, 2), (4, 3)]),
    'house X (Ghox)': (5, [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3), (0, 4), (1, 4)]),
    'K5-e (GKme)': (5, [e for e in itertools.combinations(range(5), 2) if e != (3, 4)]),
    'K5 (GKfi)': (5, list(itertools.combinations(range(5), 2))),
}

MULTI = {
    'P4 first edge doubled (Gafa)': (4, [(0, 1), (0, 1), (1, 2), (2, 3)]),
    'P4 middle edge doubled (Gafb)': (4, [(0, 1), (1, 2), (1, 2), (2, 3)]),
    'P4 outer edges doubled (Gafc)': (4, [(0, 1), (0, 1), (1, 2), (2, 3), (2, 3)]),
    'paw, triangle edge doubled (Gppw)': (4, [(0, 1), (0, 2), (0, 3), (2, 3), (2, 3)]),
    'C4 one edge doubled (Gcfa)': (4, [(0, 1), (0, 1), (1, 2), (2, 3), (3, 0)]),
    'C4 opposite edges doubled (Gcfb, NSh Ex.2)': (4, [(0, 1), (0, 1), (1, 2), (2, 3), (2, 3), (3, 0)]),
    'C4 adjacent edges doubled (NSh Ex.2)': (4, [(0, 1), (0, 1), (1, 2), (1, 2), (2, 3), (3, 0)]),
    'triangle one edge doubled': (3, [(0, 1), (0, 1), (1, 2), (2, 0)]),
    'E2 (double edge)': (2, [(0, 1)] * 2),
    'E3 (theta)': (2, [(0, 1)] * 3),
    'E4': (2, [(0, 1)] * 4),
}


def all_graphs():
    out = dict(NAMED)
    # check that NAMED is exactly the list of connected simple graphs on 2..5 vertices
    keys = {canon(n, E) for n, E in NAMED.values()}
    assert len(keys) == len(NAMED)
    count = 0
    for nv in range(2, 6):
        for E in connected_graphs(nv):
            count += 1
            assert canon(nv, E) in keys, (nv, E)
    assert count == len(NAMED) == 30, count
    out.update(MULTI)
    return out


if __name__ == '__main__':
    G = all_graphs()
    print(len(G), 'graphs; the 30 connected simple graphs on 2..5 vertices are all named')
