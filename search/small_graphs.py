"""Nonisomorphic trees and connected unicyclic simple graphs on few vertices.

Standard library only.  Graphs are returned as edge lists on the vertices
0, ..., n-1.  Isomorphism classes are separated by complete canonical forms:

* a tree is encoded by the AHU string of the tree rooted at its center
  (the smaller of the two strings if there are two centers);
* a unicyclic graph is encoded by the cyclic sequence of AHU strings of the
  rooted trees hanging from the cycle vertices, minimized over the 2k
  rotations and reflections of the cycle.

Both encodings determine the graph up to isomorphism, and isomorphic graphs
have equal encodings.  New graphs on n vertices are obtained by attaching a
leaf to every vertex of every graph on n-1 vertices (every tree, and every
unicyclic graph other than a cycle, has a leaf); cycles are added directly.

Running this file checks the counts against OEIS A000055 (trees) and
A001429 (connected unicyclic graphs).
"""
from collections import defaultdict


def adjacency(n, edges):
    adj = defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    for v in range(n):
        adj[v]
    return adj


def _rooted(adj, v, parent, blocked=()):
    """AHU string of the subtree at v away from parent, avoiding blocked vertices."""
    kids = sorted(_rooted(adj, w, v, blocked) for w in adj[v] if w != parent and w not in blocked)
    return '(' + ''.join(kids) + ')'


def tree_key(n, edges):
    if n == 1:
        return '()'
    adj = adjacency(n, edges)
    deg = {v: len(adj[v]) for v in adj}
    layer = [v for v in adj if deg[v] == 1]
    left = n
    while left > 2:
        nxt = []
        for v in layer:
            left -= 1
            for w in adj[v]:
                deg[w] -= 1
                if deg[w] == 1:
                    nxt.append(w)
        layer = nxt
    return min(_rooted(adj, c, None) for c in layer)


def cycle_of(n, edges):
    """Vertices of the unique cycle, in cyclic order."""
    adj = adjacency(n, edges)
    deg = {v: len(adj[v]) for v in adj}
    alive = set(adj)
    stack = [v for v in adj if deg[v] == 1]
    while stack:
        v = stack.pop()
        alive.discard(v)
        for w in adj[v]:
            if w in alive:
                deg[w] -= 1
                if deg[w] == 1:
                    stack.append(w)
    start = min(alive)
    order, prev, cur = [start], None, start
    while True:
        nxt = [w for w in adj[cur] if w in alive and w != prev]
        if prev is None:
            nxt = nxt[:1]
        w = nxt[0]
        if w == start:
            break
        order.append(w)
        prev, cur = cur, w
    return adj, order


def unicycle_key(n, edges):
    adj, cyc = cycle_of(n, edges)
    on_cycle = set(cyc)
    labels = [_rooted(adj, c, None, on_cycle - {c}) for c in cyc]
    k = len(labels)
    best = None
    for seq in (labels, labels[::-1]):
        for s in range(k):
            cand = tuple(seq[s:] + seq[:s])
            if best is None or cand < best:
                best = cand
    return best


def trees(nmax):
    """Yield (n, edges) for every tree with 1 <= n <= nmax vertices, one per class."""
    level = {tree_key(1, []): []}
    yield 1, []
    for n in range(2, nmax + 1):
        nxt = {}
        for edges in level.values():
            for v in range(n - 1):
                new = edges + [(v, n - 1)]
                nxt.setdefault(tree_key(n, new), new)
        level = nxt
        for key in sorted(level):
            yield n, level[key]


def unicycles(nmax):
    """Yield (n, edges) for every connected unicyclic simple graph with 3 <= n <= nmax vertices."""
    level = {}
    for n in range(3, nmax + 1):
        nxt = {}
        cycle = [(i, (i + 1) % n) for i in range(n)]
        nxt[unicycle_key(n, cycle)] = cycle
        for edges in level.values():
            for v in range(n - 1):
                new = edges + [(v, n - 1)]
                nxt.setdefault(unicycle_key(n, new), new)
        level = nxt
        for key in sorted(level):
            yield n, level[key]


A000055 = {1: 1, 2: 1, 3: 1, 4: 2, 5: 3, 6: 6, 7: 11, 8: 23, 9: 47, 10: 106, 11: 235, 12: 551, 13: 1301}
A001429 = {3: 1, 4: 2, 5: 5, 6: 13, 7: 33, 8: 89, 9: 240, 10: 657, 11: 1806}

if __name__ == '__main__':
    from collections import Counter
    ct = Counter(n for n, _ in trees(12))
    cu = Counter(n for n, _ in unicycles(10))
    for n in sorted(ct):
        assert ct[n] == A000055[n], (n, ct[n])
    for n in sorted(cu):
        assert cu[n] == A001429[n], (n, cu[n])
    print('trees     n=1..12:', [ct[n] for n in sorted(ct)], 'matches OEIS A000055')
    print('unicycles n=3..10:', [cu[n] for n in sorted(cu)], 'matches OEIS A001429')
