#!/usr/bin/env python3
"""Exact boundary-state recursion for finite trees.

This implements messages across edges; it does not build the full 2^|E|
edge algebra. Coefficients are exact fractions. The local operation is the
rank calculation in boundary_local.root_states.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import product
from typing import Iterable, Mapping
from boundary_local import root_states, close


def message_function(edges: Iterable[tuple[int,int]], coefficients: Mapping[int, object]):
    edges=list(edges)
    adjacency=defaultdict(list)
    for u,v in edges:
        if u==v:raise ValueError('A tree has no loops.')
        adjacency[u].append(v);adjacency[v].append(u)
    if len(edges)!=len(adjacency)-1:raise ValueError('The input must be a nonempty tree.')
    seen=set();stack=[next(iter(adjacency))]
    while stack:
        v=stack.pop()
        if v in seen:continue
        seen.add(v);stack.extend(adjacency[v])
    if len(seen)!=len(adjacency):raise ValueError('The input must be connected.')
    coeff={d:Fraction(a) for d,a in coefficients.items()}
    if coeff.get(1)!=1:raise ValueError('Normalize the linear coefficient to 1.')
    @lru_cache(None)
    def local(deltas):return root_states(deltas,coeff)
    @lru_cache(None)
    def message(v,parent):
        child=[w for w in adjacency[v] if w!=parent]
        spaces=[message(w,v) for w in child]
        answer=Counter()
        for terms in product(*(tuple(S.items()) for S in spaces)):
            base=sum(pair[0] for pair,multiplicity in terms)
            deltas=tuple(sorted(pair[1]-pair[0] for pair,multiplicity in terms))
            weight=1
            for pair,multiplicity in terms:weight*=multiplicity
            for (p,q),multiplicity in local(deltas).items():
                answer[p+base,q+base]+=weight*multiplicity
        assert all(p>=0 and q>=0 for p,q in answer)
        return answer
    return adjacency,message


def hilbert_tree(edges,coefficients):
    edges=list(edges)
    if not edges:return [1]
    adjacency,message=message_function(edges,coefficients)
    u,v=edges[0]
    return close(message(u,v),message(v,u))


def boundary_states_at_leaf(edges,leaf,coefficients):
    adjacency,message=message_function(edges,coefficients)
    if len(adjacency[leaf])!=1:raise ValueError('The marked vertex must be a leaf.')
    return message(adjacency[leaf][0],leaf)

if __name__=='__main__':
    from direct_hilbert import GRAPHS
    for name in ['T2_1','T2_2']:
        print(name,hilbert_tree(GRAPHS[name]+[(6,11)],{1:1,2:1,3:2}))
