# Exact Hilbert-stratum certificates

Statement numbers refer to the first arXiv version of the paper;
`../README.md` lists the script for each statement. This directory
contains the certificates for Proposition 6.4, Lemma 6.6, Theorem 6.7 and
the pairs T1, T2 and G8 of Proposition 6.15, and the example of Remark 6.19
(`tree_messages.py`). The pair T3 is checked in
`../checks/t3_profile.py`, and the partition for maximum degree four
(Theorem 6.10) is computed in `../strata/`.

These calculations use the definition **Y_v = f(X_v)** in the manuscript,
not the presentation obtained by substituting f where f^{-1} belongs.
The normalization is f(u) = u + a u^2 + b u^3 + c u^4.
All mathematical arithmetic in the certificates is exact. Timings, when
printed, are the only floating-point quantities.

## Run

From this directory:

```sh
python3 verify_parameter_results.py
```

Requirements: Python 3 with SymPy; a C++17 compiler and GMP development
headers/libraries for the unicycle calculations. The Python script builds
the C++ accelerator itself. `--skip-unicycles` runs the tree and small-graph
certificates without the compiler. Recorded complete output is included in
`verified_results.txt` and `verified_output/`.

`tree_messages.py` also supplies an exact boundary-state algorithm for
individual trees, without constructing their full edge algebras.

Files: `verify_parameter_results.py` is the driver. It uses
`boundary_local.py` (local boundary states over a field of rational
functions), `direct_hilbert.py` (the filtration computed directly in the
edge algebra, with the edge lists of G8, T1 and T2),
`direct_symbolic_degree.py` (the same over Q(b), recording every inverted
pivot) and `certify_unicycle.cpp`. `Hilbert_profiles_exact.json`
summarizes the complete Hilbert profiles of G8, T1 and T2 that these
certificates establish (Tables 5–7 of the paper).

## What makes the calculations exhaustive in the parameters?

### Boundary calculation

For a vertex with r incoming edges x_1,...,x_r and outgoing edge z, use

    E = k[x_1,...,x_r,z]/(x_1^2,...,x_r^2,z^2),
    y = f(z - x_1 - ... - x_r).

An incoming state [p_j,q_j] contributes either 1 at cost p_j or x_j at
cost q_j. After subtracting the common shift sum p_j, put d_j=q_j-p_j.
The local filtration is spanned by x_S y^k with cost sum_{j in S} d_j+k.
Only 0 <= k <= r+1 need be considered. Negative normalized costs are handled
by a common shift; actual message birth degrees remain nonnegative.

Let P_i be that local filtration, Q_i its image modulo z, and z R_i its
intersection with (z). The output states are recovered from

    dim(Q_i intersection R_j)
      = rank(z P_i) + rank(P_j) - rank(P_j + z P_i)

by taking two successive finite differences. `boundary_local.py` computes
these spaces directly by exact Gaussian elimination.

The triangular changes of basis underlying the leaf-splitting lemma are
linear over the algebra in the local edge variables. They therefore commute
with multiplication by y, so independent input states may be treated
separately and the resulting state multiplicities summed.

Two sides of an edge close by the bilinear rule

    E([p,q],[p',q']) = t^(p+p') + t^min(q+p',p+q').

The sign reversal required for the shared edge on one side does not alter
this rule or either birth degree.

### Subcubic trees: all parameters, all tree sizes

For a != 0 set a=1 by scaling and write b for b/a^2. Starting with a leaf,
all boundary states have gaps in {0,1,2,3}. The verifier computes the ten
symmetric binary-vertex transitions over Q(b). Every pivot factor divides
b-2, and all ten transitions are separately computed at b=2. The output
gaps remain in {0,1,2,3}. The bivalent path-extension rule preserves this set as well.
Induction over any rooted subcubic tree proves that all its Hilbert series
are constant on the two cells a!=0,b!=2a^2 and a!=0,b=2a^2.

On a=0, scaling after a characteristic-zero field extension gives two
possibilities: b=0 or b=1. Such extensions do not change ranks. Therefore
four test functions suffice for every subcubic tree:

    u, u+u^3, u+u^2, u+u^2+2u^3.

### The ten-vertex tree pairs T1 and T2

The first pair, T1, has a degree-four and a degree-three vertex joined by a
path of length two. The arm lengths off this path are respectively

    (1,1,1) | (2,2),     (2,1,1) | (2,1).

The second pair, T2, has two degree-three vertices joined by a single edge,
with arm lengths

    (1,1) | (3,3),       (1,3) | (1,3).

For a=1 the only local matrices have dimensions 8 and 16. For the first
pair, the generic symbolic computation over Q(b,c) excludes only

    b-2,  c-5b+5,  c-b^2-b+1.

Each of the three curves is recomputed over a one-variable rational
function field. These recomputations invert only c-5 on b=2 and only b-2
on the other two curves, and these factors vanish on the curves only at
(b,c)=(2,5), where all three curves meet. That point is recomputed over
Q, and this computation inverts no factor. There are no other excluded
parameters. The two apparent extra curves have the generic state data
away from (2,5).

For a=0 the exhaustive cases are b=1 with c arbitrary, b=0,c=1, and the
origin. Their small boundary calculations have no exceptional factors.
Bivalent vertices on terminal arms contribute factors 1+t on this locus.
Thus the displayed six-cell classification for the first tree pair and
four-cell classification for the second are exhaustive, not sampled.

### The eight-vertex unicycles

For a=1, `certify_unicycle.cpp` computes the filtration in the full
square-zero edge algebra over Q(b). Its Gaussian pivots certify the ranks
on the complement of the following seven irreducible factors:

    b-2, b-2/3, b-4/3, b-8/3, b-10/3,
    b^2-(10/3)b+8/9, b^2-(8/3)b-4/9.

The code then repeats the rank computation over Q[b]/(p) for EACH factor
p, including the two quadratic number fields. Thus the irrational and
non-rational parameter values are checked exactly as well. All factors
except b-2 and b-2/3 are artifacts of the chosen pivots, not genuine rank
drops. The locus a=0 has the two scaling cases b=0 and b!=0, each checked
over Q. This proves equality for every normalized nondegenerate function,
over any characteristic-zero field.

A basis of F_i is built from a basis of F_(i-1) and products of its newly
added basis vectors with the Y_v. Products of older basis vectors already
belong to F_(i-1), so the frontier optimization is exact. The known
unfiltered dimension is 240 for both unicycles; attaining it certifies
exhaustion. Generic pivots give a finite exceptional set because every
basis operation is invertible after localizing at those pivots. At each
exceptional irreducible factor the computation starts again from the
original generators, not from the possibly singular generic basis.

## General bounded-degree consequence (proof, not an experimental claim)

Let a graph of maximum non-loop degree D have a marked leaf edge r, with
other endpoint of degree d. For its open boundary filtrations Q_i and R_i,

    Q_i is contained in R_(i+D),
    R_i is contained in Q_(i+d-2).

Indeed r is, up to sign, the sum of f^{-1}(Y_v) over all vertices except
the marked leaf, so r belongs to P_D. For the other inclusion, coefficient
extraction in r maps P_i into f'(X_bar) Q_(i-1), where X_bar is the generator
at the neighbor after deleting r. Since X_bar^d=0, f'(X_bar) is a polynomial
of degree at most d-1 in f(X_bar). Therefore every boundary state satisfies

    2-d <= q-p <= D,

and in particular belongs to one of the 2D-1 gap types from 2-D through D.
Shifts of both birth degrees are recorded by powers of the Hilbert variable.

Consequently all trees of maximum degree at most D can be evaluated using
finitely many local transition tables, with local edge-algebra dimension
at most 2^D, independent of tree size. The coefficient matrices for those
finitely many prototypes have polynomial entries in the deformation
parameters. Their simultaneous rank strata give a finite locally closed
parameter partition on which the Hilbert series of EVERY such tree is
constant. This is a common partition for the whole class, not a separate
existence statement for each individual tree.

The explicit complete common partition is computed here for D=3
(Theorem 6.7). The partition for D=4 (Theorem 6.10) is computed in
`../strata/`; here only the local prototypes required for the pair T1 are
classified.

## Small graphs

The verifier also recomputes the complete parameter dependence for K4 and
for the triangle with a doubled edge (Proposition 6.4). For a=1, the
exceptional finite b-values are:

    K4: b=2;
    triangle with a doubled edge: b=2 and b=2/3.

Here b is the coefficient of f itself, not of the series f^(-1) that
appears in the defining relations; the two are related by

    f^(-1)(u) = u-a u^2+(2a^2-b)u^3+(-5a^3+5ab-c)u^4+...

## Optional attachment check (Remark 6.19)

```sh
python3 tree_messages.py
```

At f=u+u^2+2u^3, attach one edge at vertex 6 of each tree of the pair
T2. The resulting Hilbert functions are

    1,11,59,192,392,369
    1,11,59,192,386,375.

Their difference is 6t^4(1-t). Thus pointwise equality of the original
Hilbert profiles does not imply equality of marked boundary profiles or
preservation under attaching a tail.

## Scope

None of the certificates concerns degenerate series f with f'(0)=0
(see `../search/degenerate_check.py` for those). All higher
coefficients beyond the maximum vertex degree are irrelevant by
nilpotence, and arbitrary nonzero linear and constant terms reduce to the
normalization above.
