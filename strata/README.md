# Common parameter partitions for trees of maximum degree 4 and 5

Computed with the paper's own local algorithm (`root_states` in `certificates/boundary_local.py`), in exact arithmetic. Theorem 6.10 (maximum degree 4) rests on the computations for D = 4 below; the results for D = 5 (Remark 6.11) are data, and the recomputation on their intersections is not done. The collision statistics at the end are those of Tables 9 and 10 and Section 6.9.

## Method

Fix D and normalize α = 1, so that f = u + u² + b u³ + c u⁴ + d u⁵ (as far as D requires).

1. Start from the leaf state, which has gap 1. Repeatedly apply the local vertex rule to every multiset of at most D−1 incoming gaps, until the set of gaps that occur is closed.
2. Do this by exact elimination over ℚ(b, c, …), recording every factor inverted in a pivot.
3. Recompute the whole closure on each recorded hypersurface, over its own function field. Compare its tables with the generic ones and record the new pivot factors.
4. Treat the α = 0 cases after scaling: β = 1, or β = 0 with γ = 1, and so on.
5. As a final check, search all trees of maximum degree ≤ D up to a given size for a tree whose Hilbert function separates two cells (`tree_messages`-style recursion, exact rationals).

If every local table agrees on a set of parameters, every tree of maximum degree ≤ D has the same Hilbert function there. A separating tree shows that two cells really differ.

Notation: g = f⁻¹, and [u^k]g is its k-th coefficient. For α = 1:

- [u³]g = 2 − b
- [u⁴]g = 5b − c − 5
- [u⁵]g = 3b² − 21b + 6c − d + 14

## D = 3 (consistency check)

The closure has gaps 0–3 and 14 local tables. The only pivot factor is b − 2, that is β = 2α². This agrees with Theorem 6.7.

## D = 4: 8 cells

The closure has gaps 0–4 and 55 tables (8 s, `d4_closure.log`). Its pivot factors are:

- b = 2, i.e. [u³]g = 0;
- c = 5b − 5, i.e. [u⁴]g = 0;
- c = b² + b − 1.

On the curve c = b² + b − 1, which is g₂g₄ = g₃², no table changes, so this factor is spurious.

The recomputation on each of the three curves inverts a single new factor: c − 5 on b = 2, and b − 2 on the other two curves. Each vanishes on its curve only at the common point (2, 5), where [u³]g = [u⁴]g = 0, and the computation at (2, 5) inverts no factor. So every point of the three curves is covered. On the α = 0 family (β = 1, γ free) the only factor is γ, and γ = 0 is the single scaling orbit of u + u³. Log: `d4_cells.log`.

The cells are the 8 patterns of vanishing of [u²]g, [u³]g, [u⁴]g:

| Cell | Representative | Vanishing among [u²]g, [u³]g, [u⁴]g |
|---|---|---|
| 1 | u | all three |
| 2 | u + u⁴ | [u²], [u³] |
| 3 | u + u³ | [u²], [u⁴] |
| 4 | u + u³ + u⁴ | [u²] |
| 5 | generic α ≠ 0, e.g. (b, c) = (3, 7) | none |
| 6 | b = 2, c ≠ 5 | [u³] |
| 7 | c = 5b − 5, b ≠ 2 | [u⁴] |
| 8 | (b, c) = (2, 5), i.e. f⁻¹ = u − u² + O(u⁵) | [u³], [u⁴] |

Every pair of cells is separated by a tree with at most 10 vertices (`d4_witness.log`). Notable witnesses:

- **[u³]g = 0:** E₆ = T(1,2,2). The series 1+6t+16t²+9t³ becomes 1+6t+15t²+10t³.
- **[u⁴]g = 0:** the 10-vertex tree with edges 01 02 03 04 15 16 17 28 39. It has two adjacent degree-4 vertices, and 1+10t+47t²+136t³+234t⁴+84t⁵ becomes 1+10t+47t²+135t³+235t⁴+84t⁵.
- **u + u³ + u⁴ against u + u³:** a 9-vertex tree.

For D ≤ 4 this suggests: *the Hilbert function of a tree depends only on which coefficients [u^k]f⁻¹, 2 ≤ k ≤ D, vanish.* D = 5 shows this is not the whole story (below).

## D = 5: candidate loci (partial)

The generic closure has gaps 0–5 and 209 tables (216 s, `d5_closure.log`), with 9 pivot factors. Each was recomputed over its own function field (`surf_*.log`). Five of them change local tables. Each of the five is seen by a tree (`d5_witness.log`, all 4588 trees with Δ ≤ 5 and at most 14 vertices):

| Locus | Witness | Change |
|---|---|---|
| [u³]g = 0 | 6 vertices | drop at t² |
| [u⁴]g = 0 | 10 vertices | drop at t³ |
| [u⁵]g = 0 | 14 vertices, edges 0–1…0–5, 1–6…1–9, 2–10, 3–11, 4–12, 6–13: two adjacent degree-5 vertices with arms of length 2 | drop at t⁴ |
| **Star locus** α²δ − 3αβγ + 2β³ = 0 | the star K₁,₅ | 1+6t+16t²+9t³ becomes 1+6t+15t²+10t³ |
| **S8** | the spider with a degree-5 centre and three arms of length 2 (9 vertices) | 1+9t+37t²+93t³+116t⁴ becomes 1+9t+37t²+92t³+117t⁴ |

On the star locus:

- This is exactly the determinant expressing that f, xf, f² are linearly dependent in K[x]/(x⁶) modulo 1, x, x².
- It is invariant under f ↔ f⁻¹, as the star inversion symmetry of Section 5.3 requires.
- It contains the inverse quadratic and u/(1−u).

On S8: after inversion (normalizing −g(−u)), it becomes the Hankel condition g₃g₅ = g₄².

The other four factors change no table and are spurious: b = 3, c = b² + b − 1, and two cubic surfaces.

**Not yet done for D = 5:** recomputation on all pairwise and triple intersections, and the α = 0 part. Each step takes minutes.

## Running

From this folder, with the path to `certificates/` as the first argument:

```
python3 d4_closure.py ../certificates 4     # or 3
python3 d4_cells.py ../certificates
python3 d4_witness.py ../certificates 10
python3 d5_closure.py ../certificates
python3 d5_surfaces.py ../certificates generic   # writes d5_generic_tables.pkl
python3 d5_surfaces.py ../certificates 0         # ... up to 8, one hypersurface each
python3 d5_witness.py ../certificates 14
```

Requirements: Python 3 and SymPy.

## Collision search among trees (complete test via the partitions)

Theorem 6.7 makes comparison at the four points u, u+u³, u+u²+2u³, u+u² a *complete* test of universal Hilbert equivalence for subcubic trees. Every collision found this way is therefore proved, not sampled. For Δ ≤ 4 the eight-point test is complete by Theorem 6.10.

The tree enumeration reproduces the known counts:

- subcubic trees: 1, 1, 2, 2, 4, 6, 11, 18, 37, … , 11020;
- trees with Δ ≤ 4: the alkane numbers 1, 1, 2, 3, 5, 9, 18, 35, 75, … , 1858.

Logs: `collisions_D3.log` (n ≤ 18), `collisions_D4.log` (n ≤ 14).

**Subcubic trees.**

| n | trees | profile classes | trees in collisions | largest class |
|---|---|---|---|---|
| ≤ 9 | | all distinct | 0 | |
| 10 | 37 | 35 | 4 | 2 |
| 11 | 66 | 65 | 2 | 2 |
| 12 | 135 | 133 | 4 | 2 |
| 13 | 265 | 259 | 12 | 2 |
| 14 | 552 | 542 | 20 | 2 |
| 15 | 1132 | 1108 | 48 | 2 |
| 16 | 2410 | 2368 | 83 | 3 |
| 17 | 5098 | 5005 | 185 | 3 |
| 18 | 11020 | 10850 | 339 | 3 |

At n = 10 the two classes are the T2 pair and a second pair. The second pair is the T2 move (swapping a leaf and a pendant path of length 3 between adjacent trivalent vertices) applied to the tree with edges 01 02 03 14 15 26 27 38 89.

**Δ ≤ 4.** At n = 10 there are 3 classes (6 trees): the two subcubic pairs and one more. At n = 14 there are 30 trees in collisions among 1858.

**How much each single parameter separates** (`info_D3.log`, `info_D4.log`). The full profile is far stronger than any single f.

| | trees | full profile | u | u+u³ | u+u²+2u³ | u+u² (generic) |
|---|---|---|---|---|---|---|
| Subcubic, n = 18 | 11020 | 10850 | 1 | 805 | **3202** | 2298 |

For Δ ≤ 4 at n = 14 (1858 trees, 1843 profile classes), the most informative single point is again an exceptional one: [u³]f⁻¹ = 0 gives 634 classes, against 245 for generic f.

**Is there a simple move behind the collisions?** (`moves_D3.log`, `swaps_D3.log`)

- Up to n = 15, every subcubic collision pair differs by swapping two pendant paths of different lengths. At n = 16, four pairs need more than one swap.
- No swap type preserves the profile in general. For example, the T2 move (leaf against a length-3 path, across an edge) preserves it in 16 of the 250 cases with n ≤ 14.
- So the collisions look sporadic, not the result of a local operation.

## Usage of the collision scripts

`python3 collisions.py ../certificates D NMAX`, and the same arguments for
`collision_moves.py`, `info_per_point.py` and `swap_test.py`. The logs were
produced with D = 3, NMAX = 18 and D = 4, NMAX = 14 (`collisions_*`,
`info_*`), D = 3, NMAX = 16 (`moves_D3.log`) and D = 3, NMAX = 14
(`swaps_D3.log`). `../run_all.sh strata` reruns everything in this
directory.
