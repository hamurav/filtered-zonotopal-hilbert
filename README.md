# Code and data for "Hilbert functions of filtered graphical zonotopal algebras"

Boris Shapiro, Ilya Smirnov, Arkady Vaintrob.

This repository contains the programs behind every computer-assisted
statement of the paper, with their recorded output, tables of Hilbert
functions of small graphs (`tables/`, formerly Appendix B of the paper),
and the data for its conjectures. Statement numbers refer to the first
arXiv version of the paper. The statements are also named, so that the
table below can be matched with other versions.

## What kind of computation each statement rests on

Every computation in this repository is one of three kinds.

- **Exact.** Integer, rational, rational-function-field or number-field
  arithmetic. These computations are used in proofs.
- **Certified modular.** A rank computed modulo a prime p equals a rank over
  Q that is proved, in the paper, to be at most the same number. Since a rank
  modulo p is at most the rank over Q, equality over Q follows; see (6.21).
  This is used for Proposition 6.26 (Conjecture 6.25 for graphs with at most
  seven vertices) and for the annihilators of Section 7.
- **Modular evidence.** Ranks modulo a prime, which are lower bounds only.
  This is used for the cycle series of Conjecture 7.8 (n ≤ 17), the scans
  behind the proposed loci of Remark 6.21, the six-vertex data of
  Remark 6.23, the unimodality count of the Outlook, and one consistency
  check of Example 4.17.

## Layout

| Directory | Contents | Requirements |
|---|---|---|
| `certificates/` | parameter certificates: the subcubic branching table, small graphs, the collision pairs T1, T2 and G8 | Python 3, SymPy; C++17 and GMP for `certify_unicycle.cpp` |
| `strata/` | the parameter partition for maximum degree 4, the loci for degree 5, collision statistics | Python 3, SymPy |
| `search/` | exhaustive searches over small trees and unicyclic graphs; exact checks of the four collision pairs, of degenerate functions and of the exponential on five vertices | Python 3.9 or later, standard library only |
| `checks/` | the pair T3 by boundary states, paths and tails, necks, exponential examples, D6/E6, C6, the loci identities of Remarks 6.21–6.22, subdivision families | Python 3, SymPy for `loci_algebra.py` |
| `checks/generic_exp/` | generic and exponential series of small graphs modulo a prime; the certificate of Proposition 6.26 | Python 3, NumPy, networkx |
| `cycles/` | the cycle series H_n (n ≤ 17), the cycle formulas, the annihilators of Section 7 | Python 3, NumPy |
| `tables/` | Tables B1–B4 (formerly Appendix B, Tables 12–15) as PDF, LaTeX, JSON and CSV, with the scripts that recompute them | Python 3, SymPy for `complete_profiles.py` |

`search/` was written from scratch, independently of the other directories,
and imports nothing from them. `certificates/`, `strata/`, `cycles/`, `tables/`
and `checks/generic_exp/` have their own READMEs with details; the scripts in
`search/` and `checks/` are described below and in their headers.

## Running

```
python3 -m pip install -r requirements.txt
./run_all.sh --quick          # about 25 minutes on four cores
./run_all.sh                  # everything: about 3 hours on four cores
./run_all.sh search strata    # selected directories
```

`run_all.sh` runs every command of the table below from the directory it
names, writes the output to the file listed in the table (replacing the
recorded output) and appends one line per command, with the time taken, to
`run_all.log`. `--quick` skips the long runs: the subdivision families, the
complete cubic profiles, the strata for degree 5, the collision statistics,
and the modular scans, their classification, the six-vertex data and the
summary of Section 6.10. It also replaces the tree search with direct
verification by the search without it, the exact table check by a check
modulo a prime, and the cycles up to n = 17 by n ≤ 13. The recorded output
in this repository is from a complete run; see [Recorded run](#recorded-run).

## Statements and scripts

Run each command from the directory in the second column. "Exact",
"certified" and "evidence" refer to the three kinds above.

| Statement | Directory | Command | Kind | Output |
|---|---|---|---|---|
| Theorem 3.20 (degenerate forests), Lemma 3.23 and the recursion of Remark 3.27; checks only, the statements are proved in the paper | `search/` | `python3 degenerate_check.py trees --nmax 12`; `python3 degenerate_check.py lemma` | exact | `results/degenerate_trees.txt`, `results/degenerate_lemma.txt` |
| Proposition 3.29 (degenerate f on cycles, checked for n ≤ 12) and Remark 3.30 | `search/` | `python3 degenerate_check.py cycles` | exact | `results/degenerate_cycles.txt` |
| Table 1, Corollary 4.14 (paths), Remark 4.13; Theorem 4.9 for two marked graphs (consistency) | `checks/` | `python3 paths_tails.py` | exact | `paths_tails.log` |
| Example 4.17 and Table 2 (gap polynomials A_j of eight marked graphs, neck constants) | `checks/` | `python3 neck_states.py ../certificates` | exact | `neck_states.log` |
| Example 4.17, direct values of L(H_n) for θ∗θ (n ≤ 8), K_1,3∗K_1,3 and K_3∗K_3 (n ≤ 6); a consistency check | `checks/` | `python3 neck_direct.py theta 8` (also `K13 6`, `K3 6`) | evidence (mod 2³¹−1) | `neck_theta8.log`, `neck_K13.log`, `neck_K3.log` |
| Proposition 6.4 (K4, doubled triangle); Lemma 6.6 and Theorem 6.7 (subcubic branching table, four cells); Proposition 6.15 and Tables 5–7 for T1, T2, G8 | `certificates/` | `python3 verify_parameter_results.py` (`--skip-unicycles` without the C++ part) | exact | `verified_results.txt`, `verified_output/` |
| Proposition 6.15 and Table 8 for T3 (boundary states) | `checks/` | `python3 t3_profile.py ../certificates ../strata/d4_witness_lib.py` | exact | `t3_profile.log` |
| Theorem 6.10 and Table 4 (eight cells for maximum degree 4) | `strata/` | `python3 d4_closure.py ../certificates 4`, `python3 d4_cells.py ../certificates`, `python3 d4_witness.py ../certificates 10` | exact | `d4_closure.log`, `d4_cells.log`, `d4_witness.log` |
| Remark 6.11 (loci for maximum degree 5) | `strata/` | `python3 d5_closure.py ../certificates`; `python3 d5_surfaces.py ../certificates generic`, then `0` to `8`; `python3 d5_witness.py ../certificates 14` | exact (intersections not recomputed) | `d5_closure.log`, `surf_*.log`, `d5_witness.log` |
| Remark 6.13 and Tables B1–B3 (complete cubic profiles) | `tables/` | `python3 complete_profiles.py tables345_clean.tex` | exact | `complete_profiles.log` |
| Propositions 6.15 and 6.16, Remark 6.17: the four pairs at ten nondegenerate functions (direct, and boundary states for trees), degenerate functions, bizonotopal series | `search/` | `python3 verify_pairs.py --jobs 4` | exact | `results/pairs.txt` |
| Proposition 6.18 (i): trees with at most 10 vertices | `search/` | `python3 search.py trees --nmax 10 --direct --jobs 4` | exact | `results/trees_n10.{txt,json}` |
| Proposition 6.18 (ii): unicyclic graphs with at most 8 vertices | `search/` | `python3 search.py unicycles --nmax 8 --jobs 4` | exact | `results/unicycles_n8.{txt,json}` |
| Remark 6.19 (attachment example) | `certificates/` | `python3 tree_messages.py` | exact | `tree_messages.log` |
| Section 6.8: the exponential jet on four small graphs; Table 3 (D6, E6) | `checks/` | `python3 exp_check.py ../certificates`; `python3 d6e6_check.py ../certificates` | exact | `exp_check.log`, `d6e6_check.log` |
| [NSh, Example 2], quoted in Section 6.8 | `search/` | `python3 nsh_example2.py` | exact | `results/nsh_example2.txt` |
| Tables 9 and 10 (collision statistics) and the exchange statistics of Section 6.9 | `strata/` | `python3 collisions.py ../certificates 3 18` and `4 14`; `python3 info_per_point.py ../certificates 3 18` and `4 14`; `python3 collision_moves.py ../certificates 3 16`; `python3 swap_test.py ../certificates 3 14` | exact | `collisions_D*.log`, `info_D*.log`, `moves_D3.log`, `swaps_D3.log` |
| Proposition 6.20 and Table 11 (the exponential on six graphs with five vertices); the rational points of the loci of Remark 6.21 | `search/` | `python3 exp_small.py --jobs 4` | exact | `results/exp_small.log` |
| Remark 6.21 (loci from complete scans over GF(101) and GF(61)) | `checks/generic_exp/` | `python3 scan.py 101 ALL3`, `... 101 ALL4`, `... 61 ALL3 ALL4`; `python3 classify.py 101`, `... 61` | evidence | `results/scan*`, `results/classify_*` |
| Remarks 6.21 and 6.22 (the equations of the loci, their incidences, the conic Q, the quartic W) | `checks/` | `python3 loci_algebra.py` | exact (SymPy) | `loci_algebra.log` |
| Remark 6.23 (six vertices) and the unimodality count of the Outlook | `checks/generic_exp/` | `python3 phase1.py`, `python3 six.py 6 15`, `python3 bonus_low.py`, `python3 summary.py` | evidence | `results/six_6*`, `results/summary.txt` |
| Proposition 6.26: Conjecture 6.25 for all graphs with 4 to 7 vertices: part (i) for 968 graphs, part (ii) for 557 graphs, for generic f; the exponential attains the bounds on 968 and 555 of them, and on the other two (atlas G150, G152) h_3 is one below the bound, by an exact computation | `checks/generic_exp/` | `python3 low_certificate.py --jobs 4` | certified; exact for G150, G152 | `results/low_certificate.{txt,json}` |
| Theorems 7.1 and 7.3, Proposition 7.7, Conjecture 7.8, Proposition 7.9, Remark 7.10: the series H_n of C_n for n ≤ 17 and the formulas | `cycles/` | `python3 cycle_hilbert.py 3 17 --primes 2 --out results/cycle_series.json`; `python3 cycle_formulas.py` | H_n: evidence (two primes); formulas: exact | `results/cycle_series.json`, `results/cycle_formulas.log` |
| Lemma 7.4 (r ≤ 5), (7.17), Proposition 7.5 (m ≤ 5) and the six annihilators of C_6 in Section 7.5 | `cycles/` | `python3 annihilators.py` | certified | `results/annihilators.log` |
| the C_6 annihilators, by a second, independent script | `checks/` | `python3 c6_annihilator.py` | exact | `c6_annihilator.log` |
| Conjecture 8.2 (subdivision): Θ(1,1,n) for n ≤ 8 and K_4 with one edge subdivided, n ≤ 5 | `checks/` | `python3 subdivision.py` | exact | `subdivision.log` |
| Tables B1–B4 (formerly Appendix B), all 189 entries | `tables/` | `python3 check_tables.py tables345_clean.tex --certificates ../certificates --json tables.json` | exact | `check_tables_output.txt`, `tables.json`, `tables.csv` |

## The searches of Proposition 6.18

`search/small_graphs.py` generates trees and connected unicyclic simple
graphs up to isomorphism, by adding a leaf to every vertex of every graph
with one vertex fewer (and adding the cycles), with complete canonical
forms (AHU strings at the center of a tree; the cyclic sequence of rooted
AHU strings, minimized over rotations and reflections, for a unicyclic
graph). Running it checks the counts against OEIS A000055 (trees, n ≤ 12)
and A001429 (unicyclic graphs, n ≤ 10).

`search/hilbert_direct.py` computes the filtration directly from the
definition: Y_v = f(X_v) in Q[φ_e]/(φ_e²), F_d spanned by products of at
most d of the Y_v. Spans are kept in fraction-free integer echelon form, so
all ranks are exact ranks over Q. It works for every polynomial f with
f(0) = 0, including degenerate f (then the Y_v generate a proper
subalgebra).

`search/hilbert_boundary.py` implements the boundary-state calculus of
Section 4 for trees: leaf state [0,1], the local vertex rule, and the
closing rule t^(p+p') + t^min(q+p', p+q').

`search/search.py` computes the Hilbert functions of every graph at the
test functions

    u,  u+u^2,  u+u^2+u^3,  u+u^2+2u^3,  u+u^3,  u+u^4

and reports the classes of graphs with the same number of vertices that no
test function separates. Two graphs separated by one nondegenerate f are
not universally Hilbert equivalent, so every universal collision lies in a
reported class. With `--direct`, every tree value is also computed directly
and compared with the boundary-state value; the run stops if they differ.
For each nondegenerate f it also checks that the dimensions add up to the
number of spanning forests.

The function u+u^4 is needed for trees with a vertex of degree at least
four. Without it

    python3 search.py trees --nmax 10 --tag _without_u4 --functions u u+u^2 u+u^2+u^3 u+u^2+2u^3 u+u^3

leaves eight unseparated pairs with 5 to 9 vertices, each consisting of a
tree of maximum degree four and a tree of maximum degree three (the
smallest is the star with four leaves against the tree with degree
sequence (3,2,1,1,1)), and two such pairs with 10 vertices besides T1, T2,
T3. The output is in `search/results/trees_n10_without_u4.txt`.

Results (`search/results/`):

* trees, n ≤ 10: all 201 trees are separated except the pairs T1, T2, T3
  (10 vertices); all boundary-state values agree with the direct values;
* unicyclic graphs, n ≤ 8: all 143 graphs are separated except the pair G8.

These pairs are universally Hilbert equivalent by Proposition 6.15. This
proves Proposition 6.18.

## Degenerate functions (Section 3.7, Proposition 6.16)

For degenerate f (no linear term) the elements f(X_v) generate a proper
subalgebra of C_G. Theorem 3.20 says that for a forest and f of order
k ≥ 2 its Hilbert series is s_k, the number of partial score vectors with
all entries divisible by k, counted by weight/k. The theorem is proved in
the paper. `search/degenerate_check.py` only confirms it:

* `trees --nmax 12`: for all 986 trees with 2 to 12 vertices, the Hilbert
  functions at u^2, u^2+u^3, u^2-2u^3+5u^4, u^3, u^3+u^4, 2u^3-u^5, u^4,
  u^4+3u^5, computed directly in the edge algebra, equal s_k; the
  recursion of Remark 3.27 agrees with the enumeration of score vectors;
* `lemma`: Lemma 3.23 on 470 random families of exponent sets;
* `cycles`: Proposition 3.29 for the cycles C_3, …, C_12 (also proved
  in the paper), and the data of Remark 3.30 (unicyclic graphs, small
  simple graphs and multigraphs), where s_k is only an upper bound.

`search/verify_pairs.py` computes the Hilbert functions of the four
collision pairs for f = u^2, u^2+u^3, u^2+2u^3, u^2+u^3+u^4, u^2+u^4,
u^3, u^3+u^4, u^4. The T3 rows agree with s_2 and s_3 as used in the
proof of Proposition 6.16 (i), and the u^2 rows give part (ii). An
independent argument is recorded in its docstring: for maximum degree
three, X_v^4 = 0, and after the scaling φ_e ↦ λφ_e every f without
linear term is proportional to u^2, u^2+u^3 or u^3 modulo u^4, or is zero.

## Notes

* `tables/tables345_clean.tex` is the table source of the paper, from which
  `check_tables.py` and `complete_profiles.py` read the graphs.
* `checks/neck_direct.py` computes modulo the prime 2³¹−1. It is a
  consistency check of Example 4.17, which is proved with the exact
  boundary states of `neck_states.py`.
* The cycle series are computed modulo two primes. Theorems 7.1 and 7.3
  are proved in the paper; the agreement of the data with them is a check
  of the program, and the agreement with Conjecture 7.8 is evidence.

## Recorded run

The recorded output was produced on 25 October 2026 on a Linux machine
with four cores (Intel Xeon, 2.1 GHz), with Python 3.13.16, NumPy 2.5.3,
SymPy 1.14.0 (with python-flint 0.9.0, so that SymPy uses its flint ground
types), networkx 3.7, g++ 13.3.0 and GMP 6.3.0. `run_all.log` lists every
command once, with its time. The commands were run in several pieces, at
most two at a time, so the times include waiting for a core. They add up
to about 3.1 hours: certificates 1 min, checks 7 min, tables 3 min,
search 26 min, strata 39 min, cycles 49 min, generic_exp 61 min.

Every output agrees with the output recorded before (17–20 October 2026),
where there was one, apart from timings, with these exceptions:

- `tables/check_tables_output.txt` calls the tables B1–B4 instead of 1–4;
  the entries are the same.
- `tables/complete_profiles.log` lists the values of λ for each graph in a
  fixed order; before, the order depended on Python's hashing.
- `checks/generic_exp/results/six_6.log` and `six_6_15.json` list the
  graphs in a fixed order; before, the order depended on which process
  finished first.
- `checks/generic_exp/results/classify_p101.txt` and `classify_p61.txt`
  also list the exceptional points with α = 0, which were already in the
  JSON files.
- `strata/d4_cells.py` first stopped with an error: with SymPy's flint
  ground types, plain rational numbers also have a `numer` method.
  `certificates/boundary_local.py` now tests explicitly for elements of a
  rational function field, and the output of `d4_cells.py` is unchanged.

`checks/generic_exp/low_certificate.py` (Proposition 6.26) and the check
of Proposition 7.7 in `cycles/cycle_formulas.py` are new. The cycle
series agree with an earlier run of `cycle_hilbert.py` with the same
primes.

On 27 October 2026 the statement numbers in the scripts and READMEs were
updated to the first arXiv version of the paper, and
`cycles/results/cycle_formulas.log`, which prints statement numbers, was
regenerated. Apart from these numbers it is identical to the recorded run.

## Licence and citation

The code and data are released under the MIT licence (`LICENSE`). If you
use them, please cite the paper; `CITATION.cff` has the details.
