# Cycles (Section 7)

Requirements: Python 3 with NumPy. Run the commands from this directory.

| Script | What it computes | Arithmetic | Output | Time |
|---|---|---|---|---|
| `cycle_hilbert.py 3 17 --primes 2 --out results/cycle_series.json` | the Hilbert series H_n of the cycle C_n for 3 ≤ n ≤ 17 (data for Conjecture 7.8) | ranks modulo two primes p ≡ 1 mod n (evidence) | `results/cycle_series.json`, `results/cycle_hilbert.log` | about 50 min in the recorded run, of which 41 min for n = 17 and 2 min for n ≤ 15 |
| `cycle_formulas.py` | Theorems 7.1 and 7.3, Conjecture 7.8, the recurrence (7.28), the closed formulas (7.29)–(7.30) and the initial values of Proposition 7.9 on the data; the trace identities of the transfer matrix (4.13) and Remark 7.10 for n ≤ 40; Proposition 7.7 (the weighted trace and independent sets) for 3 ≤ n ≤ 18 | exact (integers) | `results/cycle_formulas.log` | seconds |
| `annihilators.py` | Lemma 7.4 (r ≤ 5), (7.17) and F_(m+2) = CC_(2m+1) for odd cycles (n ≤ 9), Proposition 7.5 and F_(m+1) = CC_(2m) for even cycles (n ≤ 10), the six annihilators of C_6 in Section 7.5 | certificate: exact pairings, ranks modulo 1048573 (see below) | `results/annihilators.log` | 10 s |

## The series H_n

By Proposition 3.14 every f with a nonzero quadratic coefficient gives the
same series, and by Proposition 7.12 the space F_d is spanned by the product
vectors E_m = Π_i (1 + (m_i − m_(i+1)) φ_i) with m ∈ {0,1,2}^n and
|m| ≤ d (the integer matrix A_(n,d) of (7.31), up to the reflection
i ↦ −i). `cycle_hilbert.py` adds these vectors degree by degree, testing
only E_(m+δ_i) for the m that gave new pivots in the previous degree,
and splits the edge algebra into the eigenspaces of the rotation over
F_p with p ≡ 1 (mod n); the reflection φ_i ↦ −φ_(−1−i) identifies the
eigenspaces j and −j, so about n/2 blocks of size 2^n/n are reduced.
The method is described in the docstring.

Ranks modulo p are lower bounds for ranks over Q, so each computed
coefficient sequence is evidence, not a certificate. The two primes give
the same series for every n, the series add up to 2^n − 1 (the dimension
of the algebra, (7.9)), and `cycle_formulas.py` confirms Theorems 7.1
(n ≥ 7) and 7.3, which are proved in the paper, and the conjectured
formulas, which are not:

| n | H_n |
|---|---|
| 3 | 1, 3, 2, 1 |
| 4 | 1, 4, 7, 3 |
| 5 | 1, 5, 14, 10, 1 |
| 6 | 1, 6, 20, 31, 5 |
| 7 | 1, 7, 27, 63, 28, 1 |
| 8 | 1, 8, 35, 96, 106, 9 |
| 9 | 1, 9, 44, 138, 243, 75, 1 |
| 10 | 1, 10, 54, 190, 405, 346, 17 |
| 11 | 1, 11, 65, 253, 627, 891, 198, 1 |
| 12 | 1, 12, 77, 328, 921, 1620, 1103, 33 |
| 13 | 1, 13, 90, 416, 1300, 2691, 3159, 520, 1 |
| 14 | 1, 14, 104, 518, 1778, 4200, 6237, 3466, 65 |
| 15 | 1, 15, 119, 635, 2370, 6258, 11070, 10935, 1363, 1 |
| 16 | 1, 16, 135, 768, 3092, 8992, 18288, 23328, 10786, 129 |
| 17 | 1, 17, 152, 918, 3961, 12546, 28662, 44064, 37179, 3570, 1 |

## The annihilator certificates

For each statement, `annihilators.py` lists all products S of at most d
generators (exactly, in integers) and the proposed spanning set A of the
annihilator under the pairing (7.8). It checks that every pairing
⟨a, s⟩ is zero, exactly, and computes the ranks r_S and r_A modulo
p = 1048573. Ranks modulo p are at most ranks over Q, and the pairing is
nondegenerate, so r_S + r_A = 2^N gives

    2^N = r_S + r_A ≤ rk_Q S + rk_Q A ≤ 2^N.

Hence dim F_d = r_S and F_d^⊥ = span A exactly. The script also checks
r_A against the dimension stated in the paper (3^(r−1), 1 + 3^(r−2),
2^(m−1) + 2, …).
