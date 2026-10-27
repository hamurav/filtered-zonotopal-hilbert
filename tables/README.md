# Tables of Hilbert functions (formerly Appendix B)

These are Tables B1–B4, formerly Tables 12–15 (Appendix B) of the paper:
the Hilbert functions of the filtered zonotopal algebras of the connected graphs with four or five vertices and
maximum degree three (complete cubic profiles, Tables B1–B3) and of the
connected simple graphs with five vertices and maximum degree four (nine
parameters, Table B4).

| File | Contents |
|---|---|
| `appendix_B.pdf` | the tables with their introduction, as a 3-page document |
| `appendix_B.tex` | its source; compile with `pdflatex appendix_B.tex` (it inputs `tables345_clean.tex`) |
| `tables345_clean.tex` | the table source used by the paper, with the graphs drawn in TikZ |
| `tables.json`, `tables.csv` | all 189 entries: table, graph (macro name in the source), number of vertices, edge list, f, Hilbert function |
| `check_tables.py` | recomputes every entry from the edge lists read off the TikZ drawings |
| `complete_profiles.py` | computes the complete cubic profile of each graph of Tables B1–B3 (Remark 6.13): every special value of β/α² it finds is 2 or 2/3, the parabolas β = 2α² and 3β = 2α² of the last two columns |

Run from this directory (Python 3; SymPy for `complete_profiles.py`):

```
python3 check_tables.py tables345_clean.tex --certificates ../certificates --json tables.json
python3 complete_profiles.py tables345_clean.tex
```

`check_tables.py` computes every Hilbert function exactly, in rational
arithmetic, directly in the square-free edge algebra with the generators
Y_v = f(X_v) (`../certificates/direct_hilbert.py`); the output is
`check_tables_output.txt` (189 entries, 0 mismatches). With `--mod P` it
computes modulo the prime P instead, which takes about 20 seconds and
serves as a quick check. `complete_profiles.py` computes the Hilbert function over the
field Q(λ) for f = u + u² + λu³, records every factor inverted during the
elimination, and recomputes at each of its roots; the output is
`complete_profiles.log`.

In `tables.json`, vertices are numbered 1, …, n in the order of the
`\node[vertex]` commands of the drawing, an edge drawn twice (as two arcs)
is a pair of parallel edges, and f is given by its coefficients
{k: coefficient of u^k}.
