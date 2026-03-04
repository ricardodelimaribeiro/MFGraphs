# BooleanConvert vs DNFReduce Comparison

Generated: Wed 4 Mar 2026 12:03:47
Mathematica: 14.3.0 for Mac OS X ARM (64-bit) (July 8, 2025)

## Summary

| Case | Input Disjuncts | BC Time (s) | DNF Time (s) | Speedup | BC Disjuncts | DNF Disjuncts | BC Leaves | DNF Leaves | Equivalence |
|---|---|---|---|---|---|---|---|---|---|
| 7 | 1 | 0.000063 | 0.441266 | 0.000142771x | 3 | 1 | 158 | 10 | EQUIVALENT |
| 8 | 1 | 0.000094 | 0.00585 | 0.0160684x | 8 | 1 | 1215 | 24 | EQUIVALENT |
| 12 | 1 | 0.588052 | 0.015339 | 38.337x | 65536 | 1 | 44269569 | 55 | EQUIV_TIMEOUT |
| 11 | 1 | 40.4654 | 3.90456 | 10.3636x | 3359232 | 2977 | 6094206721 | 1199868 | EQUIV_TIMEOUT |
| 14 | 1 | 0.002211 | 0.014944 | 0.147952x | 384 | 1 | 139009 | 33 | EQUIVALENT |
| Braess congest | 1 | 1.241 | 1.37349 | 0.903533x | 147456 | 179 | 166035457 | 18944 | EQUIV_TIMEOUT |
| Paper example | 1 | 0.000023 |      -6
5. 10 | 4.6x | 1 | 1 | 3 | 3 | EQUIVALENT |
| 20 | 1 | 0.005263 | 0.031009 | 0.169725x | 1024 | 8 | 469505 | 669 | EQUIVALENT |

## Interpretation

- **Speedup > 1**: BooleanConvert is slower than DNFReduce (our method is faster)
- **Speedup < 1**: BooleanConvert is faster than DNFReduce
- **Disjuncts**: Fewer disjuncts = simpler output (fewer solution branches)
- **Leaves**: Smaller LeafCount = more compact expression

## Method details

- **BooleanConvert**: Mathematica's built-in `BooleanConvert[expr, "DNF"]`. Treats the expression as a pure boolean formula and converts to disjunctive normal form without exploiting equation structure.
- **DNFReduce**: Our custom solver that (1) solves equalities by substitution, (2) caches Solve/Reduce results, (3) prunes inconsistent branches early, (4) short-circuits when a branch already covers the expression.

## Key differences

DNFReduce is *not* a pure boolean conversion. It interleaves algebraic solving with DNF conversion:
- Equalities are solved and substituted, reducing the number of free variables
- Branches that become False after substitution are pruned immediately
- The output is a simplified DNF where each disjunct has had equalities resolved

BooleanConvert treats the entire expression as a boolean formula and does not perform algebraic simplification. This can lead to:
- More disjuncts (redundant branches not eliminated)
- Larger expressions (equalities not resolved)
- Potentially faster for small problems where the boolean structure is simple