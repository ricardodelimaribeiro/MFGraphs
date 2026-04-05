# BooleanConvert vs DNFReduce Comparison

Generated: Thu 2 Apr 2026 14:11:29
Mathematica: 14.3.0 for Mac OS X ARM (64-bit) (July 8, 2025)

## Summary

| Case | Input Disjuncts | BC Time (s) | DNF Time (s) | Speedup | BC Disjuncts | DNF Disjuncts | BC Leaves | DNF Leaves | Equivalence |
|---|---|---|---|---|---|---|---|---|---|
| 7 | - | - | - | - | - | - | - | - | TRIVIAL |
| 8 | 1 | 0.000104 | 0.006419 | 0.0162019x | 4 | 1 | 313 | 16 | EQUIVALENT |
| 12 | 1 | 0.000032 | 0.00037 | 0.0864865x | 1 | 1 | 23 | 6 | EQUIVALENT |
| 11 | 1 | 0.023437 | 0.051185 | 0.457888x | 12 | 20 | 5763 | 6869 | EQUIV_TIMEOUT |
| 14 | 2 | 0.000033 | 0.001142 | 0.0288967x | 2 | 1 | 82 | 6 | EQUIVALENT |
| Braess congest | 1 | 0.000202 | 0.006403 | 0.0315477x | 4 | 1 | 649 | 67 | EQUIVALENT |
| Paper example | 1 | 0.000035 |      -6
6. 10 | 5.83333x | 1 | 1 | 3 | 3 | EQUIVALENT |
| 20 | 1 | 0.000254 | 0.006305 | 0.0402855x | 8 | 4 | 1289 | 240 | EQUIVALENT |

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