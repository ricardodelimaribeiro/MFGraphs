# BooleanConvert vs DNFReduce Comparison

Generated: Mon 6 Apr 2026 10:56:49
Mathematica: 14.0.0 for Mac OS X ARM (64-bit) (December 13, 2023)

## Summary

| Case | Input Disjuncts | BC Time (s) | DNF Time (s) | Speedup | BC Disjuncts | DNF Disjuncts | BC Leaves | DNF Leaves | Equivalence |
|---|---|---|---|---|---|---|---|---|---|
| 7 | - | - | - | - | - | - | - | - | TRIVIAL |
| 8 | 1 | 0.000057 | 0.023938 | 0.00238115x | 4 | 1 | 313 | 16 | EQUIVALENT |
| 12 | 1 | 0.000041 | 0.000434 | 0.09447x | 1 | 1 | 23 | 6 | EQUIVALENT |
| 11 | 1 | 0.022839 | 0.045009 | 0.507432x | 12 | 20 | 5763 | 6869 | EQUIV_TIMEOUT |
| 14 | 2 | 0.000046 | 0.00117 | 0.0393162x | 2 | 1 | 82 | 6 | EQUIVALENT |
| Braess congest | 1 | 0.000163 | 0.005943 | 0.0274272x | 4 | 1 | 649 | 67 | EQUIVALENT |
| 20 | 1 | 0.000267 | 0.005819 | 0.0458842x | 8 | 4 | 1289 | 240 | EQUIVALENT |

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