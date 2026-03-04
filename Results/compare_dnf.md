# BooleanConvert vs DNFReduce Comparison

Generated: Wed 4 Mar 2026 08:12:13
Mathematica: 14.3.0 for Mac OS X ARM (64-bit) (July 8, 2025)

## Summary

| Case | Input Disjuncts | BC Time (s) | DNF Time (s) | Speedup | BC Disjuncts | DNF Disjuncts | BC Leaves | DNF Leaves | Equivalence |
|---|---|---|---|---|---|---|---|---|---|
| 7 | - | - | - | - | - | - | - | - | ERROR |
| 8 | - | - | - | - | - | - | - | - | ERROR |
| 12 | - | - | - | - | - | - | - | - | ERROR |
| 11 | - | - | - | - | - | - | - | - | ERROR |
| 14 | - | - | - | - | - | - | - | - | ERROR |
| Braess congest | - | - | - | - | - | - | - | - | ERROR |
| Paper example | - | - | - | - | - | - | - | - | ERROR |
| 20 | - | - | - | - | - | - | - | - | ERROR |

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