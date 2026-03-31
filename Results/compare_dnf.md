# BooleanConvert vs DNFReduce Comparison

Generated: Tue 10 Mar 2026 10:50:50
Mathematica: 14.3.0 for Mac OS X ARM (64-bit) (July 8, 2025)

## Summary

| Case | Input Disjuncts | BC Time (s) | DNF Time (s) | Speedup | BC Disjuncts | DNF Disjuncts | BC Leaves | DNF Leaves | Equivalence |
|---|---|---|---|---|---|---|---|---|---|
| 7 | 1 | 0.000175 | 0.002347 | 0.0745633x | 3 | 1 | 158 | 10 | EQUIVALENT |
| 8 | 1 | 0.000108 | 0.006674 | 0.0161822x | 8 | 1 | 1215 | 24 | EQUIVALENT |
| 12 | 1 | 0.71028 | 0.017474 | 40.6478x | 65536 | 1 | 44269569 | 55 | EQUIV_TIMEOUT |
| 11 | 1 | 48.5742 | 4.38017 | 11.0896x | 3359232 | 2977 | 6094206721 | 1199868 | EQUIV_TIMEOUT |
| 14 | 1 | 0.005027 | 0.018959 | 0.265151x | 384 | 1 | 139009 | 33 | EQUIVALENT |
| Braess congest | 1 | 1.35026 | 1.65131 | 0.817689x | 147456 | 179 | 166035457 | 18944 | EQUIV_TIMEOUT |
| Paper example | 1 | 0.00003 |      7. 10^(-6) | 4.28571x | 1 | 1 | 3 | 3 | EQUIVALENT |
| 20 | 1 | 0.005886 | 0.035008 | 0.168133x | 1024 | 8 | 469505 | 669 | EQUIVALENT |

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
