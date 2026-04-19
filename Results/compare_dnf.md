# BooleanConvert vs DNFReduce Comparison

Generated: Sun 19 Apr 2026 10:16:34
Mathematica: 14.3.0 for Mac OS X ARM (64-bit) (July 8, 2025)

## Summary

| Case | Input Disjuncts | BC Time (s) | DNF Time (s) | Pipe Time (s) | **Speedup (Pipe)** | BC Disjuncts | DNF Disjuncts | Pipe Disjuncts | BC Leaves | DNF Leaves | Pipe Leaves | Equiv |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 7 | 1 | 0.001915 | 0.024675 | 0.252005 | 10.2906x | 12 | 1 | 1 | 1491 | 20 | 18 | ✓ |

## Interpretation

- **Speedup (Pipe)**: How much faster `DNFReduce` is compared to the `BooleanConvert` + `Reduce` pipeline.
- **BC Time**: Time for pure boolean conversion.
- **Pipe Time**: Time for boolean conversion followed by algebraic simplification (`Reduce`). This is the fair baseline for `DNFReduce`.
- **Disjuncts**: Fewer disjuncts = simpler output (fewer solution branches)
- **Leaves**: Smaller LeafCount = more compact expression

## Method details

- **BooleanConvert**: Mathematica's built-in `BooleanConvert[expr, "DNF"]`. Treats the expression as a pure boolean formula.
- **DNFReduce**: Our custom solver that interleaves algebraic solving (substitution and pruning) with DNF conversion.
- **Pipe (Method 3)**: `BooleanConvert` followed by `Reduce[..., Reals]`. This isolates whether the benefit of `DNFReduce` comes from the simplification itself or the interleaved strategy.

## Key differences

DNFReduce is *not* a pure boolean conversion. It interleaves algebraic solving with DNF conversion:
- Equalities are solved and substituted, reducing the number of free variables
- Branches that become False after substitution are pruned immediately
- The output is a simplified DNF where each disjunct has had equalities resolved

BooleanConvert treats the entire expression as a boolean formula and does not perform algebraic simplification. Method 3 (Pipe) adds that simplification as a post-processing step.