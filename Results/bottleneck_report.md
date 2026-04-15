# MFGraphs Bottleneck Analysis Report

Generated: Thu 16 Apr 2026 00:24:57
Mathematica: 14.3.0 for Mac OS X ARM (64-bit) (July 8, 2025)

## Summary

| Label | Case | Total (s) | DNFReduce (s) | DNFReduce Calls | DNFReduce Depth | TripleStep (s) | TripleStep Calls | GradProj (s) | Status |
|---|---|---|---|---|---|---|---|---|---|
| small_simple | 3 | 0.592 | 0.000 | 0 | 0 | 0.429 | 3 | 0.000 | OK |
| medium_nosw | 7 | 0.078 | 0.000 | 0 | 0 | 0.003 | 3 | 0.000 | OK |
| medium_sw | 8 | 0.290 | 0.008 | 1 | 1 | 0.003 | 9 | 0.000 | OK |
| medium_attraction | 12 | 0.127 | 0.000 | 0 | 0 | 0.003 | 3 | 0.000 | OK |
| medium_triangle | 14 | 0.086 | 0.008 | 1 | 1 | 0.003 | 9 | 0.000 | OK |
| large_multipath | 20 | 0.143 | 0.000 | 0 | 0 | 0.003 | 3 | 0.000 | OK |
| large_braess | Braess congest | 0.333 | 0.092 | 1 | 1 | 0.007 | 9 | 0.000 | OK |

## Detailed Profiles

### small_simple (case 3)

| Function | Total Time (s) | Call Count | Avg Time (s) |
|---|---|---|---|
| DNFReduce | 0.0000 | 0 | 0.000000 |
| DNFReduce_Solve | 0.0000 | 0 | 0.000000 |
| DNFReduce_Reduce | 0.0000 | 0 | 0.000000 |
| DNFReduce_Simplify | 0.0000 | 0 | 0.000000 |
| TripleStep | 0.4289 | 3 | 0.142957 |
| TripleClean | 0.4289 | 1 | 0.428946 |
| DNFSolveStep | 0.0000 | 0 | 0.000000 |
| GradientProjection | 0.0000 | 0 | 0.000000 |
| PseudoInverse | 0.0000 | 0 | 0.000000 |
| MFGSystemSolver | 0.0000 | 0 | 0.000000 |
| MFGPreprocessing | 0.0000 | 0 | 0.000000 |
| DNFReduce_MaxDepth | Missing[N/A] | 0 | Missing[N/A] |
| TripleClean_MaxIters | Missing[N/A] | 0 | Missing[N/A] |

### medium_nosw (case 7)

| Function | Total Time (s) | Call Count | Avg Time (s) |
|---|---|---|---|
| DNFReduce | 0.0000 | 0 | 0.000000 |
| DNFReduce_Solve | 0.0000 | 0 | 0.000000 |
| DNFReduce_Reduce | 0.0000 | 0 | 0.000000 |
| DNFReduce_Simplify | 0.0000 | 0 | 0.000000 |
| TripleStep | 0.0026 | 3 | 0.000858 |
| TripleClean | 0.0027 | 1 | 0.002651 |
| DNFSolveStep | 0.0000 | 0 | 0.000000 |
| GradientProjection | 0.0000 | 0 | 0.000000 |
| PseudoInverse | 0.0000 | 0 | 0.000000 |
| MFGSystemSolver | 0.0000 | 0 | 0.000000 |
| MFGPreprocessing | 0.0000 | 0 | 0.000000 |
| DNFReduce_MaxDepth | Missing[N/A] | 0 | Missing[N/A] |
| TripleClean_MaxIters | Missing[N/A] | 0 | Missing[N/A] |

### medium_sw (case 8)

| Function | Total Time (s) | Call Count | Avg Time (s) |
|---|---|---|---|
| DNFReduce | 0.0082 | 1 | 0.008200 |
| DNFReduce_Solve | 0.0000 | 0 | 0.000000 |
| DNFReduce_Reduce | 0.0000 | 0 | 0.000000 |
| DNFReduce_Simplify | 0.0000 | 0 | 0.000000 |
| TripleStep | 0.0027 | 9 | 0.000302 |
| TripleClean | 0.0029 | 5 | 0.000588 |
| DNFSolveStep | 0.0111 | 1 | 0.011087 |
| GradientProjection | 0.0000 | 0 | 0.000000 |
| PseudoInverse | 0.0000 | 0 | 0.000000 |
| MFGSystemSolver | 0.0000 | 0 | 0.000000 |
| MFGPreprocessing | 0.0000 | 0 | 0.000000 |
| DNFReduce_MaxDepth | Missing[N/A] | 1 | Missing[N/A] |
| TripleClean_MaxIters | Missing[N/A] | 0 | Missing[N/A] |

### medium_attraction (case 12)

| Function | Total Time (s) | Call Count | Avg Time (s) |
|---|---|---|---|
| DNFReduce | 0.0000 | 0 | 0.000000 |
| DNFReduce_Solve | 0.0000 | 0 | 0.000000 |
| DNFReduce_Reduce | 0.0000 | 0 | 0.000000 |
| DNFReduce_Simplify | 0.0000 | 0 | 0.000000 |
| TripleStep | 0.0030 | 3 | 0.000990 |
| TripleClean | 0.0031 | 1 | 0.003076 |
| DNFSolveStep | 0.0000 | 0 | 0.000000 |
| GradientProjection | 0.0000 | 0 | 0.000000 |
| PseudoInverse | 0.0000 | 0 | 0.000000 |
| MFGSystemSolver | 0.0000 | 0 | 0.000000 |
| MFGPreprocessing | 0.0000 | 0 | 0.000000 |
| DNFReduce_MaxDepth | Missing[N/A] | 0 | Missing[N/A] |
| TripleClean_MaxIters | Missing[N/A] | 0 | Missing[N/A] |

### medium_triangle (case 14)

| Function | Total Time (s) | Call Count | Avg Time (s) |
|---|---|---|---|
| DNFReduce | 0.0084 | 1 | 0.008368 |
| DNFReduce_Solve | 0.0000 | 0 | 0.000000 |
| DNFReduce_Reduce | 0.0000 | 0 | 0.000000 |
| DNFReduce_Simplify | 0.0000 | 0 | 0.000000 |
| TripleStep | 0.0030 | 9 | 0.000337 |
| TripleClean | 0.0033 | 5 | 0.000655 |
| DNFSolveStep | 0.0092 | 1 | 0.009202 |
| GradientProjection | 0.0000 | 0 | 0.000000 |
| PseudoInverse | 0.0000 | 0 | 0.000000 |
| MFGSystemSolver | 0.0000 | 0 | 0.000000 |
| MFGPreprocessing | 0.0000 | 0 | 0.000000 |
| DNFReduce_MaxDepth | Missing[N/A] | 1 | Missing[N/A] |
| TripleClean_MaxIters | Missing[N/A] | 0 | Missing[N/A] |

### large_multipath (case 20)

| Function | Total Time (s) | Call Count | Avg Time (s) |
|---|---|---|---|
| DNFReduce | 0.0000 | 0 | 0.000000 |
| DNFReduce_Solve | 0.0000 | 0 | 0.000000 |
| DNFReduce_Reduce | 0.0000 | 0 | 0.000000 |
| DNFReduce_Simplify | 0.0000 | 0 | 0.000000 |
| TripleStep | 0.0029 | 3 | 0.000957 |
| TripleClean | 0.0030 | 1 | 0.002987 |
| DNFSolveStep | 0.0000 | 0 | 0.000000 |
| GradientProjection | 0.0000 | 0 | 0.000000 |
| PseudoInverse | 0.0000 | 0 | 0.000000 |
| MFGSystemSolver | 0.0000 | 0 | 0.000000 |
| MFGPreprocessing | 0.0000 | 0 | 0.000000 |
| DNFReduce_MaxDepth | Missing[N/A] | 0 | Missing[N/A] |
| TripleClean_MaxIters | Missing[N/A] | 0 | Missing[N/A] |

### large_braess (case Braess congest)

| Function | Total Time (s) | Call Count | Avg Time (s) |
|---|---|---|---|
| DNFReduce | 0.0924 | 1 | 0.092364 |
| DNFReduce_Solve | 0.0000 | 0 | 0.000000 |
| DNFReduce_Reduce | 0.0000 | 0 | 0.000000 |
| DNFReduce_Simplify | 0.0000 | 0 | 0.000000 |
| TripleStep | 0.0066 | 9 | 0.000732 |
| TripleClean | 0.0070 | 5 | 0.001404 |
| DNFSolveStep | 0.1082 | 1 | 0.108190 |
| GradientProjection | 0.0000 | 0 | 0.000000 |
| PseudoInverse | 0.0000 | 0 | 0.000000 |
| MFGSystemSolver | 0.0000 | 0 | 0.000000 |
| MFGPreprocessing | 0.0000 | 0 | 0.000000 |
| DNFReduce_MaxDepth | Missing[N/A] | 1 | Missing[N/A] |
| TripleClean_MaxIters | Missing[N/A] | 0 | Missing[N/A] |

## Key Observations

- **Highest DNFReduce call count**: medium_sw with 1 calls, taking 0.008s
- **Highest TripleStep call count**: medium_sw with 9 calls, taking 0.003s

## Optimization Targets

1. **DNFReduce Solve Memoization**: Cache results of `Solve` and `Reduce` to avoid redundant symbolic solves across Or-branches.
2. **DNFReduce Branch Pruning**: Early-exit when `xp === False` to avoid processing inconsistent branches.
3. **PseudoInverse Caching**: Cache the PseudoInverse result in `GradientProjection` to avoid O(n^3) recomputation per ODE step.
4. **M Interpolation**: Pre-compute `M[j, x, edge]` on a grid and interpolate instead of calling `FindRoot` per point.