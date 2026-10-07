# Block-ordering sweep — results summary

Data: [`../results/block_ordering/benchmark.csv`](../results/block_ordering/benchmark.csv),
produced by [`benchmark_orderings.wls`](benchmark_orderings.wls) on 2026-10 after the
issue #222 fix made `DisjunctOrdering` live. The pre-fix sweep was retracted in 8b88781.

## Setup

- 37 registered scenarios (Grid1020 skipped, kernel crash), 5 modes, 60 s cap per run.
- Modes: `Lexicographic` (default path, routes through `dnfReduce`), `Lex-BranchState`
  (default ordering forced onto the branch-state path, the same-path baseline), and
  `Block-Vertex`, `Block-Edge`, `Block-SCC` (always branch-state).
- One untimed warm-up solve per scenario, then 3 timed repetitions per (scenario, mode),
  interleaved across modes with the order rotated per rep and per scenario. A mode that
  times out is not re-run in later reps.
- Every `OK` result is validity-checked with `isValidSystemSolution` (10 s cap).

457 rows; all 408 `OK` rows validate.

## Measurement quality

| Check | Value |
|---|---:|
| Median rep1 / rep3 wall, same (scenario, mode) | 1.01 |
| Median coefficient of variation across 3 reps | 0.01 |
| Normalised wall by run position (slots 1–15) | 1.02, 1.01, then 1.00 throughout |

The kernel warm-up effect that produced the retracted July numbers is absent.

## Finding

Geometric-mean wall-time ratios over the 26 scenarios that complete in every mode
(ratio < 1 means faster):

| Comparison | Ratio |
|---|---:|
| Lex-BranchState / Lexicographic (default path) | 2.80 |
| Block-Vertex / Lexicographic (default path) | 2.94 |
| Block-Edge / Lexicographic (default path) | 2.86 |
| Block-SCC / Lexicographic (default path) | 2.97 |
| Block-Vertex / Lex-BranchState (same path) | 1.05 |
| Block-Edge / Lex-BranchState (same path) | 1.02 |
| Block-SCC / Lex-BranchState (same path) | 1.06 |

Within the branch-state path the three Block orderings are indistinguishable from the
lexicographic default. The branch-state path itself is about 3× slower than the default
DNF path. The gap is extreme on Grid0303: 0.03 s on the default path against about 14 s on
every branch-state mode.

Per-scenario medians (s) where any mode exceeds 50 ms:

| Scenario | Lexicographic | Lex-BranchState | Block-Vertex | Block-Edge | Block-SCC |
|---|---:|---:|---:|---:|---:|
| Achdou_2023_junction | 0.069 | 0.115 | 0.114 | 0.118 | 0.117 |
| Big_Braess_congest | 0.140 | 1.157 | 0.983 | 1.110 | 0.977 |
| Braess_congest | 0.137 | 1.102 | 0.895 | 0.777 | 0.899 |
| Grid0303 | 0.030 | 14.090 | 13.921 | 15.349 | 13.970 |
| Inconsistent_attraction_shortcut | 0.341 | 1.714 | 2.685 | 2.648 | 2.663 |
| case_11 | 0.324 | 1.637 | 2.553 | 2.511 | 2.571 |
| case_20 | 0.083 | 0.626 | 0.370 | 0.538 | 0.362 |

## Timeouts (60 s)

| Mode | OK | Timeout |
|---|---:|---:|
| Lexicographic | 90 | 7 |
| Lex-BranchState | 84 | 9 |
| Block-Vertex / Block-Edge / Block-SCC | 78 each | 11 each |

- Camilli_2015_general, Grid0404, Jamaratv9: finish on the default path, time out on every
  branch-state mode.
- case_23: times out on the default path, finishes on Lex-BranchState (matches the
  BENCHMARKS.md note that forced lexicographic solving rescues this case).
- case_22: the only scenario where ordering matters. Both lexicographic modes finish; all
  three Block modes time out.
- Grid0505, Grid0707, Grid0710, Grid1010, HRF_Scenario_1, case_21: time out in every mode.

## Conclusion

`DisjunctOrdering` is live and harmless but earns no default change. Block orderings do not
beat lexicographic ordering on the branch-state path, and the branch-state path does not beat
the default DNF path. The option stays opt-in.
