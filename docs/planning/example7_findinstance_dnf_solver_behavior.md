# Example 7 Solver Behavior: `dnfReduceSystem` vs `findInstanceSystem`

Date: 2026-04-27

This note records an exploratory comparison after adding `findInstanceSystem`.
The goal was to compare solver behavior on edge flows only, ignoring transition
flows, and to identify which `u` values differ when `findInstanceSystem`
chooses a concrete feasible point.

Update 2026-04-28: `dnfReduceSystem` now harvests equality solutions from the
surviving final DNF branches. A later 2026-04-28 change relaxes each physical
edge Hamiltonian equation behind a zero-total-flow alternative:

```wolfram
j[a, b] + j[b, a] == 0 || j[a, b] - j[b, a] + u[a, b] - u[b, a] == 0
```

That means current `dnfReduceSystem` may return branch-family associations where
earlier runs returned one concrete harvested branch.

## Benchmark

Environment: Wolfram 14.3.0 for Mac OS X ARM, local `wolframscript`.
Timeout per solve: 60 seconds.

| Case | `dnfReduceSystem` repeated | `findInstanceSystem` repeated | Validity |
|---|---:|---:|---|
| `chain-2v` | 0.60 ms | 0.59 ms | both valid |
| `chain-3v-1exit` | 0.92 ms | 0.90 ms | both valid |
| `chain-3v-2exit` | 1.88 ms | 2.51 ms | both valid |
| `example-7` | 3.52 ms | 21.61 ms | both valid |
| `chain-5v-1exit` | 1.74 ms | 1.71 ms | both valid |
| `example-12` | 258.69 ms | timed out at 60s | `dnfReduceSystem` valid; `findInstanceSystem` timed out |

Conclusion from this run: `findInstanceSystem` is useful as a concretizer/oracle
for one feasible point, but `dnfReduceSystem` is faster on branching examples.
After the 2026-04-28 harvesting fix, `dnfReduceSystem` also returns concrete
rules when the final DNF residual has a unique feasible branch.

## Edge-Flow Comparison

Transition flows `j[_, _, _]` were intentionally ignored. Edge flows were
compared using `systemData[sys, "Js"]`.

Before the harvesting fix, all non-timeout cases except `example-7` had matching
returned edge-flow values. For `example-7`, the difference was not a
contradiction: `dnfReduceSystem` left `j[2, 4]` free, while
`findInstanceSystem` chose:

```wolfram
j[2, 4] -> 45
j[2, 3] -> 55
j[3, auxExit3] -> 55
j[4, auxExit4] -> 45
```

The corresponding `dnfReduceSystem` parametric expressions are:

```wolfram
j[3, auxExit3] -> 100 - j[2, 4]
j[4, auxExit4] -> j[2, 4]
j[2, 3] -> 100 - j[2, 4]
```

Current behavior after zero-flow Hamiltonian relaxation: `dnfReduceSystem`
returns the common parametric rules and a residual branch family. The concrete
values above are still one branch, but no longer the only branch represented by
the public result.

## Flow Oracle Check

An oracle check asserted the edge-flow values chosen by `findInstanceSystem`
back into the accumulated full system and called `FindInstance` for feasibility.

| Case | Flow assignments | Oracle status |
|---|---:|---|
| `chain-2v` | 4 | feasible |
| `chain-3v-1exit` | 6 | feasible |
| `chain-3v-2exit` | 7 | feasible |
| `example-7` | 9 | feasible |
| `chain-5v-1exit` | 10 | feasible |
| `example-12` | 0 | skipped because `findInstanceSystem` timed out |

This supports interpreting the `example-7` flow difference as one feasible
instantiation of the `dnfReduceSystem` family, not as a solver disagreement.

## Different `u` Values

For `chain-3v-2exit`, differences are due to unassigned symbolic values in the
DNF result:

```wolfram
u[2, 3]: dnf unassigned, findInstanceSystem -> 0
u[3, 2]: dnf -> u[2, 3], findInstanceSystem -> 0
```

Before the harvesting fix, `dnfReduceSystem` left several `example-7` value
variables symbolic, while `findInstanceSystem` chose the following concrete
values:

```wolfram
u[auxEntry1, 1]: dnf -> 100 + u[1, 2], findInstanceSystem -> 155
u[1, 2]:        dnf unassigned,       findInstanceSystem -> 55
u[2, 3]:        dnf unassigned,       findInstanceSystem -> 0
u[2, 4]:        dnf unassigned,       findInstanceSystem -> 10
u[2, 1]:        dnf -> 100 + u[1, 2], findInstanceSystem -> 155
u[3, 2]:        dnf -> 100 - j[2, 4] + u[2, 3], findInstanceSystem -> 55
u[4, 2]:        dnf -> j[2, 4] + u[2, 4],       findInstanceSystem -> 55
```

Substituting the concrete choices

```wolfram
j[2, 4] -> 45
u[1, 2] -> 55
u[2, 3] -> 0
u[2, 4] -> 10
```

into the DNF parametric expressions reproduces the `findInstanceSystem` values.

Current behavior after zero-flow Hamiltonian relaxation: `dnfReduceSystem`
keeps the concrete `u` values above inside one residual branch instead of
collapsing the whole result to that branch.

## Chain-3v-2exit Underdetermined Residual Check

Run on 2026-04-28:

```wolfram
s = gridScenario[{3}, {{1, 5}}, {{2, 10}, {3, 0}}];
sys = makeSystem[s];
```

This is the `chain-3v-2exit` shape with inflow smaller than the exit value at
vertex 2 and zero exit value at vertex 3.

The raw DNF result, before final branch harvesting, is still useful for studying
underdetermined solver behavior:

```wolfram
rulesAcc:
  j[2, "auxExit2"] -> 5 - j[2, 3]
  j[3, "auxExit3"] -> j[2, 3]
  u[2, 1] -> 5 + u[1, 2]
  u[3, 2] -> j[2, 3] + u[2, 3]
  u["auxEntry1", 1] -> 5 + u[1, 2]

remaining vars:
  {j[2, 3], u[1, 2], u[2, 3]}

raw DNF:
  (5 - j[2, 3] == 0 && 5 == 5 + u[2, 3] && -5 + u[1, 2] == 0)
  ||
  (10 <= u[2, 3] && u[2, 3] <= 10 && u[2, 3] <= 10 &&
   u[2, 3] <= 0 && u[1, 2] == 10 && j[2, 3] == 0)
```

Using the old `parseReduceResult` path on that raw DNF keeps the result as an
underdetermined association with the residual above. The current public
`dnfReduceSystem` discards the contradictory second branch and harvests the
first branch into concrete rules:

```wolfram
j[2, "auxExit2"] -> 0
j[3, "auxExit3"] -> 5
j[2, 3] -> 5
u[1, 2] -> 5
u[2, 3] -> 0
u[3, 2] -> 5
```

`findInstanceSystem` returns one concrete branch. After the zero-flow Hamiltonian
relaxation, public `dnfReduceSystem` preserves multiple branches for this case:

```wolfram
(j[2, "auxExit2"] == 5 && u[2, 1] == 15 && u[3, 2] == 10 && u[2, 3] <= 0)
||
(j[2, "auxExit2"] == 0 && u[2, 1] == 10 && u[2, 3] == 0 && u[3, 2] == 5)
```

The first branch is the newly admitted zero-flow-edge branch where all mass exits
at vertex 2 and `u[2,3]` is bounded but not pinned by the edge equation. The
second branch is the previous all-to-vertex-3 solution.
