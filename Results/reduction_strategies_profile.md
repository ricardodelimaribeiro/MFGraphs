# Disjunct Reduction Strategy Profiling

Generated: 2026-03-31
Mathematica: 14.3.0 for Mac OS X ARM (64-bit)

## Context

After `DNFReduce` produces a disjunctive normal form, some cases have many disjuncts.
This profile tests three strategies for further reducing the number of disjuncts.

## Cases with DNF Disjuncts > 1

| Case | Network | DNF Disjuncts | DNF Leaves | DNF Time (s) |
|------|---------|--------------|------------|-------------|
| Multi_20 (case 20) | Multi-entry (2 in, 3 out) | 8 | 669 | 0.078 |
| Braess congest | Braess paradox + congestion | 179 | 18,944 | 1.559 |
| Attraction_switch (case 11) | Attraction + 4 switching costs | 2,977 | 1,199,868 | 4.380 |

## Strategies

1. **Subsumption Pruning**: Pairwise `Reduce[Implies[A, B], Reals]` — if A implies B, drop B.
   O(n^2) checks worst case, but early pruning reduces actual work.

2. **Full Reduce**: `Reduce[Or @@ branches, Reals]` — let Mathematica merge everything at once.

3. **Group & Merge**: Group branches by equality skeleton, then `Reduce` the inequality
   disjunction within each group. Cannot detect cross-group subsumption.

## Results

### Multi_20 (8 disjuncts)

| Strategy | Result Disjuncts | Time (s) | Notes |
|----------|-----------------|----------|-------|
| Subsumption | **1** | 0.11 | All 7 branches subsumed by branch 1; only 7 checks needed |
| Full Reduce | **1** | 0.02 | 47 leaves |
| Group & Merge | 4 | 0.002 | 4 equality groups of 2 branches each |

### Braess congest (179 disjuncts)

| Strategy | Result Disjuncts | Time (s) | Notes |
|----------|-----------------|----------|-------|
| Subsumption | **1** | 2.52 | Branches 1-90 subsumed by branch 1, then branch 1 subsumed by branch 91, branches 92-179 subsumed by branch 91; 178 checks |
| Full Reduce | **1** | 0.14 | 79 leaves |
| Group & Merge | 126 | 0.01 | 126 equality groups; merged within groups (179->126) but can't detect cross-group subsumption |

### Attraction_switch (2,977 disjuncts)

| Strategy | Result Disjuncts | Time (s) | Notes |
|----------|-----------------|----------|-------|
| Subsumption (sample of 50) | **2** | 35.7 | 48 of 50 branches pruned in 84 checks; extrapolating, full run may yield ~2-10 disjuncts |
| Full Reduce | TIMED OUT | >120 | 1.2M leaves too large for Reduce |
| Group & Merge | 2,237 | 129 | 2,977->2,237; many per-group Reduce timeouts; cannot detect cross-group subsumption |

## Key Findings

1. **Massive redundancy exists**: Both Multi_20 (8->1) and Braess (179->1) collapsed entirely
   to a single disjunct. DNFReduce produces correct but highly redundant output.

2. **Full Reduce is fastest**: 0.02s and 0.14s respectively — much faster than subsumption
   pruning, and achieves the same result.

3. **Subsumption works but is slower**: Linear scan suffices when one dominant branch exists
   (all others are subsumed), but the pairwise check overhead is significant.

4. **Group & Merge is incomplete**: Fast but only merges within equality groups. Cannot
   detect that branches with different equality skeletons are logically subsumed.

5. **Implication for the solver**: A post-DNFReduce `Reduce[Or @@ result, Reals]` step
   could dramatically simplify the output. For Braess, this turns 179 disjuncts (18,944 leaves)
   into 1 disjunct (79 leaves) in 0.14s — a 240x reduction in expression size.

6. **Subsumption scales to large cases**: Even for Attraction_switch (2,977 disjuncts),
   sampling 50 branches found that 48 were redundant (96% pruning rate). A full subsumption
   pass would be O(n) per surviving branch, not O(n^2), when few branches survive.

7. **Recommended approach by case size**:
   - Small DNF output (<100 disjuncts): `Reduce[Or @@ result, Reals]` — fast, optimal
   - Medium (100-1000): Subsumption pruning — linear-ish when most branches are redundant
   - Large (>1000): Subsumption with early termination, then `Reduce` on the survivors

## Equivalence Check Profiling

Separate profiling of `Reduce[bcResult <=> dnfResult, Reals]` for the three EQUIV_TIMEOUT cases:

| Case | BC Leaves | DNF Leaves | Equiv Check (600s) | Result |
|------|-----------|------------|-------------------|--------|
| Braess congest | 166,035,457 | 18,944 | 601s | EQUIV_TIMEOUT |

The equivalence check times out because BC output is enormous (166M leaves).
After applying Full Reduce to the DNF output (1 disjunct, 79 leaves), a direct
equivalence check against the reduced BC output would be far more tractable.

## Benchmark Suite Results (post-merge master)

### Small tier: 21/21 OK

All 7 cases (1-6, 27) x 3 solvers (CriticalCongestion, NonLinearSolver, Monotone) passed.

### Medium tier (partial — killed during case 8 Monotone)

| Case | D2E | CriticalCongestion | NonLinearSolver | Monotone |
|------|-----|-------------------|----------------|----------|
| 7 | OK | OK | OK | OK |
| 8 | OK | OK | OK | (running when killed) |
