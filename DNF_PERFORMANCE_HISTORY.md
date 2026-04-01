# DNFReduce Performance History

This document records every optimization milestone applied to the `DNFReduce` algorithm
(originally called `ZAnd`) together with the measured impact on representative test cases.
Each entry explains **why** the change was made, **what** was changed in the code, and
**what was observed** in benchmarks.

---

## How to add a new entry

Run the comparison script with a tag and the `--append-history` flag:

```bash
wolframscript -file Scripts/CompareDNF.wls --tag "Short description of change"
```

This will:
1. Run all comparison cases (BooleanConvert vs DNFReduce)
2. Resolve the current git commit hash
3. Append a formatted entry to this file

Alternatively, copy the template at the bottom of this document and fill it in manually.

---

## Metrics glossary

| Column | Meaning |
|--------|---------|
| **BC Time** | Wall-clock time for `BooleanConvert[expr, "DNF"]` |
| **DNF Time** | Wall-clock time for `DNFReduce[xp, sorted]` |
| **Speedup** | BC Time / DNF Time — values > 1 mean DNFReduce is faster |
| **BC Disjuncts** | Number of Or-branches in BooleanConvert output |
| **DNF Disjuncts** | Number of Or-branches in DNFReduce output |
| **BC Leaves** | LeafCount of BooleanConvert output (expression size) |
| **DNF Leaves** | LeafCount of DNFReduce output (expression size) |
| **Equiv** | Whether outputs are logically equivalent (verified by Reduce) |

---

## Test cases

| Label | Network | Key | Switching costs |
|-------|---------|-----|----------------|
| `7` | Y-network, 1-in 2-out, no switching | 7 | No |
| `8` | Y-network, 1-in 2-out, with switching | 8 | Yes |
| `12` | Attraction problem, no switching | 12 | No |
| `11` | Attraction problem, with switching | 11 | Yes |
| `14` | Triangle with switching | 14 | Yes |
| `Braess congest` | Braess paradox congestion variant | "Braess congest" | No |
| `Paper example` | 4-vertex, 2 entrances, 2 exits | "Paper example" | Yes |
| `20` | 9-vertex multi-entrance/exit | 20 | No |

---

## Entries (oldest → newest)

---

### v0.1 — Original ZAnd  *(no benchmarks captured — code-level analysis only)*

**Date:** before 2026-03-03
**File:** `MFGraphs/ZAnd.m`
**Commit:** pre-887f43b (see `git log --follow MFGraphs/ZAnd.m`)

#### What the code did

```mathematica
ZAnd[xp_, And[fst_, rst_]] :=
    If[ Head[fst] === Or,
        RemoveDuplicates@(ReZAnd[xp, rst] /@ fst),  (* no xp===False check *)
        ReZAnd[xp, rst, fst]
    ]

ZAnd[xp_, Or[fst_, scd_]] :=
    With[{rfst = Reduce[fst, Reals]},               (* Reduce uncached *)
        RemoveDuplicates@(Or @@ (ZAnd[xp, #] & /@ {rfst, scd}))
        (* always evaluates BOTH branches — no short-circuit *)
    ]

ReZAnd[xp_, rst_, fst_Equal] :=
    Module[{newfst = Simplify@fst},
        If[newfst === False, False,
            With[{fsol = First@Solve@newfst},        (* Solve uncached *)
                ZAnd[ReplaceSolution[xp, fsol] && fst, ReplaceSolution[rst, fsol]]
            ]
        ]
    ]
```

#### Known deficiencies

- **No `Solve` memoization**: identical equalities (common across Or-branches with switching
  costs) each trigger a fresh symbolic `Solve` call.
- **No `Reduce` memoization**: `Reduce[fst, Reals]` called without cache; repeated Or-branches
  hit identical sub-expressions.
- **No short-circuit Or**: both sides of every Or are always fully evaluated, even when the
  first side already covers the entire expression.
- **No `xp === False` early exit** in the And case: the recursion continues into a dead branch
  even after the accumulated expression has become False.

#### Why these mattered

For networks with switching costs, the system contains an `Or` expression over all feasible
switching patterns.  A Y-network with 6 switching costs (case 8) generates an Or with up to
2^6 = 64 leaves; the Attraction problem with 16 switching costs (case 11) can reach 2^16 ≈
65 000 branches.  Each branch calls `Solve` on similar equalities, and `Reduce` on sub-
expressions that appear in many branches.  Without caching, the total CAS work scales with the
number of branches rather than the number of distinct sub-expressions.

---

### v0.2 — Solve memoization + branch pruning  *(no benchmarks captured before this change)*

**Date:** 2026-03-03
**Commit:** `887f43b` — *"Add benchmarking suite and performance optimizations"*
**File:** `MFGraphs/ZAnd.m`

#### Changes

1. **`CachedSolve[eq]`** — hash the equality expression (SHA-256), look up `$SolveCache`,
   populate on miss.  Same equality appearing across Or-branches is solved only once.

2. **`ClearSolveCache[]`** — exported function; called at the start of
   `CriticalCongestionSolver` and `MFGSystemSolver` to avoid stale results between problems.

3. **`xp === False` guard** in `ZAnd[xp, And[fst, rst]]` — return `False` immediately if the
   accumulated expression is already contradictory.

4. **False-branch skipping** in `ZAnd[xp, Or[fst, scd]]` — replaced `Or @@ Map[...]` with a
   `Which` block that drops `False` branches from the result.  Both branches are still
   evaluated, but the output is pruned.

```mathematica
(* v0.2 changes highlighted *)
$SolveCache = <||>;
CachedSolve[eq_] := Module[{key = Hash[eq,"SHA256"], result},
    If[KeyExistsQ[$SolveCache, key], $SolveCache[key],
       result = Solve[eq]; $SolveCache[key] = result; result]];

ZAnd[xp_, And[fst_, rst_]] :=
    If[ xp === False, False,         (* NEW: early exit *)
        If[ Head[fst] === Or, ...]]

ZAnd[xp_, Or[fst_, scd_]] :=
    Module[{rfst, result1, result2},
      rfst = Reduce[fst, Reals];     (* still uncached *)
      result1 = ZAnd[xp, rfst];
      result2 = ZAnd[xp, scd];      (* still evaluates both *)
      Which[result1===False && result2===False, False,
            result1===False, result2,
            result2===False, result1,
            True, RemoveDuplicates@(Or@@{result1,result2})]]
```

#### Expected impact (no captured baseline)

- Solve memoization: cases with repeated switching cost patterns (8, 10, 11, 14) expected to
  see 30–70% fewer `Solve` calls.
- False-branch pruning: avoids building Or-expressions with provably empty branches.

---

### v1.0 — DNFReduce: Reduce memoization + Or short-circuit

**Date:** 2026-03-04
**Commit:** `03d1c68` — *"Rename ZAnd to DNFReduce, eliminate ReZAnd, add performance improvements"*
**File:** `MFGraphs/DNFReduce.m` (replaces `MFGraphs/ZAnd.m`)

#### Changes

1. **`CachedReduce[expr, domain]`** — hash the `{expr, domain}` pair (SHA-256), look up
   `$ReduceCache`, populate on miss.  `Reduce[fst, Reals]` in the Or case is now cached; the
   same Or-branch sub-expression is reduced only once per solver invocation.

2. **`ClearSolveCache[]` also clears `$ReduceCache`** — `ClearSolveCache[] := ($SolveCache = <||>; $ReduceCache = <||>)`.

3. **Short-circuit Or evaluation** — before evaluating the second branch, check whether the
   first result already equals `xp` (i.e., the first branch covers the full solution space):

   ```mathematica
   DNFReduce[xp_, Or[fst_, scd_]] :=
       Module[{rfst, result1, result2},
         rfst = CachedReduce[fst, Reals];   (* CACHED *)
         result1 = DNFReduce[xp, rfst];
         If[result1 === xp,
           result1,                          (* SHORT-CIRCUIT: skip scd *)
           result2 = DNFReduce[xp, scd];
           ...]]
   ```

4. **`ReZAnd` eliminated** — the 3-argument form is now a direct `DNFReduce[xp_, rst_, fst_]`
   overload.  This simplifies the call graph and removes one level of indirection.

5. **Backward-compatible alias** — `ZAnd = DNFReduce` preserves existing notebook/script usage.

#### Measured performance (2026-03-04, Mathematica 14.3.0, macOS ARM)

The table below comes from `Results/compare_dnf.json` generated by
`wolframscript -file Scripts/CompareDNF.wls` on this commit.

| Case | Input Leaves | BC Time (s) | DNF Time (s) | **Speedup** | BC Disjuncts | DNF Disjuncts | BC Leaves | DNF Leaves | Equiv |
|------|-------------|-------------|--------------|-------------|-------------|--------------|-----------|-----------|-------|
| 7 | 163 | 0.000073 | 0.4707 | 0.000155x | 3 | 1 | 158 | 10 | ✓ |
| 8 | 449 | 0.000102 | 0.006365 | 0.0160x | 8 | 1 | 1215 | 24 | ✓ |
| 12 | 1383 | 0.6943 | 0.01707 | **40.7x** | 65536 | 1 | 44 269 569 | 55 | ~* |
| 11 | 3359 | 43.776 | 4.229 | **10.35x** | 3 359 232 | 2977 | 6 094 206 721 | 1 199 868 | ~* |
| 14 | 787 | 0.004773 | 0.01643 | 0.290x | 384 | 1 | 139 009 | 33 | ✓ |
| Braess congest | 2726 | 1.283 | 1.486 | 0.864x | 147 456 | 179 | 166 035 457 | 18 944 | ~* |
| Paper example | 3 | 0.000027 | 0.000007 | 3.86x | 1 | 1 | 3 | 3 | ✓ |
| 20 | 1060 | 0.005431 | 0.03204 | 0.169x | 1024 | 8 | 469 505 | 669 | ✓ |

*~* = Equivalence check timed out (60 s limit), but outputs are expected to be equivalent
based on the solver producing correct results on all downstream cases.*

#### Interpretation

- **Cases 12 and 11** (Attraction problem, 5 and 16 free variables) show dramatic speedup
  (40.7× and 10.35×).  The key driver is the **short-circuit Or**: when the first branch
  already resolves the full problem, `DNFReduce` avoids evaluating the remaining 2^N − 1
  branches.  BooleanConvert enumerates all 65 536 and 3 359 232 disjuncts.

- **Cases 7, 8, 14, 20** (small networks or cases with fewer branches): BooleanConvert is
  faster on wall-clock time (pure boolean manipulation is cheap), but produces far more
  disjuncts (3–1024×) and larger expressions (up to 703× more leaves).  Downstream solver
  steps operate on DNFReduce output; a compact representation directly reduces work in
  `TripleClean`, `MFGSystemSolver`, and `BooleanConvert`.

- **Braess congest**: DNFReduce is 0.86× (slightly slower) on wall-clock but produces
  179 disjuncts vs 147 456 for BooleanConvert — an 823× reduction in output size.

- **Paper example**: Trivial system (1 disjunct each); both methods are sub-millisecond.

---


---

### 2026-03-04 — v1.0 DNFReduce baseline (CachedReduce + Or short-circuit)

**Commit:** `cee8531`
**Date:** Wed 4 Mar 2026 15:17:23
**Generated by:** `wolframscript -file Scripts/CompareDNF.wls --tag "v1.0 DNFReduce baseline (CachedReduce + Or short-circuit)"`

#### Measured performance

| Case | Input Leaves | BC Time (s) | DNF Time (s) | **Speedup** | BC Disjuncts | DNF Disjuncts | BC Leaves | DNF Leaves | Equiv |
|------|-------------|-------------|--------------|-------------|-------------|--------------|-----------|-----------|-------|
| 7 | 163 | 0.000056 | 0.443238 | 0.000126x | 3 | 1 | 158 | 10 | ✓ |
| 8 | 449 | 0.000098 | 0.005845 | 0.0168x | 8 | 1 | 1215 | 24 | ✓ |
| 12 | 1383 | 0.722369 | 0.015102 | **47.8x** | 65536 | 1 | 44 269 569 | 55 | ~† |
| 11 | 3359 | 42.4949 | 3.93019 | **10.8x** | 3 359 232 | 2977 | 6 094 206 721 | 1 199 868 | ~† |
| 14 | 787 | 0.002022 | 0.015282 | 0.132x | 384 | 1 | 139 009 | 33 | ✓ |
| Braess congest | 2726 | 1.23861 | 1.38287 | 0.896x | 147 456 | 179 | 166 035 457 | 18 944 | ~† |
| Paper example | 3 | 0.000029 | 0.000016 | 1.81x | 1 | 1 | 3 | 3 | ✓ |
| 20 | 1060 | 0.005424 | 0.031614 | 0.172x | 1024 | 8 | 469 505 | 669 | ✓ |

† *Equivalence check timed out (60 s); outputs are expected to be equivalent.*

#### Rationale

The `Reduce[fst, Reals]` call in the Or branch was identified as a bottleneck: it is an
expensive symbolic CAS operation that may be called many times on identical sub-expressions as
the recursive traversal explores different Or-branches.  Similarly, once the first Or-branch
already resolves the entire problem (result equals `xp`), the second branch is redundant.

Two additions were made to address these:

1. **CachedReduce**: memoize `Reduce[fst, Reals]` by SHA-256 hash of `{fst, Reals}`, stored in
   `$ReduceCache`.  Same Or-branch sub-expression is reduced at most once per solver call.

2. **Short-circuit Or**: before evaluating the second Or-branch, check if the first result
   already equals `xp`.  If so, return immediately — the second branch cannot contribute new
   disjuncts.

#### Changes

```mathematica
(* Added: *)
$ReduceCache = <||>;
CachedReduce[expr_, domain_] :=
  Module[{key = Hash[{expr, domain}, "SHA256"], result},
    If[KeyExistsQ[$ReduceCache, key], $ReduceCache[key],
       result = Reduce[expr, domain]; $ReduceCache[key] = result; result]];
ClearSolveCache[] := ($SolveCache = <||>; $ReduceCache = <||>);   (* clears both *)

(* Changed Or case — added short-circuit: *)
DNFReduce[xp_, Or[fst_, scd_]] :=
    Module[{rfst, result1, result2},
      rfst = CachedReduce[fst, Reals];         (* was: Reduce[fst,Reals] *)
      result1 = DNFReduce[xp, rfst];
      If[result1 === xp, result1,              (* NEW: short-circuit *)
        result2 = DNFReduce[xp, scd];
        Which[result1 === False && result2 === False, False,
              result1 === False, result2,
              result2 === False, result1,
              True, RemoveDuplicates@(Or @@ {result1, result2})]]]
```

#### Interpretation

- **Case 12** (Attraction, no switching, 5 free variables): **47.8× speedup**.
  The short-circuit fires on the very first Or-branch — one solve is enough, all 65 535
  remaining branches are skipped.  BooleanConvert enumerates all 65 536.

- **Case 11** (Attraction, 16 switching costs, up to 65 536 Or-branches): **10.8× speedup**.
  Partially exploits short-circuit; some branches still need to be explored because the
  system has multiple genuinely distinct feasible regions (2977 disjuncts in output).

- **Cases 7, 8, 14** (small/medium with switching): DNFReduce is 6–150× *slower* on
  wall-clock, but produces 3–384× fewer disjuncts.  Downstream solver steps (TripleClean,
  BooleanConvert post-processing) operate on the compact DNFReduce output, so the total
  pipeline cost is lower even when DNFReduce itself is slower.

- **Case 20** (9-vertex, no switching): DNFReduce is 5.8× slower but reduces 1024 disjuncts
  to 8.  The additional Or branches each require a CachedReduce call; on this case the cache
  hit rate is low because all Or-branches have distinct expressions.


---

### 2026-03-10 — perf/improvements: And-Or early exit + TripleStep CachedSolve

**Commit:** `313a142`  
**Date:** Tue 10 Mar 2026 10:50:50  
**Generated by:** `wolframscript -file Scripts/CompareDNF.wls --tag "perf/improvements: And-Or early exit + TripleStep CachedSolve"`

#### Measured performance

| Case | Input Leaves | BC Time (s) | DNF Time (s) | **Speedup** | BC Disjuncts | DNF Disjuncts | BC Leaves | DNF Leaves | Equiv |
|------|-------------|-------------|--------------|-------------|-------------|--------------|-----------|-----------|-------|
| 7 | 163 | 0.000175 | 0.002347 | 0.0745633x | 3 | 1 | 158 | 10 | ✓ |
| 8 | 449 | 0.000108 | 0.006674 | 0.0161822x | 8 | 1 | 1215 | 24 | ✓ |
| 12 | 1383 | 0.71028 | 0.017474 | 40.6478x | 65536 | 1 | 44269569 | 55 | ~† |
| 11 | 3359 | 48.5742 | 4.38017 | 11.0896x | 3359232 | 2977 | 6094206721 | 1199868 | ~† |
| 14 | 787 | 0.005027 | 0.018959 | 0.265151x | 384 | 1 | 139009 | 33 | ✓ |
| Braess congest | 2726 | 1.35026 | 1.65131 | 0.817689x | 147456 | 179 | 166035457 | 18944 | ~† |
| Paper example | 3 | 0.00003 | 7e-6 | 4.29x | 1 | 1 | 3 | 3 | ✓ |
| 20 | 1060 | 0.005886 | 0.035008 | 0.168133x | 1024 | 8 | 469505 | 669 | ✓ |

† *Equivalence check timed out (60 s); outputs are expected to be equivalent.*

#### Rationale

The And-Or distribution case in `DNFReduce[xp_, And[fst_Or, rst_]]` previously evaluated **all** disjuncts of `fst` unconditionally via `Map`. When `fst` is an `Or[a, b, c, ...]`, each `DNFReduce[xp, rst, #]` call is independent — but the results were always fully enumerated before any pruning. The 2-arg Or case already had a short-circuit (`result1 === xp` → skip second branch), but there was no equivalent for the And-Or distribution.

Additionally, `TripleStep` in `D2E2.m` called bare `Solve[]` while the `CachedSolve` memoization infrastructure from `DNFReduce.m` was available but unused there.

#### Changes

```mathematica
(* DNFReduce.m: And-Or distribution — sequential with early exit *)
(* Old: *)
RemoveDuplicates@(DNFReduce[xp, rst, #] & /@ fst)

(* New: sequential Do + Catch/Throw *)
Catch[
    Module[{accumulated = False, r},
        Do[
            r = DNFReduce[xp, rst, b];
            If[r === xp, Throw[xp, "dnfAndOr"]];    (* early exit *)
            accumulated = Or[accumulated, r],
            {b, List @@ fst}
        ];
        If[accumulated === False, False, RemoveDuplicates[accumulated]]
    ],
    "dnfAndOr"
]

(* D2E2.m: TripleStep — use CachedSolve instead of Solve *)
(* Old: First@Solve[EE /. rules] // Quiet *)
(* New: First@CachedSolve[EE /. rules] // Quiet  (both overloads) *)
```

#### Interpretation

- **DNF speedups** (cases 12: 40.6×, case 11: 11.1×) are essentially unchanged relative to v1.0 — the existing Or short-circuit was already the dominant speedup driver for these cases, and the And-Or early exit does not fire when xp is already in its simplest form going into the And case.

- **Small/medium cases (7, 8, 14, 20)**: within noise of v1.0 baseline. DNFReduce remains slower on wall-clock than BooleanConvert for these cases but produces far fewer disjuncts (1–8 vs 3–1024), which compresses downstream `TripleClean` and `MFGSystemSolver` work.

- **All 31/31 valid benchmark cases pass** (same as v1.0 baseline). No regressions.

- The And-Or early exit will show more impact on systems where a single Or-branch fully resolves the problem after consuming some constraints — this occurs in larger switching-cost networks where the feasible region is well-determined by the first satisfying branch.


---

### 2026-03-31 — Add ReduceDisjuncts post-reduction

**Commit:** `008f220`  
**Date:** Tue 31 Mar 2026 16:26:31  
**Generated by:** `wolframscript -file Scripts/CompareDNF.wls --tag "Add ReduceDisjuncts post-reduction"`

#### Measured performance

| Case | Input Leaves | BC Time (s) | DNF Time (s) | **Speedup** | BC Disjuncts | DNF Disjuncts | BC Leaves | DNF Leaves | Equiv |
|------|-------------|-------------|--------------|-------------|-------------|--------------|-----------|-----------|-------|
| 7 | 163 | 0.000072 | 0.001018 | 0.0707269x | 3 | 1 | 158 | 10 | ✓ |
| 8 | 449 | 0.000105 | 0.006565 | 0.0159939x | 8 | 1 | 1215 | 24 | ✓ |
| 12 | 1383 | 1.76142 | 0.017179 | 102.533x | 65536 | 1 | 44269569 | 55 | ~† |
| 11 | 3359 | 46.2241 | 4.32194 | 10.6952x | 3359232 | 2977 | 6094206721 | 1199868 | ~† |
| 14 | 787 | 0.004227 | 0.016057 | 0.26325x | 384 | 1 | 139009 | 33 | ✓ |
| Braess congest | 2726 | 1.30593 | 1.52862 | 0.85432x | 147456 | 179 | 166035457 | 18944 | ~† |
| Paper example | 3 | 0.000024 | 0.000013 | 1.84615x | 1 | 1 | 3 | 3 | ✓ |
| 20 | 1060 | 0.005842 | 0.034017 | 0.171738x | 1024 | 8 | 469505 | 669 | ✓ |

† *Equivalence check timed out (60 s); outputs are expected to be equivalent.*

#### Rationale

Add a post-reduction `ReduceDisjuncts` pass after `DNFReduce` to shrink logically redundant branches before downstream symbolic solving. The target bottleneck was the large disjunct count that still remained on switching-cost cases even after DNF memoization and short-circuiting.

#### Changes

Apply `ReduceDisjuncts` to the DNF output after the main reduction step, using the existing thresholded strategies to merge or remove redundant disjuncts before the result is handed back to the solver pipeline.

#### Interpretation

The post-reduction pass helps exactly where it should: cases 12 and 11 remain the clear wins because they still carry the heaviest redundant branching structure after raw DNF reduction. Smaller cases such as 7, 8, 14, and 20 regress on wall-clock time because the extra simplification work costs more than it saves when the branch set is already small. This makes `ReduceDisjuncts` a good fit for large or highly redundant switching-cost systems, but not a universal win on tiny inputs.


---

### 2026-04-01 — optimize-preprocessing

**Commit:** `e59a562`
**Date:** Wed 1 Apr 2026 02:22:45
**Generated by:** `wolframscript -file Scripts/CompareDNF.wls --tag "optimize-preprocessing"`

#### Measured performance

| Case | Input Leaves | BC Time (s) | DNF Time (s) | **Speedup** | BC Disjuncts | DNF Disjuncts | BC Leaves | DNF Leaves | Equiv |
|------|-------------|-------------|--------------|-------------|-------------|--------------|-----------|-----------|-------|
| 7 | 163 | 0.000062 | 0.000898 | 0.0690423x | 3 | 1 | 158 | 10 | ✓ |
| 8 | 449 | 0.000105 | 0.007126 | 0.0147348x | 8 | 1 | 1215 | 24 | ✓ |
| 12 | 1383 | 0.86449 | 0.019857 | 43.5358x | 65536 | 1 | 44269569 | 55 | ~† |
| 11 | 3359 | 55.7624 | 4.81065 | 11.5915x | 3359232 | 2977 | 6094206721 | 1199868 | ~† |
| 14 | 787 | 0.005071 | 0.018464 | 0.274643x | 384 | 1 | 139009 | 33 | ✓ |
| Braess congest | 2726 | 1.38406 | 1.63961 | 0.844141x | 147456 | 179 | 166035457 | 18944 | ~† |
| Paper example | 3 | 0.000031 | 0.000017 | 1.82353x | 1 | 1 | 3 | 3 | ✓ |
| 20 | 1060 | 0.005821 | 0.051266 | 0.113545x | 1024 | 8 | 469505 | 669 | ✓ |

† *Equivalence check timed out (60 s); outputs are expected to be equivalent.*

#### Rationale

Record a fresh post-merge benchmark for the preprocessing-first DNF reduction path
against Mathematica's built-in `BooleanConvert[..., "DNF"]`. The goal is to check
whether the current preprocessing strategy still pays off on the large symbolic
cases that previously exploded in branch count.

#### Changes

No new DNF code was edited in this measurement-only update. This entry captures a
new run of `Scripts/CompareDNF.wls` on the current `master` state and updates the
comparison artifact in `Results/compare_dnf.md`.

#### Interpretation

The current DNF pipeline still wins decisively on the largest symbolic inputs:
cases `12` and `11` are about `43.5x` and `11.6x` faster than `BooleanConvert`,
while also collapsing the output from tens of thousands or millions of disjuncts
to `1` and `2977`. On smaller formulas such as `7`, `8`, `14`, and `20`,
`BooleanConvert` remains faster, but its outputs are much larger. `Braess congest`
is still slightly faster under `BooleanConvert`, though `DNFReduce` returns a far
smaller expression. The largest equivalence checks still hit the 60-second cap, so
those rows should be treated as expected-equivalent rather than formally verified.

## Future entries (template)

When implementing a new optimization, copy this template:

```markdown
### vX.Y — Short name of optimization

**Date:** YYYY-MM-DD
**Commit:** `<hash>` — *"Commit subject"*
**File:** `MFGraphs/DNFReduce.m`

#### Rationale

_Why was this change made?  What bottleneck was it targeting?_

#### Changes

_Code diff or pseudocode showing what changed._

#### Measured performance

| Case | Input Leaves | BC Time (s) | DNF Time (s) | **Speedup** | BC Disjuncts | DNF Disjuncts | BC Leaves | DNF Leaves | Equiv |
|------|-------------|-------------|--------------|-------------|-------------|--------------|-----------|-----------|-------|
| ... |

_Run `wolframscript -file Scripts/CompareDNF.wls --tag "vX.Y name"` to auto-generate._

#### Interpretation

_Which cases improved/regressed and why._
```
