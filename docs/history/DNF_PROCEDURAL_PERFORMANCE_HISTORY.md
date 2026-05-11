# dnfReduceProcedural Performance History

This document tracks optimization milestones applied to `dnfReduceProcedural`
(the stack-based iterative DNF reducer in `solversTools`Private``). Its
counterpart `dnfReduce` is logged in `DNF_PERFORMANCE_HISTORY.md`.

`dnfReduceProcedural` exists to bypass `$RecursionLimit` on deeply nested
Or-chains. It is structurally equivalent to `dnfReduce`: same output up to
ordering, as verified by `Results/dnf_recursive_vs_procedural.md`.

---

## How to add a new entry

Run the baseline profile and record the new TSV alongside the rationale:

```bash
caffeinate -dimsu wolframscript -file Scripts/ProfileDNFProcedural.wls \
    --case all --timeout 300
```

Then append a new section below following the template at the bottom.

---

## Metrics glossary

| Column | Meaning |
|--------|---------|
| **Vars** | Variables in `buildSolverInputs[sys]` |
| **Conj** | Top-level conjuncts in the system constraint |
| **BuildMs** | `makeSystem` wall-clock (ms) — context only, not part of the reducer |
| **ProcMs** | `dnfReduceProcedural` wall-clock (ms) |
| **Branches** | Or-disjunct count in the output (1 if a single And, 0 if `True`/`False`) |
| **Leaves** | `LeafCount` of the output expression |

---

## Entries (oldest → newest)

---

### 2026-05-11 — Baseline (introduction)

**Commit:** `4648fac` — *"docs(solversTools): annotate recursive functions and add iterative alternative"*
**Hardware:** Apple Silicon, Mathematica 14.x

#### What `dnfReduceProcedural` does

Stack-based equivalent of `dnfReduce`. The recursive version recurses on
`Or` and `Equal` heads; the procedural version maintains an explicit
worklist of `<|xp, rst, fst|>` items. Or-branches are pushed with
`Prepend + Reverse` so branch priority matches the recursive ordering.
The `And`-case loop folds in `Equal`-induced substitutions in place
without recursion, then queues remaining unsolved disjuncts.

Not wired into `dnfReduceSystem` by default; called directly from
`CompareSingleCase.wls` and a few diagnostic paths.

#### Measured performance

| Case | Vars | Conj | BuildMs | ProcMs | Status | Branches | Leaves |
|---|---|---|---|---|---|---|---|
| grid-2x2 | 4  | 36  | 13.66 | 2.63      | OK | 1     | 33         |
| grid-3x2 | 9  | 63  | 10.79 | 1.92      | OK | 1     | 69         |
| grid-2x3 | 9  | 63  | 8.59  | 2.34      | OK | 1     | 69         |
| grid-2x4 | 14 | 90  | 9.23  | 3.79      | OK | 1     | 113        |
| grid-3x3 | 21 | 114 | 16.44 | 12.57     | OK | 10    | 1 661      |
| grid-3x4 | 33 | 165 | 17.17 | 50.92     | OK | 30    | 9 411      |
| grid-4x4 | 52 | 240 | 48.17 | 11 306.31 | OK | 7 200 | 3 230 861  |
| grid-4x5 | 71 | 315 | 28.96 | 89 304.62 | OK | 27 000 | 19 251 001 |

#### Cross-reference vs recursive (`dnfReduce`)

Both columns below are taken from a single same-run report
(`Results/dnf_recursive_vs_procedural.md`) so the ratio is apples-to-apples.
Canonicalized verdict was `OK` on every row — same branch sets, just
ordering differs. The `ProcMs` column in the baseline table above is from
a separate fresh run today and will differ on small cases by JIT/cache
noise.

| Case | Recursive (s) | Procedural (s) | Ratio P/R |
|---|---|---|---|
| grid-2x2 | 0.002 | 0.001 | 0.5 |
| grid-3x2 | 0.002 | 0.002 | 1.0 |
| grid-2x3 | 0.002 | 0.002 | 1.0 |
| grid-2x4 | 0.004 | 0.003 | 0.75 |
| grid-3x3 | 0.014 | 0.013 | 0.93 |
| grid-3x4 | 0.054 | 0.049 | 0.91 |
| grid-4x4 | 9.72  | 10.83 | 1.11 |
| grid-4x5 | 57.39 | 79.89 | 1.39 |

#### Interpretation

- **Procedural matches recursive on disjunct count and structure** on every
  case tested. No correctness gap.
- **Procedural is at parity or slightly faster on small cases**
  (grid-2x2 through grid-3x4) — likely because skipping the WL
  function-call dispatch for shallow inputs offsets the worklist
  bookkeeping cost.
- **Wall-clock cost grows faster than the recursive version at scale.**
  At grid-4x4 procedural is ~1.11× slower; at grid-4x5 ~1.39× slower.
  The crossover is around grid-4x4 (52 vars, 7 200 output branches).
- **Likely hot spot at scale**: the `AppendTo[results, ...]` and `stack =
  Prepend[...]` / `stack = Rest[stack]` allocations are O(n) copies each.
  At 27 000 leaf branches on grid-4x5 the worklist churn dominates.
  The recursive engine avoids this by leaning on the WL function-call stack.
- **Output leaf count** scales roughly linearly with branch count (~700
  leaves/branch on grid-4x4 / 4x5). Suggests each branch carries the full
  residual rather than being compacted.

#### Open optimization candidates (not yet attempted)

1. Replace `stack` (List with Prepend/Rest) with an internally mutable
   `CreateDataStructure["Stack"]` or a packed array with an index cursor.
2. Replace `AppendTo[results, ...]` with a `Sow`/`Reap` collector or a
   pre-allocated buffer.
3. Coalesce identical residuals before stacking (cheap branch-dedup —
   recursive uses `RemoveDuplicates` after merge; procedural doesn't).
4. Restore a `CachedSolve`-style memo around `Quiet[Solve[elem], ...]`.
   Both recursive (`solversTools.wl:1348`, `1384`) and procedural
   (`solversTools.wl:1440`, `1461`) currently call bare `Solve` per
   `Equal`. The historical `$SolveCache` / `$ReduceCache` from
   `DNF_PERFORMANCE_HISTORY.md` v0.2 / v1.0 lived in the old
   `DNFReduce.m` and did not survive the move to `solversTools`. A memo
   shared by both engines would help wherever Or-branches share equality
   atoms.

---

### 2026-05-11 — Add `cachedSolve` memo (shared with recursive engine)

**Commit:** _pending_
**File:** `MFGraphs/solversTools.wl`

#### Rationale

The baseline entry (above) flagged "restore a `CachedSolve`-style memo"
as a candidate. Instrumenting the six `Quiet[Solve[expr], Solve::svars]`
call sites (L1348, L1384, L1440, L1461, L1626, L1654) and counting
hashed-by-content keys per `dnfReduceSystem` invocation gave:

| Case     | Solve calls | Distinct | Hit rate |
|----------|-------------|----------|----------|
| grid-2x2 | 6           | 6        | 0%       |
| grid-3x2 | 16          | 16       | 0%       |
| grid-3x3 | 117         | 49       | 58%      |
| grid-3x4 | 434         | 80       | 80%      |
| grid-4x4 | **85 103**  | **140**  | **99.8%** |

The "post-substitution rewrites the Equal atoms" worry in the baseline
turned out to be empirically wrong — Or-branches reconverge on the same
constraints in massive volume on the harder cases. Original plan
proposed a shape-gated cache (trivial vs non-trivial Equals); measurement
showed trivial and non-trivial Solve costs are nearly identical at scale
(~43 vs ~46 µs), so the gate was dropped — `cachedSolve` memoizes all
Solve calls inside the DNF reducers.

#### Changes

- Added top-level state in `solversTools`Private``: `$dnfSolveCache`,
  `$dnfSolveHits`, `$dnfSolveMiss`. Cache keyed by `Hash[expr]` (default
  32-bit hash, much cheaper than the historical SHA-256).
- Counters reset per `dnfReduceSystem` invocation; cache **persists
  across calls** (content-hashing makes this correctness-safe and
  long-session cache sizes stay small — 267 entries after running all
  five grid sizes consecutively).
- Internal `clearSolveCache[]` helper available in `solversTools`Private``
  for forced eviction (not exposed publicly — cache stays small in
  practice).
- Six `Quiet[Solve[...], Solve::svars]` call sites replaced with
  `cachedSolve[...]`.

#### Measured performance

`Scripts/ProfileDNFProcedural.wls --case all --timeout 300`, same hardware
as the baseline:

| Case | Vars | Conj | BuildMs | ProcMs | Δ vs baseline | Branches | Leaves |
|---|---|---|---|---|---|---|---|
| grid-2x2 | 4  | 36  | 12.79 | 2.02      | noise  | 1     | 33         |
| grid-3x2 | 9  | 63  | 10.75 | 2.18      | noise  | 1     | 69         |
| grid-2x3 | 9  | 63  | 8.39  | 2.00      | noise  | 1     | 69         |
| grid-2x4 | 14 | 90  | 9.25  | 3.59      | noise  | 1     | 113        |
| grid-3x3 | 21 | 114 | 15.70 | 10.60     | 1.19×  | 10    | 1 661      |
| grid-3x4 | 33 | 165 | 16.79 | 36.59     | 1.39×  | 30    | 9 411      |
| grid-4x4 | 52 | 240 | 45.31 | 7 987.88  | **1.42×** | 7 200 | 3 230 861 |
| grid-4x5 | 71 | 315 | 27.49 | 68 563.01 | **1.30×** | 27 000 | 19 251 001 |

#### Interpretation

- Procedural benefits *more* than recursive in absolute terms because it
  was the slower engine on these cases — the Solve fraction of total
  time was larger there.
- Cache only kicks in on cases where the Or-branch fan-out exceeds the
  number of distinct equations; below grid-3x3 the noise floor swamps
  the savings.
- The new state lives in `solversTools`Private`` and is shared by
  `dnfReduce`, `dnfReduceProcedural`, and `dnfReduceInstrumented`.
  Other solvers (`booleanMinimizeSystem`, `findInstanceSystem`, etc.)
  do not call Solve inside the DNF reducers, so they are unaffected.

#### Verification

- `Scripts/RunTests.wls fast` — 216/216 pass.
- `Scripts/RunSingleTest.wls MFGraphs/Tests/reduce-system.mt` — 88/88
  pass with cache active.

---

## Future entries (template)

```markdown
### YYYY-MM-DD — short name

**Commit:** `<hash>` — *"subject"*

#### Rationale
_Bottleneck targeted._

#### Changes
_Diff or pseudocode._

#### Measured performance

| Case | Vars | Conj | BuildMs | ProcMs | Status | Branches | Leaves |
|---|---|---|---|---|---|---|---|
| ... |

#### Interpretation
_Which cases improved / regressed and why._
```
