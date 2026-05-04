# MFGraphs Issue Fix Implementation Plan

**Author:** Manus AI  
**Date:** 2026-05-03  
**Scope:** Execution plan for the prioritized repository issues in `docs/planning/PRIORITIZED_ISSUE_LIST.md`

## Executive summary

This plan is designed to fix the current MFGraphs issues in a way that is **incremental, reversible, benchmark-aware, and low-risk**. The repository is already in a strong architectural position: it has a staged kernel centered on typed scenarios, symbolic unknown generation, system construction, and symbolic solving.[1] The immediate need is not to expand capability breadth, but to **stabilize the public workflow**, **protect it with tests**, and **reduce symbolic execution cost** in the most concentrated hotspots.[1] [2] [3]

The central principle of this plan is that every milestone should leave the repository in a **cleaner and more defensible state than before**, even if later milestones are postponed. That means the work is intentionally split into short vertical slices rather than one large refactor.

| Milestone | Goal | Risk level | User input needed |
|---|---|---|---|
| A | Stabilize the public scenario-first workflow | Low | Yes, one policy decision |
| B | Protect the workflow with tests and benchmarks | Low-to-moderate | No, unless benchmark thresholds should be stricter |
| C | Refactor symbolic preprocessing and block structure | Moderate | No |
| D | Simplify data ownership and exactification policy | Moderate | Yes, if public options are exposed |
| E | Documentation and repo-ops hardening | Low | No |

## Operating principles

The plan assumes the repository should continue to honor the current core scenario-kernel phase. That means the active surface remains centered on `makeScenario`, `getExampleScenario`, `makeSymbolicUnknowns`, `makeSystem`, `solveScenario`, and the active solver/tooling stack described in the repository guidance.[1]

The second principle is that **performance work must stay measurement-backed**. The existing benchmark history already shows topology-sensitive behavior in `dnfReduceSystem`, including a roughly 4.4× difference between `grid-2x3` and `grid-3x2`, and raw `reduceSystem` timeout behavior on `example-12` even at 160 seconds.[3] That means any solver or system-structure change must be evaluated against a fixed representative case set before and after the edit.

The third principle is that **reversibility matters**. Each milestone should be implementable as a small branch or small PR sequence so that if a change introduces ambiguity or benchmark regression, the repository can roll back the last isolated step rather than unwind a mixed batch.

## Milestone A — Stabilize the public scenario-first workflow

This is the highest-leverage milestone because it aligns the public API with the repository’s actual architecture.

### A1. Confirm and freeze the public workflow target

| Task | Detail |
|---|---|
| Goal | Freeze the canonical user-facing workflow for the current phase |
| Target workflow | `Needs["MFGraphs`"]` → `makeScenario` → `makeSymbolicUnknowns` → `makeSystem` → `solveScenario` and `SolveMFG` as the compatibility-facing top-level wrapper |
| Why first | It prevents implementation work from drifting while policy is still implicit |
| Output | A short policy note or section in the planning docs stating which entrypoints are canonical, compatibility-only, and inactive |

**Implementation notes:**

The code and guidance already imply this direction, but it should be made explicit before touching entrypoint routing. The key clarification is whether `SolveMFG` is intended to become a fully supported scenario-native top-level entrypoint or remain only a thin compatibility alias around the lower-level staged path.[1] [2]

### A2. Implement `SolveMFG[s_scenario, opts___]`

| Task | Detail |
|---|---|
| Goal | Add direct scenario dispatch to the top-level compatibility solver entrypoint |
| Main change | Add a typed scenario overload that delegates through the same canonical build-and-solve path as `solveScenario` |
| Files likely touched | `MFGraphs/orchestrationTools.wl`, tests, and documentation |
| Risk | Low, if implemented as delegation rather than duplicated logic |

**Detailed approach:**

The implementation should be intentionally thin. The new overload should not create a parallel scenario-only solving stack. It should normalize options and route through the same underlying path already used by `solveScenario`, so there is only one actual execution model to maintain.

> Preferred pattern: add the new public dispatch, keep behavior centralized, and avoid copying preprocessing or solver-selection logic into multiple wrappers.

**Validation for A2:**

| Validation type | Check |
|---|---|
| Unit | Scenario input returns the same result shape as `solveScenario[s]` on a representative small case |
| Compatibility | Existing association-based `SolveMFG` calls still work |
| Structural | No duplicate or shadowed dispatch definitions are introduced |
| Documentation | Usage text and quick-start examples remain consistent |

**Rollback point:**

If the new dispatch causes ambiguity or breaks compatibility behavior, revert only the dispatch addition and its related docs/tests. Do not bundle it with performance work.

### A3. Resolve the example-surface policy

| Task | Detail |
|---|---|
| Goal | Decide which example helper is canonical in the active phase |
| Current ambiguity | The repository has a legacy example path, a transitional typed constructor path, and a scenario-first example factory path in planning/review materials.[2] |
| Why this milestone | Tests and docs cannot be made stable until this policy is explicit |
| Output | One policy table documenting canonical, compatibility-only, and deprecated/example-fixture surfaces |

**Decision choices:**

| Choice area | Recommended default |
|---|---|
| Quick-start examples | Prefer `getExampleScenario` |
| Explicit transitional construction | Keep only if it still serves active tests or migration needs |
| Legacy raw example data | Keep compatibility-only or archive if active code no longer depends on it |

**User decision needed:** yes. I recommend that you confirm whether you want the repository to present **`getExampleScenario` as the canonical user-facing example entrypoint** for the current phase.

## Milestone B — Protect the workflow with tests and benchmark gates

Once the public story is frozen, the next step is to make regression visible immediately.

### B1. Add scenario-first routing tests

| Task | Detail |
|---|---|
| Goal | Ensure the newly canonical path is covered in the active test suite |
| Main checks | `SolveMFG[s_scenario]`, `solveScenario[s]`, and package-loading expectations for the chosen active public symbols |
| Files likely touched | `MFGraphs/Tests/orchestration.mt`, `MFGraphs/Tests/package-loading.mt`, possibly `scenario-kernel.mt` |
| Risk | Low |

**Detailed test additions:**

1. A small scenario fixture should be solved via `solveScenario[s]` and `SolveMFG[s]`, with result kinds compared.
2. A direct dispatch test should verify that scenario input is accepted without conversion by the caller.
3. Package-loading tests should check only the symbols that are truly part of the active public contract after the example-policy decision.

### B2. Add direct benchmark regression gates

| Task | Detail |
|---|---|
| Goal | Turn existing benchmark evidence into a repeatable comparison baseline |
| Main scripts | `Scripts/BenchmarkSystemSolver.wls`, `Scripts/ProfileScenarioKernel.wls`, and `BENCHMARKS.md` |
| Benchmark set | Keep it small and stable |
| Risk | Low-to-moderate |

**Recommended benchmark fixture set:**

| Case class | Suggested case | Purpose |
|---|---|---|
| Tiny | `chain-3v-1exit` | Sanity and fast iteration |
| Small branched | `chain-3v-2exit` or `example-7` | Branching behavior |
| Medium topology-sensitive | `grid-3x2` and `grid-2x3` | Topology sensitivity |
| Hard representative | `example-12` | Timeout/solver stress |

**Recommended metrics:**

| Metric | Why keep it |
|---|---|
| `BuildMs` | Detects system-construction regressions |
| `WarmupMs` and `RepMs` | Detects solve regressions |
| `Status` | Makes timeouts explicit |
| `Kind` | Detects shifts in residual/branched/rules output shape |
| `Valid` | Protects correctness |
| Byte count or structural counts | Helps interpret memory-oriented regressions |

**Acceptance rule:**

Do not enforce a hard timing ceiling at first. Instead, require that each solver-sensitive change reports before/after results on the fixed case set and flags any increase above an agreed review threshold, such as 20–25% on non-timeout cases.

**User input needed:** optional. If you already have preferred regression thresholds or preferred representative examples, I can align the plan to them; otherwise I recommend starting with relative review thresholds rather than hard fail thresholds.

### B3. Make staged profiling part of solver-change review

This should be treated as a lightweight operational rule rather than a heavy new system. For every solver-sensitive change, run the staged profiler on one easy case and one hard case and record where the cost moved: scenario build, unknown generation, system build, preprocessing, or solve.[3]

## Milestone C — Refactor the main symbolic execution hotspots

This is the highest technical-risk milestone, so it should happen only after Milestones A and B are complete.

### C1. Introduce a structured preprocessing path for `reduceSystem`

| Task | Detail |
|---|---|
| Goal | Stop repeatedly rebuilding and rescanning a monolithic Boolean constraint during preprocessing |
| Main change | Introduce an internal block-oriented preprocessing path, such as `reduceSystemBlocks` or an equivalent internal helper |
| Why it matters | The review identified repeated whole-expression scans and repeated substitution into a giant base constraint as the main time hotspot.[2] |
| Risk | Moderate |

**Detailed execution steps:**

1. Preserve the current public `reduceSystem` signature.
2. Internally split constraints into deterministic linear equations, inequalities, complementarity blocks, and residual nonlinear/disjunctive blocks.
3. Perform linear elimination against structured blocks rather than against a repeatedly rebuilt `And[...]` expression.
4. Reassemble only the residual symbolic problem that truly requires `Reduce` or DNF-style handling.

**Validation strategy:**

| Check type | Requirement |
|---|---|
| Correctness | Existing `reduce-system.mt` tests still pass |
| Result-shape stability | `solutionResultKind` is unchanged on the core benchmark set unless an intentional improvement is documented |
| Performance | `grid-3x2`, `grid-2x3`, and `example-12` do not regress; at least one representative case improves measurably |
| Maintainability | New internal helpers are documented clearly enough for future debugging |

**Rollback point:**

Keep the old preprocessing helper reachable behind a temporary internal switch until the new path proves stable on tests and sample benchmarks.

### C2. Delay flattening of complementarity-heavy blocks

| Task | Detail |
|---|---|
| Goal | Reduce symbolic inflation before the solver stage |
| Main change | Store block data as grouped lists or typed subrecords longer, and flatten only at the actual solver boundary |
| Risk | Moderate |

**Detailed execution steps:**

1. Audit `buildFlowData` and `buildComplementarityData` for eager `And@@` construction.
2. Preserve the human-readable grouping structure in the system object.
3. Add a narrow conversion helper that produces the flattened Boolean form only when a solver requests it.
4. Reuse that helper consistently across solver entrypoints so flattening happens in one place.

**Expected benefit:**

This should reduce `ByteCount` growth in `mfgSystem` and make block-level profiling easier to interpret.[2]

### C3. Re-benchmark after each hotspot change

Do **not** batch C1 and C2 together without measurement. After each sub-change, run the small benchmark set and record the result next to the branch or PR. This keeps causality clear.

## Milestone D — Simplify data ownership and exactification policy

This milestone is important, but it should come only after the main symbolic hotspots are under better control.

### D1. Reduce duplicated payload in `mfgSystem`

| Task | Detail |
|---|---|
| Goal | Stop the system object from behaving like a warehouse for already-owned upstream data |
| Main change | Define canonical ownership boundaries between `scenario`, unknowns, and `mfgSystem` |
| Risk | Moderate |

**Detailed execution steps:**

1. Inventory which top-level `mfgSystem` fields duplicate scenario topology and unknown-family data.
2. Mark each field as one of: essential local system field, convenience compatibility field, or redundant.
3. Remove or de-emphasize redundant fields in internal storage.
4. Retain accessor compatibility through typed accessors or flattening helpers where needed.

**Validation:** compare representative `ByteCount` before and after, and ensure plotting and solver utilities still function.

### D2. Make exactification configurable or delayed

| Task | Detail |
|---|---|
| Goal | Avoid paying the full symbolic exactness cost when it is not necessary |
| Main change | Introduce an option or internal control point for exactification timing |
| Risk | Moderate, because symbolic semantics matter |

**Recommended path:**

Start internally rather than exposing a public user option immediately. First prove that delayed or selective exactification helps on benchmark cases without changing correctness expectations. Only then decide whether it should become public.

**User input needed:** only if you want this to become a documented public option. Otherwise I recommend keeping it internal at first.

## Milestone E — Documentation and operational hardening

This milestone should land after the behavior is stable enough that documentation will not churn immediately again.

### E1. Harden guard scripts and path portability

| Task | Detail |
|---|---|
| Goal | Make policy-enforcement scripts portable and reliable |
| Main change | Remove machine-specific path assumptions and fold the checks into normal repo workflows |
| Risk | Low |

### E2. Consolidate active-phase architecture documentation

| Task | Detail |
|---|---|
| Goal | Publish one concise source of truth for the current active architecture |
| Main change | Add a current-architecture document and link it from README and planning docs |
| Risk | Low |

## Proposed PR structure

To keep the work robust, I recommend using a small-PR sequence rather than one broad implementation branch.

| PR | Scope | Contents |
|---|---|---|
| PR-1 | Public workflow | `SolveMFG[s_scenario]`, docs touch-up, minimal direct tests |
| PR-2 | Example policy and package-surface cleanup | Example policy docs, package-loading tests, any symbol-surface adjustments |
| PR-3 | Benchmark gate setup | Benchmark fixture set, reporting convention, optional lightweight guard script integration |
| PR-4 | Preprocessing refactor part 1 | Structured preprocessing path behind stable public API |
| PR-5 | Complementarity/block-structure refactor | Deferred flattening and block-boundary cleanup |
| PR-6 | Data ownership cleanup | `mfgSystem` payload trimming and compatibility access refinements |
| PR-7 | Exactification follow-up | Internal delayed/configurable exactification experiment and measured conclusion |
| PR-8 | Docs and ops hardening | Current architecture doc, guard-script portability, contributor guidance |

## Validation checklist per PR

Every PR in this plan should satisfy the same standard checklist.

| Category | Required check |
|---|---|
| Correctness | Run the active fast suite |
| API stability | Confirm the active public contract affected by the PR |
| Benchmarking | Run the small representative benchmark set if solver or system behavior changed |
| Profiling | Run staged profiling for solver-sensitive changes |
| Documentation | Update any usage or architecture note changed by the PR |
| Rollback safety | Keep changes scoped so they can be reverted independently |

## What input I need from you

I **do not need more input** to draft or begin the technical execution plan itself. I can proceed with the implementation sequence using the current repository guidance.

However, there are **three decision points** where your confirmation would improve stability and avoid policy churn.

| Decision | Needed? | My recommendation |
|---|---|---|
| Should `getExampleScenario` be the canonical example entrypoint for the active phase? | **Yes** | **Yes**, make it canonical for quick starts. |
| Should exactification remain an internal behavior change first, rather than a public option? | Optional | **Yes**, keep it internal first. |
| Do you want benchmark regression enforcement to begin as reporting-only or with hard thresholds? | Optional | Start as **reporting-only with review thresholds**. |

If you do not want to spend time deciding those now, I can still proceed under the recommended defaults above.

## Immediate next action

The best next move is to start **Milestone A / Task A2**: implement `SolveMFG[s_scenario, opts___]` in the thinnest possible way, then add direct scenario-routing tests immediately afterward. That gives the repo the most architectural value for the least risk.[1] [2]

## References

[1]: /mnt/desktop/MFGraphs/CLAUDE.md "MFGraphs canonical workflow guidance"
[2]: /home/ubuntu/mfgraphs_simplicity_performance_review.md "MFGraphs simplicity and performance review"
[3]: /mnt/desktop/MFGraphs/BENCHMARKS.md "MFGraphs benchmark policy and history"
[4]: /mnt/desktop/MFGraphs/docs/planning/PRIORITIZED_ISSUE_LIST.md "MFGraphs prioritized issue list"
