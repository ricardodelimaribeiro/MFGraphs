# MFGraphs Master Execution Plan

**Status**: Active — authoritative over all other planning docs in this directory.
**Last updated**: 2026-04-22

---

## Program Goals and Success Metrics

The overriding goal is to move MFGraphs from a research-grade codebase with an accidental public surface toward a **typed, scenario-first package** with an intentional API, contract-grade documentation, and validation proportional to risk.

| Goal | Measurable success criterion |
|---|---|
| Typed scenario kernel in production | `makeScenario → validateScenario → completeScenario → scenario[...]` pipeline works end-to-end on the 34 built-in examples |
| Scenario-first entrypoints | `DataToEquations` and `SolveMFG` accept `scenario[...]` directly without requiring external symbol juggling |
| Intentional public surface | Exported symbol count is reduced from 104 to a deliberate keep/advanced/compatibility set with no accidental exports (`alpha$`, `j$` fixed) |
| Paper readiness maintained | HRF Scenario 1 continues to produce `Feasibility → Feasible` on the critical solver within 150s |
| Test coverage on surviving API | Every "core" and "advanced-retained" symbol has at least one direct MUnit test reference |
| No accidental regression | Fast test suite (11 files per RunTests.wls) passes at zero failures after each phase |

---

## Current-State Baseline (2026-04-22)

### Branch and solver surface

- Active development branch: `codex/reduce-test-redundancy` (critical-only)
- `NonLinearSolver.wl`, `Monotone.wl`, and `TimeDependentSolver.wl` do **not** exist on this branch (deleted by PR #107).
- `rollback/pre-critical-only-surface` (at `c3711fa`) is a preserved local branch with the full solver suite — not the active branch.
- `master` is pinned at `5509995` (critical-only surface). **Do not force-push master.**

### Package structure (post-hygiene-fix refactor, 2026-04-22)

| File | Role |
|---|---|
| `MFGraphs/MFGraphs.wl` | Bootstrap: ::usage declarations, single Get block, shared private helpers, `MakeSolverResult` |
| `MFGraphs/FictitiousPlayBackend.wl` | Phase 5 Fictitious Play numeric backend — internal only |
| `MFGraphs/SolveMFGDispatch.wl` | Unified `SolveMFG` entrypoint, `Options[SolveMFG]`, cascade logic |
| `MFGraphs/Solvers.wl` | Critical-congestion solver suite |
| `MFGraphs/DataToEquations.wl` | Network topology → equation compiler |
| `MFGraphs/DNFReduce.wl` | Boolean algebra / disjunctive normal form solver |
| `MFGraphs/Graphics.wl` | Public visualization helpers |
| `MFGraphs/Scenario.wl` | Typed scenario kernel (`makeScenario`, `validateScenario`, `completeScenario`, `scenarioQ`) |

### Exported API surface

Current exported surface: **104 symbols** (see `public_api_symbol_inventory.md` for full table).

Known hygiene issues — status after Phase 2 (2026-04-22):
- `alpha$` and `j$` — neither exists in the current codebase (confirmed by grep; gate already satisfied)
- Duplicate `::usage` declarations — **resolved**: removed from `Solvers.wl`, `DataToEquations.wl`, `SolveMFGDispatch.wl`; canonical declarations now only in MFGraphs.wl BeginPackage
- Many "advanced" symbols exported without tests (60/104 have zero direct test references) — deferred to Phase 6/7

### Validation evidence

- Small tier (cases 1–6, 27): 21/21 passed across all three solvers (2026-04-14).
- Paper tier (HRF Scenario 1): CriticalCongestion = Feasible / ~122s (2026-04-14).
- Scenario-onward test subset: 74/74 passed (2026-04-22, post-organization refactor).

---

## Resolved Cross-Doc Conflicts

The following conflicts existed across the input planning documents. These decisions are now canonical.

| Conflict | Decision |
|---|---|
| **Solver scope for current release line** | Critical-congestion solver is the only public-API-grade solver for the current release. NonLinear, Monotone, and TimeDependentSolver remain in the codebase on `rollback/pre-critical-only-surface` and are tested but not part of the primary user-facing narrative. They may be re-exposed as `Method` options under `SolveMFG` once scenario routing is stable. |
| **Scenario-first vs. document-current-symbols-first** | Scenario-first wins. Do not freeze the current 104-symbol surface as the final API. Classify symbols relative to their role in the scenario-first workflow, then document only what survives the classification. (`public_api_full_surface_plan.md` is authoritative on this.) |
| **Paper claims** | Only validated artifacts may be cited as paper-ready. HRF Scenario 1 on CriticalCongestion is validated. NonLinear/Monotone results for that case are validated as of 2026-04-14. No new paper-readiness claims without re-running validation on the current branch. |
| **Fictitious Play backend exposure** | Internal only until the four documented vulnerabilities (GitHub #72–#75) are resolved. Do not wire it into `CriticalCongestionSolver` or expose it in `SolveMFG` until those issues are closed. |
| **MFGraphs.wl as bootstrap vs. as definition file** | Bootstrap only (resolved by 2026-04-22 refactor). Function definitions belong in submodules; MFGraphs.wl only holds package flags, small public helpers, shared private envelope builders, and `Get` calls. |

---

## Workstreams and Dependencies

```
Stream A: API Hygiene & Symbol Policy
  → blocks Stream C (documentation cannot be written until the surface is frozen)
  → blocks Stream D (tests should target the post-policy surface)

Stream B: Scenario-First Architecture
  → blocks Stream A (symbol classification depends on scenario workflow decisions)
  → blocks Stream C
  → blocks Stream D

Stream C: Documentation
  → requires Stream B and Stream A to be at least at Phase 3 gate

Stream D: Verification & Maintenance Gates
  → requires Stream A Phase 2 (symbol transition matrix exists)
  → can proceed in parallel with Stream C for new symbols in B

Stream E: Paper Readiness
  → independent of A/B/C/D but must not break validation baseline
  → re-validation required whenever solver or example changes land
```

---

## Phase-by-Phase Implementation Sequence

### Phase 1 — Scenario workflow specification ✅ COMPLETED 2026-04-22

**Gate**: written artifact exists, no code changes.

Deliverables:
- Canonical end-to-end user workflow — **done**: `makeScenario → validateScenario → completeScenario → SolveMFG → plot`
- Scenario transition matrix — **done**: all 104 symbols classified (14 core, 41 advanced-retained, 5 compatibility, 9 internal-blocked, 5 orphaned, 22 examples-only, 2 retire); see `phase1_scenario_transition_matrix.md`
- Quick-start outline for future README rewrite — **done**: in `phase1_scenario_transition_matrix.md`
- Validation boundary policy — **done**: structural errors at scenario boundary; infeasibility at solver result

Source inputs: `mfgraphs_scenario_schema_proposal.md`, `revised_mfgraphs_plan_response.md`, `public_api_full_surface_plan.md` (Phase 1 section), `issue_86_scenario_kernel.md`.

Acceptance gate: scenario transition matrix covers all 104 symbols with no blank decisions. ✅ Satisfied.

---

### Phase 2 — Fix export hygiene and duplicate usage declarations ✅ COMPLETED 2026-04-22

**Gate**: clean package load with zero shadow warnings on the core symbol set.

Deliverables:
- ~~Remove accidental exports `alpha$` and `j$`~~ — neither existed in current codebase; gate already satisfied
- Resolve duplicate `::usage` declarations — **done**: removed from `Solvers.wl`, `DataToEquations.wl`, `SolveMFGDispatch.wl`
- Fix `MFGraphs.wl` structure: single Get block (correct files only), removed duplicate function definitions (`EncodeFlowAssociation`, `DecodeFlowVector`, `ComputeSignedEdgeFlowsFast`, `ComputeKirchhoffResidualFast`, `SolveMFGCompiledInputQ`)
- Fixed `SolveMFGDispatch.wl` `Options[SolveMFG]` — was referencing non-existent `NonLinearSolver`/`MonotoneSolverFromData`

Source inputs: `public_api_gap_report.md` (Suspicious exports section).

Acceptance gate: `Names["MFGraphs`*"]` does not include `alpha$` or `j$`; no `::shdw` warnings when loading via `wolframscript -file` (not `-code` — inline invocation produces false positives). Fast suite: 0 failures.

---

### Phase 3 — Adapt core entrypoints to accept scenarios

**Gate**: `DataToEquations[s_scenario]` and `SolveMFG[s_scenario]` both work end-to-end.

Deliverables:
- `DataToEquations` accepts a `scenario[assoc]` object: reads `"Model"` and `"Data"` blocks, produces compiled equation association
- `SolveMFG` accepts a `scenario[assoc]` object as primary input while preserving backward compatibility with raw associations
- `GetExampleData` results can be wrapped with `makeScenario` without manual symbol juggling
- Scenario validation failures surface before solver execution (not deep inside symbolic reduction)
- MUnit tests cover the new scenario-input paths for `DataToEquations` and `SolveMFG`
- Fast suite still passes at zero failures

Source inputs: `public_api_full_surface_plan.md` (Phase 3 section), `revised_mfgraphs_plan_response.md` (Solver-boundary section), `issue_86_scenario_kernel.md`.

Acceptance gate:
```mathematica
s = makeScenario[GetExampleData[12] /. {I1 -> 100, U1 -> 0}];
r = SolveMFG[s];
IsFeasible[r]  (* must return True *)
```

---

### Phase 4 — Scenario catalog and representative migration

**Gate**: `scenarioCatalog[]` returns at least three named scenarios; one core, one inconsistent-switching.

Deliverables:
- `scenarioCatalog[]` registry concept: queryable by name, tier, tags, eligibility
- `scenarioByName[...]` resolution function
- `deriveScenario[parent, overrides]` with recursive merge semantics and lineage metadata (parent identity + scenario hash)
- Wave 1 migration: one core 4×4-type case and one inconsistent-switching case migrated to scenario objects
- Wave 2 migration: HRF Scenario 1 paper case migrated (only after Wave 1 assumptions are stable)
- MUnit tests for catalog lookup, lineage resolution, and merge semantics

Source inputs: `issue_87_scenario_catalog.md`, `revised_mfgraphs_plan_response.md` (Revised migration priority section).

Acceptance gate: `scenarioByName["HRF Scenario 1"]` returns a completed `scenario[...]` object that produces `Feasible` when passed to `SolveMFG`.

---

### Phase 5 — Extract scenario-oriented helper layers around DataToEquations

**Gate**: `DataToEquations` internal monolith is factored into named stages without changing its external contract.

Deliverables (provisional names — exact names decided during implementation):
- `NormalizeScenarioModel` — canonical model-plus-substitution bundle from scenario blocks
- `ValidateScenarioModel` — pre-compilation admissibility checks
- `BuildGraphData`, `BuildSwitchingCostData`, `BuildTransitionData`, `BuildKirchhoffData`, `BuildConstraintBlocks` — internal staged construction
- All stages covered by direct unit tests
- External `DataToEquations` contract unchanged

Source inputs: `public_api_full_surface_plan.md` (Phase 4 section).

Acceptance gate: fast suite passes; `DataToEquations` existing test coverage unchanged.

---

### Phase 6 — Freeze the post-transition public surface

**Gate**: revised symbol inventory with explicit keep/deprecate/retire decisions for all 104 original symbols.

Deliverables:
- Updated `public_api_symbol_inventory.md` with scenario-relevance column and post-transition decisions
- Compatibility aliases carry explicit deprecation language and migration notes
- Retired symbols removed from `BeginPackage` export list
- Example-parameter symbols (`I1`–`I3`, `U1`–`U3`, `S1`–`S16`) classified and either moved to examples-only namespace or explicitly exported with rationale

Source inputs: `public_api_full_surface_plan.md` (Phase 5 section), `public_api_gap_report.md`.

Acceptance gate: `Length[Names["MFGraphs`*"]]` is less than 80 (retired symbols removed); zero accidental exports in `Names["MFGraphs`*"]`.

---

### Phase 7 — Documentation and trust hardening

**Gate**: every surviving "core" and "advanced-retained" symbol has complete `::usage` and at least one MUnit test.

Deliverables:
- Contract-grade `::usage` for every surviving public symbol: input contract, output contract, defaults, failure behavior, stability note
- Generated `API_REFERENCE.md` reflects the stabilized surface
- `docs/planning/README.md` updated to reflect which planning docs are still active vs. superseded
- Quick-start README section demonstrates the scenario-first workflow end-to-end
- Maintenance gate: CI-equivalent check that catches accidental exports, missing `::usage`, and zero-test public symbols

Source inputs: `public_api_full_surface_plan.md` (Phase 6–7 sections), `PAPER_GENERATION_PLAN.md`, `PAPER_INTEGRATION_CHECKLIST.md`.

Acceptance gate: `wolframscript -file Scripts/RunTests.wls fast` reports zero failures; `GenerateDocs.wls` produces clean output with no missing-usage warnings.

---

## Acceptance Criteria and Release Gates

| Criterion | Phase gate |
|---|---|
| Scenario primacy: new user solves and inspects without ambient symbolic substitutions | Phase 3 |
| Legacy-symbol clarity: `alpha`, `j`, `u`, `z`, example parameters all have explicit decisions | Phase 6 |
| Stable public surface: inventory reflects post-transition API | Phase 6 |
| Documentation completeness: every surviving public function has meaningful usage text | Phase 7 |
| Trust completeness: every surviving public function has proportional test coverage | Phase 7 |
| Maintenance discipline: accidental exports and missing usage caught by checks | Phase 7 |
| Paper baseline maintained: HRF Scenario 1 still Feasible / <150s | Every phase |
| Fast suite zero-failure: 16-file fast suite passes throughout | Every phase |

---

## Risk Register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Scenario-first adaptation breaks existing raw-association workflows | Medium | High | Preserve backward-compatible overloads in `DataToEquations` and `SolveMFG`; add alias test in Phase 3 |
| Symbol retirement breaks downstream user code | Low–Medium | Medium | Compatibility symbols kept for at least one release cycle with deprecation warnings; migration notes in Phase 6 |
| Fictitious Play backend leaks into public surface prematurely | Low | Medium | FictitiousPlayBackend.wl has `(* INTERNAL ONLY *)` banner; no `::usage` for its symbols until issues #72–#75 resolved |
| NonLinear/Monotone re-exposure on master conflicts with critical-only policy | Medium | Medium | Re-expose only as `Method` options under `SolveMFG`; do not remove tests from the rollback branch before policy decision |
| Performance regression in DNFReduce from structural changes | Low | High | Run `CompareDNF.wls` before and after any change to `DataToEquations.wl` or `DNFReduce.wl`; record in `DNF_PERFORMANCE_HISTORY.md` |
| Scenario hash or lineage identity format becomes a breaking change | Low | Medium | Define hash format in Phase 4 with explicit versioning; treat format as `"Stability" -> "internal"` until Phase 6 freeze |

---

## Rollback and Compatibility Notes

- `master` at `5509995` is the last merged state. Never force-push master.
- The rollback branch `rollback/pre-critical-only-surface` at `c3711fa` is the working development base. All Phase 1–7 work should be done on feature branches off this branch or off a clean branch cut from it.
- The stash created before the rollback (labeled `pre-rewind WIP: MFGParallelMap + exit-pass-through test loosening + debug scripts`) is recoverable via `git stash pop`; it applies against the pre-refactor solver surface.
- Backward-compatible overloads for `DataToEquations` and `SolveMFG` must be retained through at least Phase 6 so existing scripts that pass raw associations continue to work.
