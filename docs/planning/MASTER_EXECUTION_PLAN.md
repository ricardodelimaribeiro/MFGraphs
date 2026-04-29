# MFGraphs Master Execution Plan

**Status**: Active — authoritative over all other planning docs in this directory.
**Last updated**: 2026-04-28

---

## Program Goals and Success Metrics

The overriding goal is to move MFGraphs from a research-grade codebase with an accidental public surface toward a **typed, scenario-first package** with an intentional API, contract-grade documentation, and validation proportional to risk.

| Goal | Measurable success criterion |
|---|---|
| Typed scenario kernel in production | `makeScenario → validateScenario → completeScenario → scenario[...]` pipeline works end-to-end on the 34 built-in examples |
| Scenario-first entrypoints | High-level orchestration wrappers accept `scenario[...]` directly without requiring manual unknown or system construction |
| Intentional public surface | Exported symbol count is reduced from 104 to a deliberate keep/advanced/compatibility set with no accidental exports (`alpha$`, `j$` fixed) |
| Paper readiness maintained | HRF Scenario 1 continues to produce `Feasibility → Feasible` on the critical solver within 150s |
| Test coverage on surviving API | Every "core" and "advanced-retained" symbol has at least one direct MUnit test reference |
| No accidental regression | Fast test suite (11 files per RunTests.wls) passes at zero failures after each phase |

---

## Current-State Baseline (2026-04-28)

### Branch and solver surface

- Active development branch: `docs/cleanup` (this branch)
- `NonLinearSolver.wl`, `Monotone.wl`, and `TimeDependentSolver.wl` do **not** exist on this branch (archived or deleted).
- `rollback/pre-critical-only-surface` (at `c3711fa`) is a preserved local branch with the full solver suite — not the active branch.
- `master` is pinned at `5509995` (critical-only surface). **Do not force-push master.**

### Package structure (modular tools refactor, 2026-04-28)

| File | Role |
|---|---|
| `MFGraphs/MFGraphs.wl` | Bootstrap: ::usage declarations, single Get block, shared private helpers |
| `MFGraphs/scenarioTools.wl` | Typed scenario kernel (`makeScenario`, `validateScenario`, `scenarioQ`) |
| `MFGraphs/unknownsTools.wl` | Symbolic unknown bundle construction (`makeSymbolicUnknowns`) |
| `MFGraphs/systemTools.wl` | Structural equation system kernel (`makeSystem`) |
| `MFGraphs/solversTools.wl` | Core solvers (`reduceSystem`, `dnfReduceSystem`, `booleanReduceSystem`) |
| `MFGraphs/graphicsTools.wl` | Public visualization helpers |
| `MFGraphs/primitives.wl` | Geometric and graph primitives |
| `MFGraphs/examples.wl` | Scenario factories and example generators |
| `MFGraphs/archive/` | Legacy files: `DataToEquations.wl`, `SolveMFGDispatch.wl`, `Solvers.wl`, `DNFReduce.wl`, `FictitiousPlayBackend.wl` |

### Exported API surface

Current exported surface: **104 symbols** (see `public_api_symbol_inventory.md` for full table).

Known hygiene issues — status after Phase 2 (2026-04-22):
- `alpha$` and `j$` — neither exists in the current codebase (confirmed by grep; gate already satisfied)
- Duplicate `::usage` declarations — **resolved**: removed from sub-packages; canonical declarations now only in MFGraphs.wl BeginPackage
- Many "advanced" symbols exported without tests (60/104 have zero direct test references) — deferred to Phase 6/7

### Validation evidence

- Small tier (cases 1–6, 27): 21/21 passed across all three solvers (2026-04-14).
- Paper tier (HRF Scenario 1): CriticalCongestion = Feasible / ~122s (2026-04-14).
- Scenario-onward test subset: 74/74 passed (2026-04-22).

---

## Resolved Cross-Doc Conflicts

The following conflicts existed across the input planning documents. These decisions are now canonical.

| Conflict | Decision |
|---|---|
| **Solver scope for current release line** | Critical-congestion solver is the only public-API-grade solver for the current release. NonLinear, Monotone, and TimeDependentSolver remain in the codebase on `rollback/pre-critical-only-surface` and are tested but not part of the primary user-facing narrative. |
| **Scenario-first vs. document-current-symbols-first** | Scenario-first wins. Do not freeze the current 104-symbol surface as the final API. Classify symbols relative to their role in the scenario-first workflow. |
| **Paper claims** | Only validated artifacts may be cited as paper-ready. HRF Scenario 1 on CriticalCongestion is validated. |
| **Fictitious Play backend exposure** | Internal only until vulnerabilities (GitHub #72–#75) are resolved. |
| **MFGraphs.wl as bootstrap vs. as definition file** | Bootstrap only. Function definitions belong in submodules; MFGraphs.wl only holds package flags, usage, and `Get` calls. |

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
```

---

## Phase-by-Phase Implementation Sequence

### Phase 1 — Scenario workflow specification ✅ COMPLETED 2026-04-22

**Gate**: written artifact exists, no code changes.

Deliverables:
- Canonical end-to-end user workflow — **done**: `makeScenario → solveScenario → plot`
- Scenario transition matrix — **done**: all 104 symbols classified.
- Quick-start outline for future README rewrite — **done**.
- Validation boundary policy — **done**: structural errors at scenario boundary; infeasibility at solver result.

---

### Phase 2 — Fix export hygiene and duplicate usage declarations ✅ COMPLETED 2026-04-22

**Gate**: clean package load with zero shadow warnings on the core symbol set.

Deliverables:
- Resolve duplicate `::usage` declarations — **done**.
- Fix `MFGraphs.wl` structure: single Get block (correct files only).
- Segregate legacy monoliths into `archive/` and establish `*Tools.wl` modular architecture.

---

### Phase 3 — Orchestration layer for Scenario-to-Solution

**Gate**: A single entrypoint (e.g., `solveScenario[s_scenario]`) chains the modular pipeline.

Deliverables:
- `solveScenario` (or similar high-level wrapper) that automatically sequences `makeSymbolicUnknowns`, `makeSystem`, and `reduceSystem`.
- High-level orchestration for multi-population scenarios.
- Backward compatibility layer for legacy `SolveMFG` calls that routes to the new modular pipeline.
- MUnit tests verifying the full chain from `makeScenario` to a feasible result.
- Fast suite still passes at zero failures.

Acceptance gate:
```mathematica
s = makeScenario[GetExampleData[12] /. {I1 -> 100, U1 -> 0}];
r = solveScenario[s]; (* chained makeSymbolicUnknowns -> makeSystem -> reduceSystem *)
IsFeasible[r]  (* must return True *)
```

---

### Phase 4 — Scenario catalog and representative migration

**Gate**: `scenarioCatalog[]` returns at least three named scenarios; one core, one inconsistent-switching.

Deliverables:
- `scenarioCatalog[]` registry concept: queryable by name, tier, tags, eligibility.
- `scenarioByName[...]` resolution function.
- `deriveScenario[parent, overrides]` with recursive merge semantics and lineage metadata.
- Wave 1 migration: one core 4×4-type case and one inconsistent-switching case migrated to scenario objects.
- Wave 2 migration: HRF Scenario 1 paper case migrated.

---

### Phase 5 — Refine and Validate Modular Helpers

**Gate**: Modular tools (`unknownsTools.wl`, `systemTools.wl`, `solversTools.wl`) are fully validated against legacy outputs.

Deliverables:
- Comprehensive test suite for `makeSymbolicUnknowns` and `makeSystem` across all topology types.
- Performance parity check between archived `DataToEquations` and new modular path.
- All stages of compilation covered by direct unit tests.

---

### Phase 6 — Freeze the post-transition public surface

**Gate**: revised symbol inventory with explicit keep/deprecate/retire decisions for all 104 original symbols.

Deliverables:
- Updated `public_api_symbol_inventory.md` with scenario-relevance column.
- Compatibility aliases carry explicit deprecation language.
- Retired symbols removed from `BeginPackage` export list.

---

### Phase 7 — Documentation and trust hardening

**Gate**: every surviving "core" and "advanced-retained" symbol has complete `::usage` and at least one MUnit test.

Deliverables:
- Contract-grade `::usage` for every surviving public symbol.
- Generated `API_REFERENCE.md` reflects the stabilized surface.
- Quick-start README section demonstrates the scenario-first workflow end-to-end.
- Maintenance gate: CI-equivalent check that catches accidental exports, missing `::usage`, and zero-test public symbols.

---

## Acceptance Criteria and Release Gates

| Criterion | Phase gate |
|---|---|
| Scenario primacy: new user solves and inspects without ambient symbolic substitutions | Phase 3 |
| Legacy-symbol clarity: `alpha`, `j`, `u`, `z`, example parameters all have explicit decisions | Phase 6 |
| Stable public surface: inventory reflects post-transition API | Phase 6 |
| Documentation completeness: every surviving public function has meaningful usage text | Phase 7 |
| Trust completeness: every surviving public function has proportional test coverage | Phase 7 |
| Paper baseline maintained: HRF Scenario 1 still Feasible / <150s | Every phase |
| Fast suite zero-failure: All tests pass throughout | Every phase |

---

## Risk Register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Scenario-first adaptation breaks existing raw-association workflows | Medium | High | Preserve backward-compatible overloads; add alias test in Phase 3. |
| Performance regression in symbolic reduction from structural changes | Low | High | Run `BenchmarkReduceSystem.wls` before and after changes. |
| Modular tools introduce internal symbol leakage | Medium | Medium | Use strict `Internal`` context or package-private scoping. |
