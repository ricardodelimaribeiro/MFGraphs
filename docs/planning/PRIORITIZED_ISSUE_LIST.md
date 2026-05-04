# MFGraphs Prioritized Issue List

**Author:** Manus AI  
**Date:** 2026-05-03  
**Repository phase:** Core scenario-kernel phase

## Executive summary

This issue list translates the recent repository review into a concrete implementation queue. The highest-value near-term work is **stabilization**, not feature expansion. The repository already has a coherent staged architecture centered on `makeScenario`, `makeSymbolicUnknowns`, `makeSystem`, and `solveScenario`, but the public boundary, test protection, and symbolic-performance discipline are not yet equally mature.[1] [2] [3]

The most important conclusion is that the next sprint should make the **current workflow smaller, clearer, and safer**. That means finishing the scenario-first public solve path, protecting it with targeted tests and benchmarks, and reducing the symbolic expression growth that is currently the main time-and-memory risk.[1] [2] [3]

| Priority band | Theme | Why it comes first |
|---|---|---|
| P0 | Public workflow stabilization | The repo should expose one canonical path before adding more surface area. |
| P1 | Performance and solver hardening | The symbolic solve path already shows timeout sensitivity and topology-dependent growth. |
| P2 | API and data-structure simplification | The architecture is good, but some compatibility and duplication costs remain. |
| P3 | Documentation and operational polish | Important, but lower leverage than stabilizing correctness and performance contracts first. |

## Recommended execution order

The recommended short-term sequence is shown below.

| Order | Issue | Outcome |
|---|---|---|
| 1 | Implement scenario-native public solve routing | Makes the scenario-first workflow fully real. |
| 2 | Add scenario-first regression tests | Protects the newly canonical path from drift. |
| 3 | Formalize the example-surface policy | Prevents overlapping example APIs from hardening accidentally. |
| 4 | Add benchmark regression gates for the active solver path | Converts performance observations into an enforceable contract. |
| 5 | Refactor solver preprocessing away from monolithic constraint rebuilding | Attacks the main time hotspot directly. |
| 6 | Delay flattening of complementarity and inequality blocks | Reduces symbolic growth and memory pressure before solving. |
| 7 | Trim duplicated payload in `mfgSystem` | Improves memory discipline and data ownership clarity. |
| 8 | Make exactification configurable or delayed | Opens a lighter-weight path for easy numeric cases. |
| 9 | Harden the repo guard scripts and operational checks | Improves maintainability and policy enforcement. |
| 10 | Consolidate architecture/API docs for the active phase | Keeps the docs aligned after the stabilization sprint. |

## P0 issues: canonical public workflow stabilization

### Issue 1 — Implement `SolveMFG[s_scenario, opts___]` as a first-class public path

| Field | Content |
|---|---|
| Priority | **P0** |
| Title | Implement scenario-native `SolveMFG` dispatch |
| Problem | The repository has a strong scenario kernel and a scenario-aware middle layer, but the canonical public solve path is still incomplete because `SolveMFG` does not accept a typed scenario object directly.[2] |
| Why it matters | Without this dispatch, the package remains architecturally scenario-first but operationally mixed, forcing users to know where the typed abstraction stops and where the older association-oriented boundary begins.[2] |
| Dependencies | None; this should lead the stabilization sprint. |
| Suggested scope | Add `SolveMFG[s_scenario, opts___]` and delegate through the existing system-building and solver-routing workflow rather than creating a parallel implementation path. |
| Acceptance criteria | `SolveMFG` accepts a `scenario[...]` object directly; the path is documented as supported; it delegates through the same core solver stack as `solveScenario`; existing association-based compatibility remains intact. |

This is the single highest-priority issue because it closes the most visible gap between the repository’s current architecture and its intended user experience.[1] [2]

### Issue 2 — Add regression tests for scenario-first compile and solve routing

| Field | Content |
|---|---|
| Priority | **P0** |
| Title | Add tests for `DataToEquations[s_scenario]` and `SolveMFG[s_scenario]` |
| Problem | The intended scenario-first path is implementation-leading but test-lagging. Current tests still protect mostly association-based and compatibility-oriented flows.[2] |
| Why it matters | A partially migrated public path without tests is likely to regress silently during later API cleanup or solver refactoring. |
| Dependencies | Issue 1. |
| Suggested scope | Extend active tests to cover scenario construction, scenario compilation, direct scenario solve routing, and package-loading expectations for the chosen public symbols. |
| Acceptance criteria | New tests fail before the implementation and pass after it; the active `fast` suite includes direct scenario solve coverage; package-loading expectations reflect the chosen active symbol surface. |

### Issue 3 — Freeze the example-surface policy

| Field | Content |
|---|---|
| Priority | **P0** |
| Title | Decide the canonical example API and deprecation policy |
| Problem | The repository currently has overlapping example pathways, including legacy raw example data, transitional scenario construction, and scenario-first example factories.[1] [2] |
| Why it matters | Overlapping example APIs create conceptual drag, documentation drift, and test ambiguity. |
| Dependencies | Issue 1 is independent; Issue 2 should reflect the decision made here. |
| Suggested scope | Decide which helper is canonical for quick starts, which helpers remain compatibility-only, and which ones are archived or explicitly deprecated. |
| Acceptance criteria | One short policy table exists in the docs; package tests reflect the policy; compatibility helpers are labeled consistently; inactive paths are either removed from the active surface or clearly marked as compatibility-only. |

## P1 issues: performance and solver hardening

### Issue 4 — Add benchmark regression gates for the active solver workflow

| Field | Content |
|---|---|
| Priority | **P1** |
| Title | Turn benchmark observations into regression gates |
| Problem | The repo already documents active solver benchmarks and records important topology-sensitive performance behavior, but this evidence is not yet enforced as a regression contract.[3] |
| Why it matters | The active DNF-first path already shows severe case sensitivity: `grid-2x3` is about 4.4× slower than `grid-3x2`, and raw `reduceSystem` times out on `example-12` even at 160 seconds.[3] |
| Dependencies | None. |
| Suggested scope | Formalize a small benchmark set and record wall time, status, result kind, and validation outcome for each case; use the benchmark scripts already present as the operational backbone. |
| Acceptance criteria | A lightweight benchmark policy exists; the representative case set is versioned; changes can be compared against a baseline; regressions are visible in review before merge. |

### Issue 5 — Refactor solver preprocessing away from repeated monolithic constraint rebuilding

| Field | Content |
|---|---|
| Priority | **P1** |
| Title | Replace repeated whole-system rescans in `reduceSystem` preprocessing |
| Problem | The current symbolic reduction path repeatedly rebuilds and rescans a large global Boolean constraint, which is the main time hotspot identified in the review.[1] |
| Why it matters | This design is simple to reason about, but it scales poorly as complementarity structure grows and is the most likely cause of avoidable solver time amplification.[1] [3] |
| Dependencies | None, though it should be benchmarked under Issue 4. |
| Suggested scope | Introduce an internal structured-block preprocessing path that separates deterministic linear elimination from branch-heavy residual solving. |
| Acceptance criteria | The new path preserves solver correctness on the active test suite; representative benchmark cases do not regress; internal preprocessing logic no longer relies on repeated substitution into one monolithic base constraint on every elimination round. |

### Issue 6 — Delay flattening of complementarity and inequality blocks until the solver boundary

| Field | Content |
|---|---|
| Priority | **P1** |
| Title | Keep system blocks structured longer to reduce symbolic growth |
| Problem | `makeSystem` and its builders currently materialize global conjunctions too early, especially in complementarity-heavy areas.[1] |
| Why it matters | Early `And@@` construction inflates expressions before the solver has a chance to exploit structure, increasing both time and memory pressure.[1] |
| Dependencies | Can proceed independently, but should be measured under Issue 4 and coordinated with Issue 5. |
| Suggested scope | Store raw lists or grouped records for flow inequalities, alternative flows, transition alternatives, and switching constraints, with flattened views created only when required. |
| Acceptance criteria | System objects retain readable structured blocks; solver entrypoints can consume the structured representation; byte-count growth on representative cases is reduced or at least does not worsen. |

### Issue 7 — Add staged profiling output to normal review practice for solver changes

| Field | Content |
|---|---|
| Priority | **P1** |
| Title | Make staged profiling part of solver-change review |
| Problem | The repository already has a staged profiler, but the typical review workflow still relies more on end-to-end pass/fail behavior than on structured phase-by-phase cost analysis.[3] |
| Why it matters | When symbolic behavior regresses, it is important to know whether the growth came from scenario building, unknown generation, system construction, preprocessing, or solving. |
| Dependencies | None. |
| Suggested scope | Add a documented “run this profiler before merging solver-sensitive changes” practice and a small summary template for recording results. |
| Acceptance criteria | Solver-sensitive pull requests or change notes can cite staged profiler output for at least one representative easy case and one representative hard case. |

## P2 issues: API and data-structure simplification

### Issue 8 — Reduce duplicated payload inside `mfgSystem`

| Field | Content |
|---|---|
| Priority | **P2** |
| Title | Trim duplicated topology and unknown metadata in `mfgSystem` |
| Problem | The system object currently carries substantial overlap with scenario topology and unknowns data, increasing memory use and weakening data ownership boundaries.[1] |
| Why it matters | The staged architecture is one of the repository’s main strengths; letting each stage own only the data it truly needs would preserve that strength and reduce payload size. |
| Dependencies | Best done after Issues 5 and 6, so solver-facing structure changes are already settled. |
| Suggested scope | Decide which fields belong canonically to `scenario`, which belong to unknowns, and which belong to the system record only; keep flattening helpers as compatibility bridges rather than default storage strategy. |
| Acceptance criteria | Top-level `mfgSystem` size is reduced on representative cases; field ownership is documented; downstream access still works through typed accessors. |

### Issue 9 — Make exactification configurable or delayed

| Field | Content |
|---|---|
| Priority | **P2** |
| Title | Avoid unconditional early exactification for all cases |
| Problem | Exact rationalization is currently applied early, which is symbolically safe but likely inflates expression complexity for easy numeric cases.[1] |
| Why it matters | A configurable or delayed exactification policy could preserve symbolic correctness where needed while allowing a lighter-weight path for simpler cases. |
| Dependencies | Best evaluated with Issues 4–6 in place. |
| Suggested scope | Introduce an option or internal mode that preserves the current exact path by default for symbolic review, while allowing delayed or selective exactification in safe cases. |
| Acceptance criteria | The behavior is explicitly documented; tests cover both the exact path and the configured alternative if one is exposed; no correctness regressions occur on the active suite. |

### Issue 10 — Clarify the active exported surface versus compatibility surface

| Field | Content |
|---|---|
| Priority | **P2** |
| Title | Audit exported names and classify them as active, compatibility, or archival |
| Problem | The loader and docs are cleaner than before, but some compatibility wrappers still blur the distinction between the durable public API and transitional helpers.[2] |
| Why it matters | A smaller and better-labeled public surface makes the repo easier to maintain and safer to document. |
| Dependencies | Issue 3. |
| Suggested scope | Review exported symbols, classify them, and update package-loading expectations and docs accordingly. |
| Acceptance criteria | A current public-API table exists; compatibility-only names are labeled as such; archived symbols stay outside the active loaded surface. |

## P3 issues: operational and documentation polish

### Issue 11 — Harden the active-surface guard scripts for portability and enforcement

| Field | Content |
|---|---|
| Priority | **P3** |
| Title | Make repo guard scripts path-portable and part of routine checks |
| Problem | The repository already includes policy-enforcement scripts for the active test surface, but the guard layer is still lightweight and path-local in places, which limits reliability as a repeatable repo check. |
| Why it matters | Policy scripts are most useful when they are portable, easy to run, and hard to bypass accidentally. |
| Dependencies | None. |
| Suggested scope | Remove machine-specific assumptions, document usage clearly, and consider including the guard in the standard pre-merge validation flow. |
| Acceptance criteria | The guard scripts run from repo root on a fresh checkout; they do not depend on author-specific absolute paths; their intended role is documented in the repository workflow. |

### Issue 12 — Consolidate active-phase architecture documentation

| Field | Content |
|---|---|
| Priority | **P3** |
| Title | Create a single current-architecture reference for the active phase |
| Problem | The repository now has improved README, development history, benchmark notes, and module docs, but the architecture is evolving quickly enough that drift remains a risk.[2] [3] |
| Why it matters | A concise current-architecture document would reduce onboarding cost and help contributors distinguish active pathways from compatibility and archive layers. |
| Dependencies | Best done after Issues 1–3 and 10, so the public surface is stable. |
| Suggested scope | Add one document explaining the active load graph, canonical workflow, public API, benchmark expectations, and intentionally inactive surfaces. |
| Acceptance criteria | The doc exists, is linked from the README, and remains consistent with `CLAUDE.md`, the package loader, and the active test suite. |

## Suggested milestone grouping

| Milestone | Included issues | Goal |
|---|---|---|
| Milestone A — Public workflow stabilization | 1, 2, 3 | Make the scenario-first path real, test-protected, and conceptually clear. |
| Milestone B — Solver-performance hardening | 4, 5, 6, 7 | Make symbolic performance measurable, diagnosable, and less fragile. |
| Milestone C — Data-structure simplification | 8, 9, 10 | Reduce payload duplication and clarify long-term public boundaries. |
| Milestone D — Documentation and ops polish | 11, 12 | Improve maintainability and contributor guidance after the core is stable. |

## Final recommendation

The repository should resist the temptation to expand breadth before stabilizing the current kernel. The right next move is a **stabilization sprint**: finish the scenario-first public solve workflow, protect it with tests, formalize the example policy, and then spend one focused pass reducing symbolic solver overhead. That sequence preserves the repo’s strongest current asset — a cleaner staged architecture — while addressing the two most important unresolved risks: **public-boundary incompleteness** and **solver cost growth**.[1] [2] [3]

## References

[1]: /home/ubuntu/mfgraphs_simplicity_performance_review.md "MFGraphs simplicity and performance review"
[2]: /home/ubuntu/mfgraphs_scenario_migration_code_review.md "MFGraphs scenario migration code review"
[3]: /mnt/desktop/MFGraphs/BENCHMARKS.md "MFGraphs benchmark policy and history"
