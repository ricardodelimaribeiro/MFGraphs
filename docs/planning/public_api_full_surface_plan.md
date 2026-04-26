# MFGraphs plan: scenario-first API stabilization, legacy-symbol retirement, universal documentation, and trust hardening

This plan rewrites the earlier full-public-surface roadmap around a different architectural priority. The repository should **stabilize the scenario abstraction first**, then use that stabilized abstraction to decide which current names deserve long-term public status, which should remain advanced helpers, and which should be retired. That sequence is safer than immediately documenting and hardening every currently exported name, because the repository already contains a newer scenario kernel that can absorb responsibilities now spread across raw associations, example-parameter symbols, and symbolic helper heads.

The key strategic change is simple. **Do not freeze the current symbol list as the final API.** Instead, define the future workflow around `scenario`, `makeScenario`, `validateScenario`, `completeScenario`, `ScenarioData`, `SolveMFG`, `DataToEquations`, and plotting functions that can all accept the same scenario-first model. Once that boundary is stable, the repository can document the remaining surface with much greater confidence and much less churn.

The numerics package README reinforces this direction. Its strongest ideas are not language-specific; they are **API-shaping ideas**: lead with typed and validated abstractions, give users a short canonical workflow, separate “quick start” from deeper reference material, and present the package as a layered system rather than a pile of utilities. MFGraphs should borrow that discipline for its scenario-first redesign.

## Goal

The goal is to move MFGraphs toward a **scenario-first public API** in which user-facing workflows operate on a stable scenario object, legacy symbolic names are retired or demoted deliberately, and every surviving public function is documented and test-backed.

| Goal area | End-state requirement |
|---|---|
| Primary abstraction | Scenarios become the canonical user-facing model container |
| Legacy symbols | Names such as `alpha`, `j`, `u`, and example-parameter families are either retired, internalized, or clearly classified |
| Documentation | Every surviving public function has contract-grade usage documentation |
| Trust | Every surviving public function has verification proportional to its importance and risk |
| Simplification | `DataToEquations` and related helpers become easier to use through scenario-oriented helper layers |
| Reviewability | The transition lands in small, reversible changes rather than one disruptive rewrite |

## Why the plan changed

The current codebase already contains evidence that the long-term abstraction should not be “raw symbolic fragments plus ad hoc substitutions.” The scenario module defines a typed scenario kernel with clear lifecycle operations and a canonical block structure, including a dedicated `"Model"` block for raw topology data and a `"Data"` block for substitutions. At the same time, the example-data module still exports parameter symbols such as `I1`, `U1`, and `S1` specifically so user substitution rules continue to match built-in examples. That export pattern is useful today, but it also signals that the current surface mixes **workflow semantics** with **example convenience**.

> The repository should therefore stabilize the **container and workflow** before it stabilizes every currently visible symbol.

This is especially important for symbolic names such as `alpha`, `j`, `u`, and `z`. Some of them may remain valid internal mathematical heads or advanced symbolic utilities, but if the scenario layer becomes canonical then these names should not automatically inherit first-class public status merely because they are exported today.

## Revised planning principles

The rewritten plan is governed by five principles. First, **scenario-first beats symbol-first**. Second, **public must mean intentional** rather than merely exported. Third, **documentation should follow the stabilized workflow**, not transient implementation details. Fourth, **trust effort should concentrate on the future API**, not on names that may soon be removed. Fifth, **simplification should create usable layers**, especially around `DataToEquations`.

| Principle | Practical meaning |
|---|---|
| Scenario-first | Design the user workflow around scenario objects before finalizing the public symbol surface |
| Intentional public API | A name remains public only if it serves the long-term scenario-oriented workflow |
| Contract-first documentation | Usage text should describe inputs, outputs, defaults, and failure behavior for the stabilized API |
| Verification follows stability | Tests should first protect the APIs expected to survive the transition |
| Layered simplification | Break large functions into scenario-friendly helper stages instead of exposing raw internals indiscriminately |

## Non-goals for the next pass

The next pass should not attempt a full mathematical redesign of the solver, a broad benchmark campaign, or a promise that every currently exported helper is permanent. It should also avoid prematurely polishing legacy symbols that scenarios are likely to replace. The near-term objective is **API normalization**, not solver research and not repository-wide ornamentation.

## Proposed execution in seven phases

### Phase 1 — Define the scenario-first public workflow

The first phase should make the desired user journey explicit. The repository needs a normative answer to questions such as: What should a new user construct first? Should examples become scenarios immediately? Should `SolveMFG` and `DataToEquations` accept either raw model data or a `scenario[...]` object? How should plotting functions obtain metadata, substitutions, and visualization hints?

The output of this phase should be a workflow specification rather than code. That specification should define the canonical path from raw example or user input to completed scenario, compiled equations, solver result, and plots. It should also define the **one-page quick-start story** that the README will eventually tell, because the numerics README shows the value of presenting the package through a short, repeatable sequence of typed construction steps rather than through a flat list of symbols.

| Workflow decision | Target answer |
|---|---|
| Primary user input | `scenario[...]` or raw associations wrapped immediately into scenarios |
| Parameter handling | Substitutions live in the scenario `"Data"` block rather than in ambient symbol replacement habits |
| Compilation entrypoint | `DataToEquations` accepts scenarios directly and extracts `"Model"` plus `"Data"` consistently |
| Solve entrypoint | `SolveMFG` treats scenarios as the preferred high-level input form |
| Plotting entrypoint | Plotting helpers operate on solver results plus optional scenario metadata |
| Quick-start narrative | README demonstrates one canonical scenario creation, solve, inspect, and plot flow |
| Validation surface | Invalid scenarios fail early at the scenario boundary, not deep inside solver internals |

The concrete deliverables for this phase should therefore include both an internal workflow specification and a future-facing documentation skeleton.

| Phase 1 deliverable | Purpose |
|---|---|
| Scenario workflow specification | Defines the canonical end-to-end user flow |
| Scenario transition matrix | Maps current exported symbols to keep, advanced, compatibility, or retire outcomes |
| Quick-start outline | Defines the eventual README story for scenario-first use |
| Validation boundary notes | States which errors must be caught before compilation or solve |

### Phase 2 — Build a symbol transition policy

Once the workflow is defined, the second phase should classify the existing exported symbols according to how they relate to that workflow. This phase should answer which names remain core, which remain advanced, which survive only as compatibility aliases, and which are expected to disappear from the long-term public surface.

This is the phase where names such as `alpha`, `j`, `u`, `z`, `I1`–`I3`, `U1`–`U3`, and `S1`–`S16` should be judged carefully. Some may remain as advanced symbolic constructs. Others may move to examples-only documentation. Others may become compatibility names kept only long enough to support migration.

| Symbol class | Meaning in the rewritten plan |
|---|---|
| Core | Required for normal scenario construction, validation, solving, and plotting |
| Advanced | Supported for power users but not central to the primary scenario workflow |
| Compatibility | Kept temporarily for migration, with explicit deprecation language |
| Internal or retired | No longer part of the deliberate public contract once scenario-first routing is complete |

This phase should also establish the removal rule. A symbol should not survive merely because it is exported today; it should survive only if it remains useful and conceptually clear in a scenario-first package.

A practical test borrowed from the numerics README style is this: if a symbol cannot be explained naturally in the package quick start, the higher-level getting-started guide, or a clearly labeled advanced section, then it probably should not remain part of the main public contract.

### Phase 3 — Adapt the current entrypoints to scenarios before broad documentation

The third phase is the decisive implementation phase. Core functions should be adapted so they work naturally with scenarios before the repository invests in complete documentation coverage.

The most important adaptations are listed below.

| Function family | Required scenario-first behavior |
|---|---|
| `DataToEquations` | Accept a scenario directly, read `"Model"` and `"Data"`, and return a standardized compiled representation |
| `SolveMFG` | Prefer scenarios as the high-level entrypoint while preserving compatibility with raw compiled inputs |
| Plotting functions | Accept result objects derived from scenario-based solves without forcing users back into symbolic conventions |
| Example loading | Provide a path from `GetExampleData[...]` to `makeScenario[...]` without requiring manual symbol juggling |
| Validation helpers | Make scenario validation failures explicit and uniform before solver execution |

This phase is where the repository should reduce reliance on user-facing symbolic substitution patterns. The objective is not necessarily to remove all symbolic names immediately, but to make them **non-essential** to the primary workflow.

The implementation target should be a user experience where the scenario object plays the same organizing role that typed grids, unknown collections, layouts, and operators play in the numerics README: users advance through a small number of validated objects, each with a clear purpose, instead of managing loose symbolic conventions by hand.

### Phase 4 — Extract scenario-oriented helper layers around `DataToEquations`

After scenarios work at the entrypoint level, `DataToEquations` should be simplified into stages that align with the scenario-first model. The refactor should focus on creating helpers whose responsibilities match what a scenario-aware user or advanced contributor would naturally expect.

| Proposed helper family | Responsibility |
|---|---|
| `NormalizeScenarioModel` | Convert scenario blocks into a canonical model-plus-substitution bundle |
| `ValidateScenarioModel` | Check required keys, dimensions, and parameter completeness before compilation |
| `BuildGraphData` | Construct graph structures, partitions, and adjacency-derived state |
| `BuildSwitchingCostData` | Normalize switching-cost structures and consistency metadata |
| `BuildTransitionData` | Assemble transition relations, splitting/gathering structures, and switch alternatives |
| `BuildKirchhoffData` | Produce Kirchhoff matrices, signed-edge representations, and related linear data |
| `BuildConstraintBlocks` | Assemble equation and inequality blocks in a standardized form |
| `AssembleCompiledSystem` | Return the final `DataToEquations` association or an equivalent typed compiled object |

These names are intentionally provisional. The important point is that the decomposition should be driven by **scenario semantics** rather than by the accidental shape of the current monolith.

### Phase 5 — Freeze the post-transition public surface

Only after the scenario-first path is working should the repository freeze the final public API surface. At that point the inventory can be updated with much better confidence. Symbols that remain public will do so because they still matter after the transition, not because they happened to exist before it.

The output of this phase should be a revised public symbol inventory with explicit keep, deprecate, and retire decisions.

| Inventory field | Required meaning |
|---|---|
| Symbol | Exact public name |
| Role | Core, advanced, compatibility, or retired |
| Scenario relevance | Whether the symbol is central, peripheral, or replaced by scenario workflow |
| Documentation status | Complete, partial, or missing |
| Verification status | Covered, partially covered, or missing |
| Migration note | Whether users need a compatibility path |

### Phase 6 — Write universal documentation for the stabilized surface

Once the public surface is frozen, documentation can be written with much less risk of churn. Each public function should receive usage text that describes its contract in the context of the stabilized scenario-first package.

| Required documentation field | What it should state |
|---|---|
| Input contract | Accepted argument forms, especially scenario input support |
| Output contract | Return type, association keys, and expected invariants |
| Defaults and options | Automatic behavior and option interactions |
| Failure behavior | Expected failure envelopes, messages, or unresolved symbolic outcomes |
| Stability note | Core, advanced, or compatibility status |
| Minimal example | A scenario-oriented working example |

This phase should also reorganize generated docs so example-parameter symbols, if still present at all, do not dominate the main reference narrative.

The documentation set should eventually be split into the same kinds of layers that work well in the numerics README: a concise quick start, a project-structure overview, task-oriented guides for common workflows, and a lower-level API reference for advanced helpers.

### Phase 7 — Harden trust and establish maintenance gates

The final phase should harden correctness guarantees for the surviving public surface and prevent regression into accidental exports. Verification effort should be matched to the stabilized API, not to the discarded one.

| Function type | Minimum verification requirement |
|---|---|
| Core scenario workflow function | Happy-path regression, failure-path test, return-shape contract test, and one representative example-case test |
| Advanced helper | Direct unit test plus one edge-condition test |
| Compatibility symbol | Alias-equivalence test and explicit deprecation coverage |
| Retired path with migration shim | Migration test demonstrating replacement behavior or explicit failure |

The repository should then add lightweight gates for export hygiene, usage completeness, and symbol-to-test coverage so future changes do not silently reintroduce accidental names such as `alpha$` and `j$` or re-expand the public surface without policy review.

## Recommended near-term execution order

The revised plan changes the next action materially. The next branch should not start with blanket documentation work. It should start with scenario routing and symbol-policy decisions.

| Commit wave | Scope |
|---|---|
| Wave 1 | Add the revised plan and a symbol transition matrix based on the existing inventory |
| Wave 2 | Adapt `DataToEquations` and `SolveMFG` to accept scenarios cleanly |
| Wave 3 | Add migration shims from example-data workflows to scenario workflows |
| Wave 4 | Extract scenario-oriented `DataToEquations` helper layers without changing solver mathematics |
| Wave 5 | Freeze keep/deprecate/retire decisions for the public surface |
| Wave 6 | Write contract-grade docs for the stabilized surface |
| Wave 7 | Fill verification gaps and add maintenance gates |

## Acceptance criteria

The rewritten plan should be considered complete only when the repository can satisfy all of the following statements.

| Acceptance criterion | Description |
|---|---|
| Scenario primacy | A new user can solve and inspect a model through a scenario-first workflow without depending on ambient symbolic substitutions |
| Legacy-symbol clarity | Symbols such as `alpha`, `j`, `u`, `z`, and the example-parameter families have explicit keep, deprecate, or retire decisions |
| Stable public surface | The public inventory reflects the post-transition API rather than the historical export list |
| Documentation completeness | Every surviving public function has meaningful usage text and generated reference coverage |
| Trust completeness | Every surviving public function has tests proportional to its risk and importance |
| Maintenance discipline | New exports, missing usage text, and uncovered public functions are caught by repository checks |

## Recommended immediate next step

The immediate next step should be a **scenario transition matrix**, not a documentation sweep. That matrix should map each currently exported symbol to one of four outcomes: **core in the scenario workflow, advanced but retained, compatibility-only, or retired after migration**. Once that matrix exists, implementation can begin with much less ambiguity.

> In practical terms, the repository now needs a **scenario-first normalization phase**, not another round of general publicization.
