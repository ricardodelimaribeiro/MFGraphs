# Revised MFGraphs Refactor Plan

This revision incorporates your feedback and stays strictly at the **planning** level. The main adjustment is that the earlier roadmap should lean more explicitly into the `gridtools.m` pattern by treating scenarios as **typed, validated objects** rather than only as canonical associations.

## Architectural revision

The most important refinement is that the future MFGraphs case layer should revolve around a **typed scenario object**, written conceptually as `scenario[<|...|>]`. The association remains the canonical payload, but the wrapper becomes the contract that downstream code consumes. This makes it possible for solver and benchmark entry points to accept only validated scenario objects, just as `gridtools.m` uses typed wrappers to separate raw input from completed domain objects.

| Revised principle | Planning implication |
|---|---|
| **Typed wrapper** | Downstream functions should target `scenario[...]`, not raw associations |
| **Validate then complete** | Scenario authoring remains flexible, but execution always uses canonical completed objects |
| **Lineage-aware inheritance** | Derived scenarios should preserve parent identity and override provenance |
| **Callable data support** | Data fields such as `RunningCost` should be allowed to be values or functions |
| **Solver as consumer** | `MFGSystemSolver` should consume scenarios rather than define them |

## Revised scenario-model direction

The schema proposal should be refined in three concrete ways.

First, the wrapper should be explicit. Instead of thinking only in terms of `makeScenario` returning an association, the planning target should be a constructor pipeline that produces `scenario[completedAssoc]`. This supports stricter interfaces such as `solveScenario[s_scenario]` and reduces accidental execution on partially formed case data.

Second, inheritance should carry **lineage metadata**. The earlier `Inheritance` block should be extended conceptually to preserve a parent identifier and a parent scenario hash. The purpose is not just provenance, but also reproducibility: benchmark tooling can detect whether a derived scenario is materially distinct from its parent or only cosmetically renamed.

Third, the `Data` block should explicitly allow **coordinate-aware or callable specifications**. In other words, fields like `RunningCost` should be permitted to be static values, vectors, or callables, provided the future validator checks that they are admissible for the intended discretization and evaluation pipeline.

| Revised schema area | Added planning requirement |
|---|---|
| `scenario[...]` wrapper | Typed public object, not only a completed association |
| `Inheritance` | Add lineage identity, including parent hash/provenance |
| `Data` | Permit callable or coordinate-aware inputs where appropriate |
| `Validation` | Include admissibility checks for callable fields |
| `Benchmark` | Use scenario identity and lineage in run metadata |

## Revised answer to the solver-boundary question

Your feedback sharpens the answer to the earlier architectural question about `MFGSystemSolver`.

> **`MFGSystemSolver` should be treated as a coordinator in the solver layer, not as part of scenario definition or equation-definition ownership.**

That means the planning target should place it under a solver-oriented module structure. Conceptually, it should receive a validated `scenario[...]` object, invoke `DataToEquations` to obtain the symbolic system, and then dispatch to specialized solver families according to scenario eligibility and congestion regime.

| Module area | Planned role |
|---|---|
| `Scenario` layer | Define, validate, complete, and catalog scenarios |
| `DataToEquations` | Transform a completed scenario into an equation system |
| `Solvers/MFGSystemSolver` | Coordinate stage selection and dispatch |
| `Solvers/CriticalCongestionSolvers` | Critical-congestion-specific methods |
| `Solvers/NonCriticalCongestionSolvers` | Non-critical-congestion-specific methods |

This clarifies the ownership boundary: **scenario construction defines what the problem is; solver modules decide how to solve it**.

## Revised roadmap

The work should be reordered slightly so that the typed scenario kernel and registry are established before broader solver or benchmark refactors. This makes the benchmark and solver cleanup a downstream consequence of the scenario architecture rather than a parallel redesign.

| Phase | Milestone | Revised planning focus |
|---|---|---|
| **1** | **Scenario Kernel** | Define the planning contract for `makeScenario`, `validateScenario`, `completeScenario`, and the `scenario[...]` typed wrapper |
| **2** | **Scenario Registry** | Introduce a scenario catalog concept and migrate a very small representative set first |
| **3** | **Solver Interface Cleanup** | Make `MFGSystemSolver` and related solver layers consume scenario objects |
| **4** | **Benchmark v2** | Rebuild benchmark selection around catalog queries such as tier, tags, and eligibility |
| **5** | **Graphics Reorganization** | Move visualization to a scenario-aware graphics layer |
| **6** | **Documentation Audit** | Rewrite docs after the scenario and solver boundaries are stable |

## Revised migration priority

The first migrated scenarios should be chosen to stress the architecture in different ways. I agree with the direction implied by your feedback: the initial set should include one core benchmark case and one inconsistent-switching case. I would add one paper-tier case only after the kernel and registry assumptions are stable enough to support a more complex scenario cleanly.

| Migration wave | Scenario type | Why it matters |
|---|---|---|
| Wave 1 | Core 4x4 grid | Small, representative, benchmark-relevant baseline |
| Wave 1 | Inconsistent-switching example | Exercises validation metadata and expected inconsistency behavior |
| Wave 2 | Paper-tier case | Exercises richer metadata, display labels, and benchmark defaults |

## Planning note on inheritance merge

Your note about deep merge is important and should be promoted from an implementation detail to a planning requirement. Derived scenarios should not replace whole nested blocks when only one parameter changes. The roadmap should therefore assume a **recursive association merge policy** for inheritance, especially in `Model`, `Data`, `Validation`, `Benchmark`, and `Visualization` blocks.

| Merge rule | Planning consequence |
|---|---|
| Shallow merge is insufficient | Nested blocks would be too easy to wipe accidentally |
| Deep merge required | Derived scenarios can override one field without destroying siblings |
| Revalidation after merge | The merged scenario must still pass full validation |
| Recomputed identity | Hashes and derived metadata must be refreshed after inheritance resolution |

## Revised immediate next step

The most appropriate next planning step is no longer a generic schema note, but a **small design brief for the Scenario Kernel**. That brief should define:

| Item | Planning target |
|---|---|
| Public typed object | `scenario[assoc]` |
| Constructor lifecycle | `makeScenario -> validateScenario -> completeScenario -> scenario[...]` |
| Canonical nested blocks | `Model`, `Data`, `Validation`, `Benchmark`, `Visualization`, `Inheritance` |
| Lineage contract | Parent reference plus scenario hash/provenance |
| Merge contract | Recursive inheritance merge with revalidation |
| Admissibility rules | Callable/value support for data fields |

## Final revised recommendation

The plan should now be framed around a stronger central claim:

> **MFGraphs should adopt a typed scenario kernel, with validated canonical scenario objects serving as the common input to equation generation, solver dispatch, benchmarking, visualization, and documentation.**

That refinement makes the later tasks more coherent. It also gives a cleaner answer to your original repository questions: benchmark cleanup follows from the catalog, solver segregation follows from the typed scenario consumer model, graphics cleanup follows from scenario metadata, and documentation cleanup follows once these boundaries are stable.

If you want, the next planning-only step can be a **revised Scenario Kernel brief** that narrows this further into public functions, nested blocks, and merge semantics, still without implementing any code.
