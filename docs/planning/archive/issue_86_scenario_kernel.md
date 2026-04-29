## Summary
MFGraphs should adopt a **typed scenario kernel** modeled on the validate-then-complete object pattern used in `gridtools.m`. The target is a scenario lifecycle of:

`makeScenario -> validateScenario -> completeScenario -> scenario[...]`

This should become the canonical entry point for defining cases, instead of relying on ad hoc example records and script-local assumptions.

## Why this matters
- A typed `scenario[...]` wrapper gives downstream code a stable contract and avoids running solvers on raw, partially formed associations.
- A validate-then-complete workflow creates one canonical internal representation for benchmarks, solvers, graphics, and docs.
- The current refactor roadmap depends on this kernel: benchmark cleanup, solver segregation, and graphics cleanup all become easier once scenarios are first-class objects.

## Scope
- define the public planning contract for `makeScenario`, `validateScenario`, `completeScenario`, and `scenarioQ`
- define the canonical nested blocks for a scenario object, including at least `Model`, `Data`, `Validation`, `Benchmark`, `Visualization`, and `Inheritance`
- decide which scenario fields may be values versus callables, and how validation should check admissibility
- define how completed scenarios receive stable identity metadata such as scenario hashes or canonical-version markers
- keep this issue focused on the **scenario kernel**, not full migration of all existing examples

## Deliverables
- a design brief for the scenario kernel
- a canonical scenario schema proposal aligned with MFGraphs needs
- a clear typed-object contract for downstream consumers

## Out of scope
- full migration of all benchmark cases
- full benchmark runner rewrite
- moving graphics helpers into the public API

## Related issues
- This is the architectural prerequisite for a scenario catalog and scenario-based benchmark runner.
- Related to #82, #83, and #25, but not covered by them.
