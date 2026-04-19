# Proposed MFGraphs Scenario Schema

This note proposes a **scenario-object architecture** for MFGraphs inspired by the design pattern used in `gridtools.m`: a raw specification is validated, completed into a canonical representation, and then consumed by downstream computational layers.

## Design objective

The central objective is to make every MFGraphs case a **first-class scenario object** rather than an ad hoc bundle of graph data, parameters, and benchmark-specific metadata. The scenario should be the single authoritative source for equation generation, solver dispatch, benchmarking, visualization, and documentation.

| Design principle | Meaning for MFGraphs |
|---|---|
| **Canonical schema** | Every scenario resolves to one stable internal representation |
| **Validate then complete** | Authors may write compact specs, but execution uses completed canonical objects |
| **Inheritance** | Variants can derive from parent scenarios and override only selected fields |
| **Derived artifacts are separate** | Equations, solver results, and plots should not be stored as raw case definitions |
| **Explicit metadata** | Tier, tags, assumptions, and expected properties should live with the scenario |

## Recommended object lifecycle

The most important architectural lesson from `gridtools.m` is the lifecycle:

> raw spec -> validator -> completer -> canonical object -> derived computational objects

For MFGraphs, I recommend the following analogous lifecycle.

| Stage | Suggested MFGraphs function | Output |
|---|---|---|
| Authoring | `makeScenario[assoc]` | User-facing scenario constructor |
| Validation | `validateScenario[assoc]` | Boolean plus diagnostics |
| Completion | `completeScenario[assoc]` | Canonical scenario association |
| Wrapping | `Scenario[assoc]` or equivalent typed wrapper | Stable typed object |
| Equation derivation | `makeEquationSystem[scenario]` or `DataToEquations[scenario]` | Equation-system object |
| Solver derivation | `solveScenario[scenario, stage]` | Solver-result object |
| Benchmark derivation | `benchmarkScenario[scenario, opts]` | Normalized benchmark result |
| Graphics derivation | `scenarioPlots[scenario, result]` | Graphics/visualization objects |

This lifecycle keeps case definition concerns separate from solver internals and benchmark orchestration.

## Proposed canonical scenario schema

I recommend dividing the schema into field groups. Some fields are **required authoring fields**, while others are **completed or derived fields** added automatically.

### 1. Identity and catalog metadata

These fields make the scenario addressable inside benchmarks, papers, tests, and documentation.

| Field | Required | Purpose |
|---|---|---|
| `"Name"` | Yes | Stable machine-readable identifier |
| `"DisplayName"` | No | Human-readable label for tables and plots |
| `"ScenarioGroup"` | No | Family label, such as `"Core"`, `"Paper"`, or `"Switching"` |
| `"Tier"` | No | Benchmark tier such as `"core"`, `"large"`, `"paper"` |
| `"Tags"` | No | Search/filter tags |
| `"Description"` | No | Short prose summary |
| `"Source"` | No | Provenance, e.g. paper section, example source, or benchmark family |
| `"Version"` | No | Optional semantic revision marker |

### 2. Structural model definition

These fields define the mathematical scenario itself.

| Field | Required | Purpose |
|---|---|---|
| `"Network"` | Yes | Canonical graph/network description |
| `"States"` | Optional | Explicit state list if not derivable from the network |
| `"Transitions"` | Optional | Explicit transition list if not derivable |
| `"Modes"` | Optional | Population classes, regimes, or agent types |
| `"TimeDomain"` | Optional | Horizon or temporal structure |
| `"Populations"` | Optional | Population specifications when multi-population is supported |
| `"BoundaryConditions"` | Optional | Terminal, initial, or boundary data |

The exact split between `"Network"`, `"States"`, and `"Transitions"` should depend on how much can be derived automatically. The better pattern is to store one authoritative structural representation and derive the rest during completion.

### 3. Model data and economics

These fields capture the costs and parameters that drive equation construction.

| Field | Required | Purpose |
|---|---|---|
| `"RunningCost"` | Usually | Cost structure used in HJB or flow equations |
| `"TerminalCost"` | Optional | Terminal objective data |
| `"SwitchingCosts"` | Optional | Switching or transition cost specification |
| `"CongestionData"` | Optional | Congestion coefficients, capacities, or related model parameters |
| `"Demand"` | Optional | Demand, mass, or forcing data |
| `"Parameters"` | No | General scalar/vector parameter association |
| `"Assumptions"` | No | Symbolic assumptions used by simplification or solver routines |

I strongly recommend collecting numerically meaningful constants under `"Parameters"`, while keeping large structural objects such as switching-cost maps and congestion models in dedicated fields.

### 4. Validation and expected properties

This group is where MFGraphs can encode known scenario-level invariants and expected behavior.

| Field | Required | Purpose |
|---|---|---|
| `"ExpectedConsistency"` | No | Expected consistency status for switching or structural checks |
| `"ExpectedProperties"` | No | Known properties such as symmetry, monotonicity, or inconsistency |
| `"Eligibility"` | No | Which solvers or benchmark stages apply |
| `"KnownIssues"` | No | Explicit caveats for documentation and benchmarking |
| `"ReferenceResults"` | No | Historical expected outputs, timings, or structural summaries |

This field group is especially useful for your recent inconsistent-switching examples, because the scenario can directly declare the expectation that consistency should fail and that extra alternative-flow structure should appear.

### 5. Benchmark and workflow metadata

These fields help centralize benchmark handling instead of scattering it across scripts.

| Field | Required | Purpose |
|---|---|---|
| `"Benchmark"` | No | Benchmark-specific metadata bundle |
| `"DefaultStages"` | No | Recommended stages such as `DataToEquations` or `CriticalCongestion` |
| `"Timeouts"` | No | Suggested per-stage timeouts |
| `"IsolationMode"` | No | Preferred execution mode for benchmarks |
| `"ResultLabels"` | No | Stable labels for CSV/JSON export |

A good design is for `"Benchmark"` itself to be a nested association, rather than flattening every operational field into the top level.

### 6. Visualization metadata

These fields allow plots to be driven by scenario metadata rather than by notebook-local conventions.

| Field | Required | Purpose |
|---|---|---|
| `"Visualization"` | No | Nested defaults for layout, colors, ordering, and labels |
| `"PlotLabels"` | No | Scenario-specific naming for states, edges, or controls |
| `"DisplayOrder"` | No | Stable ordering for tables and charts |
| `"PublicationStyle"` | No | Defaults for publication-ready rendering |

This helps move graphics logic out of `GettingStarted.wl` and into a reusable layer.

### 7. Inheritance and derivation metadata

This is the part most directly inspired by the parent-child pattern in `gridtools.m`.

| Field | Required | Purpose |
|---|---|---|
| `"ParentScenario"` | No | Reference to base scenario |
| `"Overrides"` | No | Explicit override record for changed fields |
| `"DerivedFrom"` | No | Provenance string or symbolic link to a scenario family |
| `"VariantKind"` | No | E.g. `"BenchmarkVariant"`, `"PaperVariant"`, `"InconsistentSwitchingVariant"` |

The key rule should be that inherited scenarios are still completed into the **same canonical schema** as standalone scenarios.

## Suggested canonical top-level shape

A practical top-level schema could look like this.

```wl
<|
  "Name" -> "Grid0404",
  "DisplayName" -> "4x4 Grid",
  "ScenarioGroup" -> "Large",
  "Tier" -> "large",
  "Tags" -> {"grid", "benchmark", "switching"},
  "Description" -> "Large benchmark grid scenario.",

  "Model" -> <|
    "Network" -> ..., 
    "Modes" -> ...,
    "TimeDomain" -> ...,
    "BoundaryConditions" -> ...
  |>,

  "Data" -> <|
    "RunningCost" -> ...,
    "TerminalCost" -> ...,
    "SwitchingCosts" -> ...,
    "CongestionData" -> ...,
    "Demand" -> ...,
    "Parameters" -> <||>,
    "Assumptions" -> True
  |>,

  "Validation" -> <|
    "ExpectedConsistency" -> True,
    "ExpectedProperties" -> {},
    "Eligibility" -> <|
      "DataToEquations" -> True,
      "CriticalCongestion" -> True,
      "NonLinearSolver" -> True,
      "Monotone" -> False
    |>,
    "KnownIssues" -> {}
  |>,

  "Benchmark" -> <|
    "DefaultStages" -> {"DataToEquations", "CriticalCongestion"},
    "Timeouts" -> <|"CriticalCongestion" -> 120|>,
    "IsolationMode" -> "Process"
  |>,

  "Visualization" -> <|
    "PlotLabels" -> <||>,
    "DisplayOrder" -> Automatic,
    "PublicationStyle" -> "Default"
  |>,

  "Inheritance" -> <|
    "ParentScenario" -> None,
    "Overrides" -> <||>,
    "VariantKind" -> None
  |>,

  "Derived" -> <|
    "ScenarioHash" -> ...,
    "CanonicalVersion" -> 1
  |>
|>
```

I recommend this **nested** structure over a flat one, because it separates semantic roles cleanly and makes documentation easier. The top level should remain compact, while each subsystem reads only the nested block it needs.

## Constructor and helper API

The schema becomes much more useful if it is paired with a small public API.

| Function | Role |
|---|---|
| `makeScenario[assoc_]` | Public constructor from raw authoring spec |
| `validateScenario[assoc_]` | Shape and semantic validation |
| `completeScenario[assoc_]` | Fill defaults, derive metadata, normalize fields |
| `scenarioQ[obj_]` | Predicate for typed scenario objects |
| `scenarioCatalog[]` | Return registry of named scenarios |
| `scenarioByName[name_]` | Resolve one scenario from registry |
| `deriveScenario[parent_, overrides_]` | Create a variant from a base scenario |

This is the direct analogue of the `makeGrid` / `validateGrid` / `completeGrid` pattern.

## Recommended inheritance model

I would keep inheritance intentionally simple at first. A derived scenario should reference a parent and provide only a shallow override association. Completion should merge parent and child definitions, then re-run validation on the completed result.

| Inheritance rule | Recommendation |
|---|---|
| Parent lookup | Resolve parent before validation of derived fields |
| Merge policy | Deep merge nested associations by key |
| Validation timing | Validate both override shape and final merged scenario |
| Provenance | Preserve `ParentScenario` and explicit `Overrides` |
| Derived metadata | Recompute hashes, derived defaults, and benchmark hints after merge |

This will let you define families of cases without duplicating full records.

## Scenario categories that would benefit immediately

Several existing MFGraphs workflows appear to map naturally onto this schema.

| Current MFGraphs need | Scenario-schema benefit |
|---|---|
| Core benchmark cases | Stable `Tier`, `Tags`, and `Benchmark` metadata |
| Paper scenarios | Better `DisplayName`, `Source`, and publication metadata |
| Inconsistent-switching examples | Direct `ExpectedConsistency -> False` and expected-structure annotations |
| Large process-isolated cases | Benchmark defaults live with the scenario instead of in runner branches |
| Getting-started examples | Clean visualization and description metadata |

## Separation from downstream artifacts

One of the most important structural rules should be that scenarios do **not** directly store full solver outputs or benchmark tables as part of the canonical definition. Those should be separate derived objects.

| Object | What it should contain |
|---|---|
| `Scenario[...]` | Canonical model definition and metadata |
| `EquationSystem[...]` | Variables, equations, and structural summaries |
| `SolverResult[...]` | Stage-specific output and diagnostics |
| `BenchmarkResult[...]` | Timings, statuses, and export-ready metadata |
| `VisualizationResult[...]` | Generated graphics or plot specifications |

That separation will help answer your later question about whether `MFGSystemSolver` belongs in equation generation or the solver layer.

## Minimal migration strategy

I recommend a staged migration rather than a one-shot rewrite.

| Step | Action |
|---|---|
| 1 | Define `makeScenario`, `validateScenario`, and `completeScenario` |
| 2 | Create the canonical schema and nested field groups |
| 3 | Migrate three representative cases: one core case, one inconsistent-switching case, and one paper case |
| 4 | Add compatibility wrappers so current benchmark and solver code can still consume old case definitions temporarily |
| 5 | Move benchmark selection to use scenario catalog queries by tier/tag/name |
| 6 | Update documentation and graphics utilities to consume scenario metadata |

This minimizes disruption and keeps the refactor testable.

## Recommended first implementation scope

To keep the first step small, I would initially implement only the following required authoring fields:

| Initial required fields | Reason |
|---|---|
| `"Name"` | Stable identifier |
| `"Model" -> "Network"` | Core structural input |
| `"Data" -> "RunningCost"` | Essential model content |
| `"Data" -> "SwitchingCosts"` | Essential for your current scenario families |
| `"Benchmark" -> "DefaultStages"` | Immediate benchmark usefulness |

Everything else can be completed by defaults or added gradually.

## Recommendation

The best MFGraphs adaptation of `gridtools.m` is **not** to copy its exact code style, but to adopt its object discipline: canonical validated records, explicit completion, optional inheritance, and strict separation between source objects and derived computational objects. If you do that, the later cleanup of benchmarks, solver boundaries, graphics, and documentation will become much easier.

## Suggested immediate next step

The next practical step should be to draft a first `Scenario` API file with:

| Item | Purpose |
|---|---|
| `makeScenario` | Public constructor |
| `validateScenario` | Schema and semantic checks |
| `completeScenario` | Defaulting and canonicalization |
| `scenarioCatalog` | Central registry of named scenarios |
| `deriveScenario` | Parent-child variant creation |

Once that exists, you can migrate a very small subset of cases and let the benchmark reorganization follow from the new schema.
