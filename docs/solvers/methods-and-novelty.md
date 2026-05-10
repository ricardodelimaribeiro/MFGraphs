# Solver Methods and Package Novelty

This note separates the algorithms provided by Wolfram Language from the
MFGraphs-specific modeling and solver pipeline built around them.

## Solver Method Map

| MFGraphs function | Method used |
|---|---|
| `reduceSystem` | Wolfram `Reduce[constraints, allVars, Reals]`. For polynomial real systems, Wolfram documents this as using cylindrical algebraic decomposition (CAD) / real quantifier elimination machinery. |
| `dnfReduceSystem` | Package-defined recursive DNF reducer: solve equalities, substitute the resulting rules, and distribute over disjunctions. Equality solving uses Wolfram `Solve`. |
| `optimizedDNFReduceSystem` | Package hybrid: branch-state DNF reduction for small residual variable sets, otherwise the package `dnfReduce` reducer. |
| `activeSetReduceSystem` | Package active-set enumeration for complementarity alternatives, with exact linear substitution. Larger residual systems fall back to `dnfReduce`. |
| `booleanReduceSystem` | Wolfram `BooleanConvert[..., "DNF"]`, followed by Wolfram `Reduce` on each disjunct. |
| `booleanMinimizeSystem` | Wolfram `BooleanMinimize[..., "DNF"]`, followed by Wolfram `Reduce` on each disjunct. This is exact Boolean minimal-DNF/SOP minimization. It is analogous to classical Quine-McCluskey/Petrick-style two-level minimization, but Wolfram does not document `BooleanMinimize` as specifically using QMC or Petrick internally. |
| `booleanMinimizeReduceSystem` | Package hybrid: prune infeasible Boolean arms with `FindInstance`, decompose the remaining formula by variable sharing, call `BooleanMinimize[..., "DNF"]` per component, then call `Reduce` per disjunct. |
| `findInstanceSystem` | Wolfram `FindInstance[constraints, allVars, Reals]` after package preprocessing. |
| `directCriticalSystem` | Package graph-distance heuristic that forbids flow directions moving farther from every exit, then solves feasibility with Wolfram `FindInstance`. |
| `flowFirstCriticalSystem` | Package flow-first path: build a linear equality system, compute a `LeastSquares` candidate, validate it, and fall back to Wolfram `FindInstance`. |

## Important Private Algorithms

| Helper | Method used |
|---|---|
| `buildSolverInputs` | MFGraphs constraint assembly plus equality preprocessing. |
| `rulesFromEqualities` | Linear elimination through `CoefficientArrays` and `RowReduce` on the augmented matrix. |
| `accumulateEqualityRules` | Fixed-point equality elimination using repeated `rulesFromEqualities`. |
| `branchStateReduce` | Package branch-state enumeration: carry rules and residual constraints directly instead of expanding the entire DNF expression first. |
| `pruneDisjunctiveArms` | Feasibility pruning with Wolfram `FindInstance`. |
| `decomposeByVariableSharing` | Connected-component decomposition on a variable co-occurrence graph using `Graph` and `ConnectedComponents`. |
| `BuildLinearSystemFromEqualities` | Linear matrix extraction via `CoefficientArrays`. |
| `LinearSolveCandidate` | Numeric least-squares solve via Wolfram `LeastSquares`. |

## What Is Novel in MFGraphs

The package should not claim novelty for Wolfram built-ins such as `Reduce`,
`FindInstance`, `BooleanMinimize`, `RowReduce`, or `LeastSquares`. Those are
backend algorithms. The package contribution is the MFG-specific construction
and orchestration around those tools.

Defensible novelty claims are:

- A typed scenario-to-system workflow for stationary mean-field games on
  networks: `makeScenario -> makeSymbolicUnknowns -> makeSystem -> solveScenario`.
- An auxiliary state-space graph that turns physical network edges and turning
  choices into explicit state pairs and transition arcs.
- A critical-congestion symbolic encoding where Kirchhoff balance, boundary
  conditions, edge-flow relations, switching optimality, and exit
  complementarity are represented as an exact logical-algebraic system.
- Custom complementarity-aware preprocessing that extracts linear equality
  rules before the expensive Boolean or real-algebraic solve stage.
- Package-defined DNF and active-set reducers that exploit the structure of the
  MFGraphs system rather than sending the full formula directly to `Reduce`.
- Hybrid Boolean/real solver variants that treat complementarity as Boolean
  structure, then solve each surviving real-algebraic branch.
- Solver diagnostics that classify residual branches, transition-flow
  determinacy, Kirchhoff residuals, and branch objective costs.
- Paired raw and rich network visualizations that expose both the physical
  infrastructure and the augmented state-space graph.

Concise statement:

> MFGraphs provides an exact symbolic workflow for critical-congestion
> stationary mean-field games on networks by transforming the MFG optimality
> and conservation laws into a structured logical-algebraic system and solving
> it through MFG-specific preprocessing, DNF/active-set reductions, Boolean
> decomposition, and Wolfram real-algebraic backends.

