# Troubleshooting

## Active Scenario-Kernel Issues

### Package does not load

**Symptom:** `Needs["MFGraphs`"]` fails.

Run from the repository root or add it to `$Path`:

```mathematica
PrependTo[$Path, "/path/to/MFGraphs"];
Needs["MFGraphs`"];
```

The active loader brings in `primitives.wl`, `scenarioTools.wl`,
`examples.wl`, `unknownsTools.wl`, `systemTools.wl`, `solversTools.wl`, and
`graphics.wl`.

### Scenario validation returns `Failure`

**Symptom:** `makeScenario[...]` or an example factory returns a `Failure[...]`.

Check the model shape:

```mathematica
model = <|
  "Vertices" -> {1, 2, 3, 4},
  "Adjacency" -> {{0,1,0,0}, {0,0,1,0}, {0,0,0,1}, {0,0,0,0}},
  "Entries" -> {{1, 200}},
  "Exits" -> {{4, 0}},
  "Switching" -> {{1, 2, 3, 5}, {3, 2, 1, 5}}
|>;

s = makeScenario[<|"Model" -> model|>];
```

Common causes:

- Missing one of `"Vertices"`, `"Adjacency"`, `"Entries"`, `"Exits"`, or `"Switching"`.
- Non-integer or non-positive vertex labels in helper factories.
- Adjacency dimensions that do not match the vertex list.
- Boundary values that are not numeric.
- Switching costs that are not numeric or `Infinity`.

### `reduceSystem` fails on a non-critical system

**Symptom:** `reduceSystem`, `dnfReduceSystem`, or `booleanReduceSystem` returns
a failure for non-critical congestion.

The current structural solvers support critical congestion only. Use
`Alpha -> 1` globally and avoid non-1 edge-specific `EdgeAlpha` values:

```mathematica
s = gridScenario[{3}, {{1, 10}}, {{3, 0}}, {}, 1];
sys = makeSystem[s];
sol = reduceSystem[sys];
```

### Hamiltonian `V` and `G` appear ignored

`V`, `G`, `EdgeV`, and `EdgeG` are validated and stored in the scenario schema,
but the active HJ-equation builder currently uses only `Alpha` / `EdgeAlpha`.
Treat those fields as schema-only until the Hamiltonian builder is extended.

## Testing Issues

### Test failures with "expected failure" messages

**Symptom:** Test output shows `FAILED` but indicates it's expected

**What this means:** The test validates that certain problems **correctly fail**. Examples:
- Switching costs that violate triangle inequality (problem has no valid equilibrium)
- Network configurations that produce infeasible systems

**Action:** No action needed — these are feature validations, not bugs

### Archive tests fail

**Symptom:** `RunTests.wls full` fails while `RunTests.wls fast` passes.

The `archive` suite contains legacy solver/API compatibility tests. It is
excluded from the active suite by design. Use it only when intentionally
working on archived solver surfaces.

### Active test command

```bash
wolframscript -file Scripts/RunTests.wls fast
```

`all` is currently an alias for `fast`.

## Benchmarking Issues

### Legacy benchmark command is missing

**Symptom:** `Scripts/BenchmarkSuite.wls` or `Scripts/BottleneckReport.wls`
cannot be found.

Those scripts now live in `Scripts/archive/` and depend on legacy solver-era
symbols. Use the active benchmark unless you are explicitly working on legacy
solver migration:

```bash
wolframscript -file Scripts/BenchmarkReduceSystem.wls --tag "my change"
```

## Getting More Help

### Enable comprehensive debugging

```mathematica
$MFGraphsVerbose = True;
s = getExampleScenario[7, {{1, 100}}, {{3, 0}, {4, 10}}];
sys = makeSystem[s];
sol = reduceSystem[sys];
isValidSystemSolution[sys, sol, "ReturnReport" -> True]
```

### Inspect typed records

```mathematica
scenarioData[s]
unknownsData[makeUnknowns[s]]
systemDataFlatten[sys]
```

### Review relevant documentation

- **CLAUDE.md** — Developer guide with architecture and debugging workflows
- **BENCHMARKS.md** — Performance profiling methodology and bottleneck details
- **API_REFERENCE.md** — Complete function signatures and options
- **README.md** — Quick start examples and configuration

### Report or discuss issues

Detailed issues can be reported on GitHub with:
1. Minimal reproducible example (data + solver call)
2. Output of `$MFGraphsVerbose = True` run
3. Mathematica version
4. OS and machine specs (optional but helpful)
