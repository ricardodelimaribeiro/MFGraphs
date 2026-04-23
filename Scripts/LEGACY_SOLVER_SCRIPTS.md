# Legacy solver-era scripts

The scripts listed here still depend on solver-era symbols such as
`DataToEquations`, `CriticalCongestionSolver`, and/or `GetExampleData`.

They are **not part of the active core-only workflow** while `MFGraphs.wl`
loads only Scenario/ExampleScenarios/Unknowns/System.

These scripts now live under `Scripts/archive/`.

## Typical affected scripts

- `archive/BenchmarkSuite.wls`
- `archive/BottleneckReport.wls`
- `archive/CompareDNF.wls`
- `archive/Profile*.wls`
- `archive/InstrumentJamarat.wls`
- `archive/GenerateReferenceHashes.wls`
- `archive/SmokeTestExactMode.wls`
- `archive/JamaratFeasibilityPhaseDiagram.wls`

## Active workflow

Use `Scripts/RunTests.wls fast` and the active test files in `MFGraphs/Tests/`.

## Re-enabling solver-era scripts later

When solver modules are restored to the default package load path, revisit and
reactivate these scripts by removing/adjusting this legacy designation.
