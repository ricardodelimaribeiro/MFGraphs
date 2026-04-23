# Legacy Scripts Using `GetExampleData`

These scripts still depend on `GetExampleData` (directly or indirectly) and are **not part of the active scenario-first workflow**.
They are archived under `Scripts/archive/`.

## Legacy Inventory

- `Scripts/archive/BenchmarkSuite.wls`
- `Scripts/archive/BottleneckReport.wls`
- `Scripts/archive/CompareDNF.wls`
- `Scripts/archive/GenerateReferenceHashes.wls`
- `Scripts/archive/InstrumentJamarat.wls`
- `Scripts/archive/JamaratFeasibilityPhaseDiagram.wls`
- `Scripts/archive/Phase0_Extract.wls`
- `Scripts/archive/Profile.wls`
- `Scripts/archive/ProfileEquivCheck.wls`
- `Scripts/archive/ProfileMFGSS.wls`
- `Scripts/archive/ProfilePreprocessing.wls`
- `Scripts/archive/ProfileReductionStrategies.wls`
- `Scripts/archive/RunBenchmarkCI.wls` (indirect via `BenchmarkSuite.wls`)
- `Scripts/archive/SmokeTestExactMode.wls`

## Migration Note

Target replacement path is scenario-based construction (e.g., `ScenarioByKey[...]`) with input-defined boundaries/switching costs.
