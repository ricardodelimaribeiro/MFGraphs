# Legacy Scripts Using `GetExampleData`

These scripts still depend on `GetExampleData` (directly or indirectly) and are **not part of the active scenario-first workflow**.
They remain available for manual runs during migration.

## Legacy Inventory

- `Scripts/BenchmarkSuite.wls`
- `Scripts/BottleneckReport.wls`
- `Scripts/CompareDNF.wls`
- `Scripts/GenerateReferenceHashes.wls`
- `Scripts/InstrumentJamarat.wls`
- `Scripts/JamaratFeasibilityPhaseDiagram.wls`
- `Scripts/Phase0_Extract.wls`
- `Scripts/Profile.wls`
- `Scripts/ProfileEquivCheck.wls`
- `Scripts/ProfileMFGSS.wls`
- `Scripts/ProfilePreprocessing.wls`
- `Scripts/ProfileReductionStrategies.wls`
- `Scripts/RunBenchmarkCI.wls` (indirect via `BenchmarkSuite.wls`)
- `Scripts/SmokeTestExactMode.wls`

## Migration Note

Target replacement path is scenario-based construction (e.g., `ScenarioByKey[...]`) with input-defined boundaries/switching costs.
