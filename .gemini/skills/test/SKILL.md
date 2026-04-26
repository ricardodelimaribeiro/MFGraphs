---
name: test
description: Run MFGraphs fast tests and report exact pass/fail totals and failing test files.
---

# /test

1. Execute:
- `wolframscript -file Scripts/RunTests.wls fast`

2. Parse and report:
- Exact totals from the runner output: `Total passed` and `Total failed`.
- Any failing test files listed under each `--- <file> ---` block where `Failed` is greater than 0.
- Any individual failed test IDs listed under `Failed tests:`.

3. Result format:
- `Passed: <n>`
- `Failed: <n>`
- `Failing files: <comma-separated list or none>`
- `Failing test IDs: <list or none>`
