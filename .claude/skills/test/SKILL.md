---
name: test
description: Run the MFGraphs fast test suite and report exact pass/fail counts
---

You are executing the `/test` workflow for the MFGraphs project.

## Steps

### 1. Run fast suite
```bash
wolframscript -file Scripts/RunTests.wls fast
```

### 2. Parse and report
Scan the output for lines matching patterns like:
- `X success` / `X failure` / `X skipped`
- `Tests passed: X` / `Tests failed: X`
- Any `FAIL` or `ERROR` lines with test names

Report a clean summary in this format:

```
Test suite: fast
  Passed:  <N>
  Failed:  <N>
  Skipped: <N>
```

If there are failures, print each failing test name and its error message.
Exit status: if any tests failed, say so clearly so the user knows the tree is red.
