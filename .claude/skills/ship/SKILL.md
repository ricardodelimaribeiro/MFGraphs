---
name: ship
description: Ship the current work — auto-branch on master, run tests, commit, push, open PR
---

You are executing the `/ship` workflow for MFGraphs. Follow each step in order and halt on failures.

## Steps

### 1. Branch setup
Run `git branch --show-current`.
- If the result is `master`, create and switch to a new branch, then continue.
- Branch name format: `chore/ship-<YYYYMMDD-HHMMSS>`.
- Command: `git switch -c chore/ship-$(date +%Y%m%d-%H%M%S)`.

### 2. Working tree status
Run `git status --short` and list pending files.

### 3. Run fast tests
```bash
wolframscript -file Scripts/RunTests.wls fast
```
Parse pass/fail counts. If any fail, report details and stop.

### 4. Commit
Stage intended changes and commit with a Conventional Commit message:
- Format: `<type>(<scope>): <summary>` or `<type>: <summary>`
- Common types: `feat`, `fix`, `refactor`, `docs`, `test`, `chore`

### 5. Push
```bash
git push -u origin HEAD
```

### 6. Open PR
Open a PR to `master` using GitHub MCP (`create_pull_request`).
Fallback: `gh pr create --fill`.
Report the PR URL.
