---
name: ship
description: Ship the current branch — guard checks, tests, commit, push, open PR
---

You are executing the `/ship` workflow for the MFGraphs project. Follow every step in order; halt and report if any step fails.

## Steps

### 1. Branch guard
Run `git branch --show-current`. If the result is `master`, stop immediately and tell the user: "You are on master — create a feature branch first."

### 2. Working tree status
Run `git status --short`. If there are unstaged or untracked files the user may want to include, list them and ask whether to stage them or proceed with only the already-staged changes.

### 3. Run fast tests
```bash
wolframscript -file Scripts/RunTests.wls fast
```
Parse the output for pass/fail counts. If any test fails, print the failure details and **halt** — do not commit or push.

### 4. Commit
Stage any files the user confirmed in step 2, then commit with a [Conventional Commits](https://www.conventionalcommits.org/) message:
- Format: `<type>(<scope>): <short summary>`
- Common types: `feat`, `fix`, `refactor`, `docs`, `test`, `chore`
- Keep subject line ≤ 72 characters
- Add a body if the change warrants explanation

### 5. Push
```bash
git push -u origin HEAD
```

### 6. Open PR
Use the GitHub MCP tools (`mcp__github__create_pull_request`) to open a pull request against `master`. Title should match the commit subject. Body should include:
- A brief summary of what changed and why
- A test plan checklist
- The standard Claude Code attribution footer

Report the PR URL when done.
