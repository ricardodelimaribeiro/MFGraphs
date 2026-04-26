---
name: ship
description: Delivery workflow for MFGraphs changes with automatic branch creation on master, mandatory fast-suite verification, and PR creation.
---

# /ship

Run this workflow in order:

1. Branch Setup
- Run `git branch --show-current`.
- If current branch is `master`, create and switch to a new feature branch before continuing.
- Branch name format: `chore/ship-<YYYYMMDD-HHMMSS>`.
- Use: `git switch -c chore/ship-$(date +%Y%m%d-%H%M%S)`.

2. Verification
- Run `wolframscript -file Scripts/RunTests.wls fast`.
- Abort immediately if any test fails.

3. Commit
- Stage all changes with `git add -A`.
- Commit using a Conventional Commits message such as `feat: ...` or `fix: ...`.

4. Delivery
- Push to origin: `git push -u origin HEAD`.
- Open a PR using the GitHub MCP `create_pull_request` tool.
- If MCP PR creation is unavailable, fallback to `gh pr create --fill`.
