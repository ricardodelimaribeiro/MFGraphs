---
name: ship
description: Branch-guarded delivery workflow for MFGraphs changes with mandatory fast-suite verification and PR creation.
---

# /ship

Run this workflow in order:

1. Branch Guard
- Run `git branch --show-current`.
- Abort if the branch is `master`.

2. Verification
- Run `wolframscript -file Scripts/RunTests.wls fast`.
- Abort immediately if any test fails.

3. Commit
- Stage all changes with `git add -A`.
- Commit using a Conventional Commits message such as `feat: ...` or `fix: ...`.

4. Delivery
- Push to origin: `git push origin HEAD`.
- Open a PR using the GitHub MCP `create_pull_request` tool.
- If MCP PR creation is unavailable, fallback to `gh pr create --fill`.
