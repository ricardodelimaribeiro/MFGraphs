# MFGraphs Development History

*A chronicle of **864 commits** across almost six years — from a research prototype to a scenario-first, benchmarked Mean Field Games toolkit.* [1]

---

## 1. What this repository became

MFGraphs began on **2020-05-06** as a small Wolfram Language research repository and, by **2026-04-30**, had accumulated **864 commits**. The most striking fact in the current history is not only the project’s age, but its acceleration: **409 commits landed in 2026 alone, including 105 in March and 304 in April**.[1] The repository’s development history is therefore not a smooth six-year line. It is a long research incubation followed by a very sharp engineering transition.

At a high level, the repository moved through four recognizable phases. The first phase built the mathematical model and the earliest solver ideas. The second expanded the example set and solver experiments but left much of the repository in research-code form. The third, in March 2026, professionalized the project through pull requests, benchmarking, test infrastructure, and aggressive cleanup. The fourth, in April 2026, reorganized the active package surface around typed scenarios, modular system builders, orchestration helpers, and a narrower, more explicit public workflow.[1] [2]

| Phase | Time window | Character |
|---|---|---|
| Research prototype | 2020 | Rapid mathematical experimentation, minimal process |
| Expansion and drift | 2021–2024 | New applications and solvers, but uneven software discipline |
| Professionalization sprint | March 2026 | PR workflow, benchmarks, CI, cleanup, context hygiene |
| Architecture reset | April 2026 | Scenario-first core, modular typed kernels, orchestration, diagnostics |

---

## 2. The research phase (2020)

The repository was born with the aptly literal commit message **`First commit`** on 2020-05-06.[1] Two days later came a burst of commits all titled **`Association`**, which in retrospect captured one of the project’s best early decisions: representing network data as a Wolfram `Association`. That choice gave the package a flexible, self-describing modeling format early enough that many later refactors could preserve the underlying conceptual model even while the surrounding solver architecture changed.

The 2020 history reads like a direct research notebook translated into git. Commit messages from the period included remarks such as *“Having heart attacks and recovering quickly!”*, *“problem with jays”*, and *“getting to the final idea. found a way to reduce the equations more efficiently.”* The value of those commits is historical rather than procedural: they show the repository learning what its hard problems actually were. Those hard problems were not generic Wolfram coding issues. They were the algebraic shape of equilibrium constraints, the behavior of flow variables `j[...]`, and the combinatorial impact of switching-cost logic.[1]

Two milestones from late 2020 remained important for years. First, the non-linear solver became viable in August 2020, establishing that the repository could move from model construction into genuine equilibrium computation. Second, November 2020 introduced the `ZAnd` family of Boolean-reduction ideas that later evolved into `DNFReduce`. That moment matters because it identified the central computational bottleneck correctly: once switching costs generate nested `Or` structure, symbolic simplification and branch control dominate performance.

| 2020 milestone | Why it mattered later |
|---|---|
| Association-based network model | Became the durable conceptual representation for scenarios and solver inputs |
| First workable non-linear solver path | Proved the package could do more than symbolic derivation |
| Early `ZAnd` / logical reduction experiments | Foreshadowed the later `DNFReduce` optimization story |
| Switching-cost support | Made the package mathematically richer and computationally harder at the same time |

The early phase also established several liabilities that would remain with the repository for years. Development was notebook-heavy, tests were absent, and every change landed directly on `master`. In other words, the mathematics was advancing faster than the software process. That is normal for research code, but it set up the large cleanup bill that would arrive later.

---

## 3. Expansion, applications, and drift (2021–2024)

The years 2021 through 2024 broadened the package’s ambitions but did not yet fully professionalize it. The git counts show that activity remained substantial in 2021 and 2022, then dropped sharply in 2023 before reviving in 2024.[1] This was the repository’s long middle: there was real intellectual progress, but process maturity lagged behind.

In 2021 the repository grew through concrete applications, especially **Braess paradox** examples and **Jamarat** crowd-routing scenarios. Those examples mattered because they forced the solver to confront network structures that were more realistic than minimal chains and Y-graphs. The first `.mt` tests also appeared in September 2021. Even though the test surface remained small, that shift from visual notebook checking to explicit assertions was an important change in what counted as “done.”

The quieter 2022–2024 period was less dramatic in commit volume, but it still produced substantive solver work. `D2E2` continued to evolve, `ZAnd` remained an object of performance concern, and June 2024 introduced the **Monotone** solver path, a mathematically distinct alternative based on gradient-flow ideas. That move is historically important because it showed the repository was no longer trying to express the entire project through one solver philosophy alone. It was becoming a family of approaches around the same modeling domain.

| Year | Commits | Historical reading |
|---|---:|---|
| 2021 | 187 | Broadest early expansion: examples, solver experiments, and first tests |
| 2022 | 89 | Continued refinement, especially around `D2E2` and reduction concerns |
| 2023 | 1 | Near-hibernation of active development |
| 2024 | 18 | Selective revival, including Monotone-solver work |

What this middle phase lacked was not intelligence or technical depth. It lacked a disciplined repository contract. There was still too much duplication, too little benchmark evidence, and too little separation between exploratory code and maintained public surface. Those omissions are exactly what made March 2026 so consequential.

---

## 4. The professionalization sprint (March 2026)

March 2026 was the repository’s turning point. The history records **105 commits** and **28 merged pull requests** in that month alone.[1] The shift was not cosmetic. It changed how the project was developed, measured, and explained.

This was the month in which MFGraphs became a repository with a repeatable engineering process. Benchmark suites were introduced and performance history began to be written down systematically rather than inferred from memory.[3] [4] The old logical-reduction machinery was renamed and reframed around `DNFReduce`, cleanup removed stale notebooks and renamed modules to modern `.wl` conventions, and a wave of fixes addressed Wolfram-specific context hygiene bugs that had likely been latent for years.

One of the most important historical lessons from this month is that **measurement came before most of the high-confidence optimization claims**. That sequencing mattered. Once benchmarking existed, speedups and regressions could be discussed concretely rather than rhetorically. The repository started preserving not just code, but evidence.

| March 2026 theme | Historical effect |
|---|---|
| Pull-request workflow becomes normal | Replaced direct-to-`master` habits with reviewable increments |
| Benchmark and profiling scripts land | Turned solver work into a measurable engineering activity |
| `DNFReduce` performance work | Made symbolic-branch handling a first-class optimization topic |
| Cleanup and file modernization | Reduced notebook debt and clarified module boundaries |
| Context-hygiene fixes | Exposed how much Wolfram package correctness depends on symbol discipline |
| Test runner and CI hardening | Made regression detection sustainable rather than ad hoc |

By the end of March 2026, the repository had stopped being merely “a codebase that can solve some cases.” It had become a codebase that could explain itself, test itself, and benchmark itself.

---

## 5. The architecture reset (April 2026)

If March 2026 professionalized the repository, **April 2026 redefined its active center of gravity**. The git history shows **304 commits** and **93 merged pull requests** in April alone, making it the single most intense month in the repository’s life.[1]

The dominant pattern in April was not a single optimization or one solver breakthrough. It was **architectural narrowing and modularization**. The history records typed scenario-kernel work, scenario-topology refactors, typed unknown bundles, extraction of the system kernel, public constructor additions such as `GridScenario`, `CycleScenario`, `GraphScenario`, and `AMScenario`, archiving of legacy components, and later the addition of orchestration helpers, graphics improvements, diagnostics, and benchmark-history documentation.[1]

In practical terms, April 2026 changed the answer to the question “What is MFGraphs right now?” Before April, that answer was still heavily framed by the older `DataToEquations` / critical-solver / nonlinear-solver lineage. After April, the active answer became much more explicit: MFGraphs is a **scenario-first toolkit with typed kernels and a narrower active public workflow**, while older components remain part of the project’s historical record and archived capability surface.[1] [2]

| April 2026 cluster | Representative history entries |
|---|---|
| Typed scenario core | `typed scenario kernel and modular load order`; `extend Hamiltonian defaults and edge params` [1] |
| Topology and unknown refactor | `Refactor scenario topology lifecycle`; `Add typed unknowns head`; `Rename symbolic unknown bundle API` [1] |
| System modularization | `Extract system kernel`; `modularize mfgSystem builders` [1] |
| Public constructor expansion | `Add GridScenario`; `Add CycleScenario, GraphScenario, AMScenario` [1] |
| Active-surface narrowing | `Archive DataToEquations`; `Reduce MFGraphs loader to scenario core`; `archive legacy components` [1] |
| Workflow layer growth | `implement orchestration layer`; graphics and diagnostics improvements [1] |
| Institutional memory | `docs: solver design notes, benchmark history, repo restructure` [1] |

This phase is historically significant because it did not simply add features. It clarified **which parts of the old package should remain active, which should become archived, and which should be re-expressed through typed scenario objects and orchestration helpers**. That kind of narrowing is a sign of maturity. It means the repository is no longer trying to expose every layer of its own past as if all of it were equally current.

---

## 6. The architecture that exists today

The current repository surface, as reflected in the active module tree and workflow guidance, is organized around a staged scenario-first pipeline.[1] [2] The active modules now center on `scenarioTools.wl`, `examples.wl`, `unknownsTools.wl`, `systemTools.wl`, `solversTools.wl`, `graphicsTools.wl`, and `orchestrationTools.wl`, with the fast suite correspondingly focused on scenario, unknown, system, solver, graphics, and orchestration tests.[1] [2]

The key historical point is not merely that these file names changed. It is that the **repository’s public story changed**. The current package no longer presents itself primarily as a loose collection of solver scripts. It presents itself as a structured progression from typed scenario definition through unknown generation and system construction to solving, visualization, and orchestration.

```text
Scenario
  → Unknown bundle
    → Structural system
      → System solver / orchestration layer
        → Visualization and diagnostics
```

| Active layer | Current role |
|---|---|
| Scenario tools | Typed scenario construction, validation, topology derivation |
| Examples | Public scenario constructors and reusable example factories |
| Unknowns | Typed symbolic unknown bundles built from scenario topology |
| System | Modular structural-equation assembly from typed inputs |
| Solvers | Active symbolic-reduction and system-solving entrypoints |
| Graphics | Plotting and result-display helpers |
| Orchestration | Higher-level solve flows, validation, and workflow glue |

Historically, this architecture is the repository’s strongest synthesis so far. It preserves the long investment in symbolic modeling while making the active package surface smaller, more explicit, and easier to test.

---

## 7. Durable lessons from the history

Several lessons recur across the entire repository history.

First, the **Association-based modeling instinct was correct very early**. Even though names and wrappers changed, the repository consistently benefited from treating model inputs as structured, inspectable data rather than positional argument soup. That early decision made later typed wrappers feel natural rather than forced.

Second, **Wolfram Language context hygiene is not a cosmetic concern**. The March 2026 fixes demonstrated that context mistakes can remain latent for years and then break solver chains in non-obvious ways. In this repository, context discipline was not just style; it was correctness.

Third, **benchmarking changed the development culture**. Once performance history documents and benchmark scripts existed, solver work became legible as engineering. That is one of the clearest positive discontinuities in the repository timeline.[3] [4]

Fourth, **narrowing the active surface was a strength, not a retreat**. April 2026 succeeded in part because it accepted that not every historical module should remain equally public. Archiving and reducing the loader made the maintained package easier to reason about.

| Lesson | Evidence in the history |
|---|---|
| Good data models survive refactors | The project kept returning to structured scenario-style representations |
| Performance claims need records | Benchmark-history documents became a core institutional-memory tool [3] [4] |
| Process matters after the research phase | PRs, tests, and CI sharply improved repository clarity in 2026 [1] |
| Narrower public APIs are healthier | April 2026 deliberately reduced and clarified the active surface [1] [2] |

---

## 8. By the numbers

The current history is now large enough that approximate counts from older versions of this document are no longer sufficient. The exact snapshot on 2026-04-30 is as follows.[1]

### Commits by year

| Year | Commits |
|---|---:|
| 2020 | 160 |
| 2021 | 187 |
| 2022 | 89 |
| 2023 | 1 |
| 2024 | 18 |
| 2026 | 409 |
| **Total** | **864** |

### 2026 acceleration

| Period | Commits | Merged pull requests |
|---|---:|---:|
| March 2026 | 105 | 28 |
| April 2026 | 304 | 93 |

### Boundary commits

| Position | Commit |
|---|---|
| First | `09cc195` — `First commit` |
| Latest in this snapshot | `c594e8d` — `Merge pull request #175 from ricardodelimaribeiro/chore/ship-20260430-182409` |

These numbers matter because they quantify the scale of the repository’s recent transformation. The old story of MFGraphs as a slowly evolving research codebase is now incomplete. A very large share of its current identity was forged in 2026, and especially in April 2026.[1]

---

## 9. Where the history points next

The current history suggests a coherent forward direction. The active package surface is now organized enough that future work can focus less on rescuing old structure and more on extending the scenario-first workflow deliberately. The most likely next stage is not a return to sprawling solver sprawl. It is deeper consolidation around orchestration, diagnostics, and scenario-driven public entrypoints.[1] [2]

That direction would be consistent with the recent history. April 2026 invested heavily in typed scenario objects, modular builders, graphics helpers, benchmark history, and solver workflow cleanup. The next logical gains are therefore likely to come from making those layers more complete and more consistent with each other, rather than from reopening every archived branch of older architecture at once.

| Plausible next direction | Why the history supports it |
|---|---|
| Scenario-first public workflows | April 2026 concentrated effort on typed scenario construction and example factories [1] |
| Stronger orchestration and validation | Late-April history added orchestration and solver-diagnostics work [1] |
| Better benchmark continuity | The repository now preserves benchmark history as maintained documentation [3] [4] |
| Cleaner active/archive boundary | Recent history repeatedly archived or narrowed superseded surfaces [1] [2] |

---

## 10. Epilogue: the shape of this repository

MFGraphs still fits the classic pattern of research software: first prototype the mathematics, then stabilize just enough to produce results, and only later discover that software-engineering debt has accumulated everywhere at once. What makes this repository distinctive is how visible that transition is in the git history. The period from March through April 2026 did not merely “clean things up.” It changed the repository’s self-understanding.

The package today is best understood as a **benchmarked, scenario-first Wolfram toolkit with typed kernels, modular system construction, solver orchestration, and a documented performance memory**. That is a different repository from the one that began in May 2020, even though it descends directly from it.[1] [2] [3] [4]

## References

[1]: ./REPO_HISTORY_SNAPSHOT.md "Repository History Snapshot"
[2]: ../../CLAUDE.md "Canonical repository workflow and active package surface"
[3]: ./DNF_PERFORMANCE_HISTORY.md "DNF performance history"
[4]: ./PARALLEL_PERFORMANCE_HISTORY.md "Parallel performance history"
