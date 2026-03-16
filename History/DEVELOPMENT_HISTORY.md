# MFGraphs: A Development History

*A chronicle of 521 commits across six years — from a research prototype to a production-grade Mean Field Games solver.*

---

## 1. What This Package Does

MFGraphs is a Wolfram Language package for solving **Mean Field Games on networks** with congestion and switching costs. Given a network topology — vertices, directed edges, entry/exit flows, and optional switching costs — it constructs a system of equations encoding equilibrium conditions, then solves for the flow distributions and value functions that describe how rational agents distribute themselves across the network.

The mathematical pipeline is deceptively simple in concept:

```
Network data → Equations → Solver → Equilibrium flows
```

Getting this right took six years, three solver implementations, one breakthrough algorithm, and more than a few "heart attacks."

---

## 2. The Research Phase (May – November 2020)

### The first week

The repository was born on **May 6, 2020** with a `First commit` — 12 files laying the groundwork for a Wolfram Language package. Two days later, eight commits landed in a single day, all with the same message: **"Association"**. This was the foundational data model being worked out in real-time: how should a network be represented?

The answer — a Mathematica `Association` with keys like `"Vertices List"`, `"Adjacency Matrix"`, `"Entrance Vertices and Flows"` — proved to be a durable choice. It survived all six years essentially unchanged.

### The emotional commits

The commit messages from this period read like a research diary:

- *"Having heart attacks and recovering quickly!"* (May 17)
- *"DONE"* (May 12)
- *"problem with jays"* (May 18) — the flow variables `j[...]` were misbehaving
- *"what is this?"* (June 22)
- *"getting to the final idea. found a way to reduce the equations more efficiently."* (July 15)

This was clearly a solo researcher working through hard mathematics, committing progress as it happened, sometimes multiple times per day, sometimes disappearing for weeks.

### The solver comes alive

**August 31, 2020** marked the first major milestone:

> *"Non-linear solver is working! Maybe we need to look carefully into the NonNegative equations so that we avoid Reduce."*

The qualifying remark — "maybe we need to look carefully" — was prophetic. The interaction between non-negative constraints and the `Reduce` function would remain a central tension in the codebase for years to come.

### ZAnd is born

**November 10, 2020** brought a pivotal commit:

> *"New ZAnd, ReplaceSolution, and NewReduce functions. The final associations are now simplified (just elementary numeric computations)"*

`ZAnd` — a function that recursively processes conjunctions and disjunctions of equations — would become the package's critical bottleneck and, eventually, its most optimized component. The name was cryptic (a portmanteau of "Z" for the zero-flow case and "And" for logical conjunction), and it would carry that name for almost six years before being renamed to `DNFReduce`.

### Switching costs enter the picture

The period from **November 16–26** saw the introduction of switching costs — the mechanism that models the price agents pay when transitioning between edges at a vertex. This was mathematically elegant but computationally devastating: switching costs generate `Or` expressions over all feasible switching patterns, and a network with N switching costs can produce up to 2^N branches.

> *"Still trying to understand the switching costs problem"* (November 25)

This combinatorial explosion would remain unaddressed for years.

### Assessment of the early phase

**Good decisions:**
- **The Association data model.** Choosing an Association with descriptive string keys (`"Vertices List"`, `"Adjacency Matrix"`, etc.) made the data self-documenting and extensible. Six years later, the format is unchanged — new keys were added, but existing ones never needed modification.
- **Symbolic parameters.** Using symbolic placeholders (`I1`, `U1`, `S1`) that get substituted with `/.` rules was a natural fit for Mathematica and allowed the same network definition to serve both symbolic analysis and numerical solving.
- **Pipeline thinking.** Even in the earliest code, there was a clear separation between "construct the equations" and "solve the equations." This would eventually crystallize into the three-stage pipeline.

**Bad decisions:**
- **Notebook-driven development.** All experimentation happened in `.nb` notebook files that mixed code, output, and commentary. These were committed to the repo, creating a growing collection of files that were impossible to diff, review, or maintain.
- **No tests.** The first unit test wouldn't appear until September 2021 — 16 months after the first commit. Correctness was verified by running notebooks and checking output visually.
- **Direct commits to master.** Every change went straight to the main branch. No branches, no pull requests, no code review.
- **Meaningless commit messages.** Eight commits titled "Association." Twenty-plus commits titled "Update D2E2.m." The git history became a timeline of activity rather than a record of decisions.

---

## 3. The Long Middle (2021 – 2024)

### Broadening the scope (January 2021)

January 2021 was an active month that expanded the package's ambitions:

- **Braess paradox networks** were implemented — a classic game theory example where adding a road to a network can make everyone worse off. Multiple commits over January 19–27 worked through the mathematics:
  > *"Now we can see the paradox!"* (January 19)
  > *"tried to solve for the braess paradox directly!"* (January 27)

- **Jamarat pilgrimage networks** appeared — a real-world application modeling crowd flow at the Jamarat Bridge in Mecca. This was clearly driven by collaboration:
  > *"running examples with Fatimah"* (January 5)

- **Plotting functions** were added for mass densities and value functions on edges.

### The solver odyssey (2021)

The non-linear solver went through multiple iterations throughout 2021:

- *"FixedStepX does not halt. But it is not reaching a good solution."* (October 25, 2020)
- *"Developed a new version of Eq(ualities)Eliminator."* (September 1, 2020)
- *"working on non-linear solver"* (October 25–26, 2021)
- *"running faster"* (November 1, 2021)

The pattern was clear: get a solver working on simple cases, discover it fails on complex ones, rewrite, repeat. The commits show a researcher systematically pushing the boundary of what the solver could handle, from 1-edge networks to Y-networks to 9-vertex grids.

### First tests (September 2021)

**September 16, 2021** marked a quiet milestone:

> *"developed the first test!"*

Followed three days later by:

> *"new unit tests"*
> *"include two functioning tests in the suit"*

The tests used Mathematica's `.mt` format and focused on critical congestion solutions. This was a step in the right direction, but the test infrastructure remained minimal — no test runner, no automation, no CI.

### The quiet years (2022 – 2023)

Activity slowed considerably. Commits became sparse, often just "Update D2E2.m" or "Update NonLinearSolver.m" with weeks or months between them. The package was being used for research but not actively developed as software.

Notable events:
- **April 2022**: D2E2 (Data-to-Equations v2) was being actively improved: *"habemus critical solver!"*
- **May 2022**: A new investigation into ZAnd: *"thinking about a different implementation of ZAnd"* — but nothing came of it yet.
- **June 2022**: *"ZAnd is 'close' to optimal."* — This assessment would prove overly optimistic.

### The Monotone solver (June 2024)

After a 16-month gap (the longest in the project's history), **June 4, 2024** brought a significant addition: `Monotone.m`, an entirely new solver based on ODE gradient flow on the Kirchhoff matrix.

This represented a fundamentally different mathematical approach — instead of fixed-point iteration, the Monotone solver formulated the problem as an ODE and used `NDSolve` to follow the gradient flow to equilibrium. It was also the beginning of the D2E2 consolidation: duplicate code from `D2E2-currents.m` was merged into the main `D2E2.m`.

Then silence again — 9 months — until March 2026.

### Assessment of the middle phase

**Good decisions:**
- **Multiple solver strategies.** Having both a fixed-point iterative solver (NonLinearSolver) and an ODE-based solver (Monotone) gave the package resilience — different problem types could favor different approaches.
- **Real-world applications.** The Jamarat and Braess examples grounded the mathematics in practical problems, which helped identify edge cases and drove solver improvements.
- **Starting to write tests.** Even minimal tests created a foundation that could be built upon.

**Bad decisions:**
- **Code duplication across files.** `ZAnd` and `ReZAnd` were implemented in multiple files (`ZAnd.m`, `D2E2.m`, `D2E2-currents.m`) with subtle differences. When one was improved, the others might not be.
- **No performance measurement.** For a package whose primary bottleneck (ZAnd/DNFReduce) could vary from milliseconds to minutes depending on the problem, there was no systematic way to measure or track performance.
- **Growing notebook debt.** By the end of this period, the repository contained 49+ example notebooks. Many were outdated, referenced old function names, or duplicated each other. They were effectively dead code but took up space and created confusion about what was current.
- **No package hygiene.** Symbol contexts (`Public` vs `Private`), variable scoping, and naming conventions were inconsistent. Variables leaked between modules. These would become serious bugs later.

---

## 4. The Great Professionalization (March 2026)

### The transformation

On **March 1, 2026**, the repository underwent a transformation so dramatic that it deserves its own chapter. In 16 days, 18 pull requests were merged — more structured development than the previous five and a half years combined.

This wasn't gradual improvement. It was a deliberate, intensive effort to convert a research codebase into a properly engineered package.

### PR #1–2: Structure and documentation (March 1)

The first pull request refactored the package structure, fixed bugs across all modules, introduced `$MFGraphsVerbose` for output control, and cleaned up global variable leaks. The second added a comprehensive README.

These set the tone: the codebase was being treated as software, not just a collection of scripts.

### PR #4: The benchmarking suite (March 3)

This was the infrastructure that made everything else possible. A full benchmarking suite with:

- **Four test tiers** by complexity: small (60s timeout), medium (300s), large (900s), vlarge (1800s)
- **31+ test cases** covering the full range of network types
- **SafeExecute wrappers** with timeout handling
- **Profiling instrumentation** that could wrap any function non-invasively

With this in place, every subsequent change could be measured against a baseline. This is the commit that turned "I think this is faster" into "this is 47x faster."

The same commit also introduced the first performance optimization: **Solve memoization** via `CachedSolve` with SHA-256 hashing. Identical equalities — common when processing Or-branches with switching costs — would now be solved only once.

### PR #7: The DNFReduce breakthrough (March 4)

This was the single most impactful change in the project's history. The cryptic `ZAnd` was renamed to `DNFReduce` (Disjunctive Normal Form Reduce), making its purpose self-documenting. But the rename was the least important part.

Two algorithmic optimizations delivered dramatic speedups:

1. **Reduce memoization** (`CachedReduce`): The `Reduce[fst, Reals]` call in the Or case was cached by SHA-256 hash. The same sub-expression across different Or-branches would be reduced only once.

2. **Short-circuit Or evaluation**: Before evaluating the second Or-branch, check if the first result already equals `xp` (the accumulated expression). If so, the second branch cannot contribute anything new — skip it entirely.

The results, measured by the newly-created `CompareDNF.wls` script:

| Case | Problem | Speedup | Output compression |
|------|---------|---------|-------------------|
| 12 | Attraction, 5 free vars | **47.8x** | 65,536 disjuncts → 1 |
| 11 | Attraction, 16 switching costs | **10.8x** | 3,359,232 → 2,977 |
| Braess | Braess congestion | 0.86x (slower) | 147,456 → 179 (823x smaller) |

The key insight: for Case 12, the short-circuit fires on the very first Or-branch. One solve is enough; all 65,535 remaining branches are skipped. `BooleanConvert` enumerates all 65,536.

For smaller networks where DNFReduce was actually slower on wall-clock time, the output was dramatically more compact — 3 to 1,024 times fewer disjuncts — which made downstream solver steps (TripleClean, MFGSystemSolver) correspondingly faster.

### PR #8–9: Performance tracking and more optimizations (March 4–10)

`DNF_PERFORMANCE_HISTORY.md` was created — a document that records every optimization milestone with measured data, code diffs, and interpretation. This is one of the project's best engineering practices: a living record of *why* changes were made and *what effect they had*.

PR #9 added two more optimizations:
- **And-Or early exit**: When distributing over Or-branches in the And case, stop immediately if any branch fully resolves the problem.
- **TripleStep caching**: The `TripleStep` function in DataToEquations now used `CachedSolve` instead of bare `Solve`.

All 31 valid benchmark cases continued to pass with no regressions.

### PR #10–12: The great cleanup (March 15)

Three pull requests on a single day:

- **PR #10**: Deleted 49 legacy notebooks and obsolete code — 574 KB removed. The example notebooks that had accumulated since 2020 were replaced by `GetExampleData[key]`, which loaded the same networks programmatically.
- **PR #11**: Standardized function naming to PascalCase (`DataToEquations`, `CriticalCongestionSolver`, `GetExampleData`). The old name `D2E2` was gone.
- **PR #12**: Renamed all package files from `.m` to `.wl` — the modern Wolfram Language extension.

### PR #13–17: The context hygiene crisis (March 15–16)

This was the most instructive episode in the project's history.

**PR #13** (March 15) discovered that `Begin["`Private`"]` statements in several modules were placed *after* symbols had already been used. This meant internal symbols leaked into the public `MFGraphs`` context, while other symbols ended up in `MFGraphs`Private`` where user substitution rules couldn't reach them.

The fix seemed comprehensive. It wasn't.

**PR #15** (March 16, next morning) found that MonotoneSolver had a cache part-assignment bug: `cache[[1]] = ...` doesn't work on function arguments in Mathematica. The cache was silently failing.

**PR #17** (March 16, afternoon) uncovered a deeper issue: `Clear[H, Cost]` in NonLinearSolver had created a `MFGraphs`Cost` symbol, but DataToEquations used `MFGraphs`Private`Cost`. When the NonLinearSolver evaluated `Cost[50, {2,3}]`, it stayed unevaluated — the wrong symbol — and crashed `TripleStep`.

**Three commits in 24 hours** to fix context issues that had been latent for years. The bugs were always there; they just hadn't been triggered because no one had run the full solver chain with the right combination of test cases.

This episode illustrates a fundamental truth about Wolfram Language development: **symbol context management is the language's most treacherous feature**, and getting it wrong produces failures that are silent, delayed, and cascading.

### PR #16: Hessian Riemannian Flow (March 16)

Amid the bug fixes, a genuine mathematical advance: the MonotoneSolver was upgraded from simple gradient projection to **Hessian Riemannian Flow**, incorporating full Hessian information for more robust convergence. This came with proper monotonicity regression tests.

### PR #18: Solver contracts and test runner (March 16)

The final PR of this period normalized the return format across all three solvers via a `"ReturnShape" -> "Standard"` option, created a comprehensive test runner (`RunTests.wls` with fast/slow/all tiers), and added solver contract tests verifying return format compliance.

### Assessment of the professionalization phase

**Good decisions:**
- **Benchmark-first approach.** Building the measurement infrastructure before optimizing meant every change had quantified impact. This prevented both premature optimization and the "I think it's faster" trap.
- **DNF_PERFORMANCE_HISTORY.md.** A document that records optimization rationale, code changes, and measured results is rare and valuable. It serves as both documentation and institutional memory.
- **Backward-compatible aliases.** `ZAnd = DNFReduce` and `DataG := GetExampleData` meant existing notebooks didn't break during the rename.
- **Tiered test cases.** Separating tests into small/medium/large/vlarge tiers with different timeouts acknowledged that solver performance varies by orders of magnitude across problem sizes.
- **Standardized solver contracts.** Having all three solvers return results in a consistent format (with opt-in normalization) made them interchangeable downstream.

**Bad decisions (or: lessons from this phase):**
- **Context hygiene was left too late.** The `Begin["`Private`"]` ordering issue affected every module and had been present since the earliest code. Fixing it required three emergency PRs because each fix revealed the next layer of context confusion.
- **No integration tests before the cleanup.** The context bugs only surfaced because the professionalization phase involved running the full solver chain on all test cases — something that apparently hadn't been done systematically before.
- **Part-assignment on function arguments.** The MonotoneSolver cache tried `cache[[1]] = result` where `cache` was a function parameter. This is a Mathematica-specific pitfall: `Part` assignment doesn't work on function arguments the way it works on module-local variables. The fix required inlining the cache variable.

---

## 5. The Architecture That Emerged

After six years of evolution, the package settled on a clean three-stage pipeline:

```
DataToEquations[data]
  → CriticalCongestionSolver[d2e]
    → NonLinearSolver[critical]   or   MonotoneSolverFromData[data]
```

Each stage has a well-defined input and output:

1. **DataToEquations** takes a network Association and produces a system decomposed into equalities (EE), inequalities (NN), and disjunctions (OR).
2. **CriticalCongestionSolver** solves the zero-flow equilibrium using `DNFReduce` for the Boolean algebra. It returns the solution under `"AssoCritical"`.
3. **NonLinearSolver** or **MonotoneSolver** takes the critical solution as a starting point and solves the full non-linear problem iteratively.

The inner solver pipeline within `DataToEquations.wl` is where most of the computational work happens:

```
CriticalCongestionSolver
  → MFGPreprocessing → TripleClean (fixed-point: solve equalities, substitute, repeat)
  → MFGSystemSolver → TripleClean → DNFSolveStep (DNFReduce on Or-branches)
```

This architecture wasn't designed up front — it emerged through iteration. The early code had everything in one file; the middle period saw functions migrate between files; the final phase consolidated everything into clean, single-responsibility modules.

---

## 6. Good Decisions: A Retrospective

### 1. The Association data model (Day 1)

Using a Mathematica Association with human-readable string keys was the right call from the start. It made the data self-documenting, extensible (new keys could be added without breaking old code), and natural for Mathematica's functional style. The same format serves symbolic analysis, numerical solving, and visualization.

### 2. Symbolic parameters with late binding

Defining networks with symbolic placeholders (`I1`, `U1`, `S1`) that get substituted before solving was a natural fit for Mathematica and enabled parametric analysis. The same network definition could be solved for different parameter values without reconstruction.

### 3. The pipeline architecture

Separating "construct equations" from "solve equations" from "iterate non-linearly" allowed each stage to be optimized independently. When `DNFReduce` got its 47x speedup, nothing downstream needed to change.

### 4. Multiple solver approaches

Having three independent solvers (Critical, NonLinear, Monotone) gave the package resilience. The Monotone solver uses a completely different mathematical framework (ODE gradient flow) than the NonLinear solver (fixed-point iteration), so problems that are pathological for one may be tractable for the other.

### 5. SHA-256 memoization

Using cryptographic hashing to detect duplicate Solve/Reduce inputs was elegant and effective. It avoided the memory bloat of storing the full expression as a key while providing reliable cache hits. The decision to clear the cache between problem instances (`ClearSolveCache[]`) prevented stale results.

### 6. Measuring before optimizing

The benchmarking suite (March 2026) was built *before* the major optimizations. This meant every change had a measured baseline and a measured result. The `CompareDNF.wls` script automated this, and `DNF_PERFORMANCE_HISTORY.md` preserved the context.

### 7. Backward-compatible aliases

When `ZAnd` was renamed to `DNFReduce`, a single line — `ZAnd = DNFReduce` — preserved compatibility with existing notebooks. Similarly for `DataG := GetExampleData`. This is the right way to rename: deprecate gracefully rather than break everything.

---

## 7. Bad Decisions: Lessons Learned

### 1. Six years without branches or pull requests

From May 2020 to February 2026, every commit went directly to `master`. No branches, no pull requests, no code review. This meant:
- No opportunity for review before changes landed
- No way to experiment without risking the main branch
- No separation between "work in progress" and "known good"
- The git history became a stream of consciousness rather than a record of deliberate changes

The switch to branch-and-PR workflow in March 2026 immediately improved code quality — not because the code was different, but because the process enforced deliberation.

### 2. Notebook-driven development

Mathematica notebooks (`.nb` files) are wonderful for exploration but terrible for version control. They contain binary formatting data, they can't be meaningfully diffed, and they mix code, output, and commentary in ways that make it impossible to tell what's "current."

By 2026, the repository contained 49 example notebooks, many of which were outdated or duplicated each other. Deleting them (PR #10) removed 574 KB of dead weight. The replacement — `GetExampleData[key]` — was cleaner, more maintainable, and testable.

**The lesson:** Use notebooks for exploration, but extract production code into `.wl` files and test cases into `.mt` files. Commit notebooks only when they serve as documentation, not as the source of truth.

### 3. No caching for six years

The `Solve` and `Reduce` calls in `ZAnd` were the dominant bottleneck from day one. The switching cost combinatorics (up to 2^16 branches) meant identical sub-expressions were being solved over and over. Yet caching wasn't added until March 2026.

The original ZAnd code:
```mathematica
With[{fsol = First@Solve@newfst},    (* Solve uncached — called thousands of times *)
    ZAnd[ReplaceSolution[xp, fsol] && fst, ReplaceSolution[rst, fsol]]
]
```

This is a common research-code pattern: correctness first, performance later. But "later" was six years, and the performance impact was 10–47x.

**The lesson:** For symbolic computation, memoization should be a reflex, not an optimization. If a pure function (same input → same output) is called repeatedly, cache it.

### 4. Ignoring symbol context hygiene

Wolfram Language's context system (`MFGraphs``, `MFGraphs`Private``, `Global``) determines which symbols are visible where. Getting this wrong is one of the most common sources of bugs in Mathematica packages, and MFGraphs got it wrong in every module.

The specific pattern:
```mathematica
Begin["`Private`"]
(* ... many lines of code ... *)
somePublicFunction[args_] := ...  (* Oops: this is now Private *)
```

When `Begin["`Private`"]` appears before symbol declarations, those symbols end up in the wrong context. User rules like `/.{I1 -> 100}` target `Global`I1` but the package expects `MFGraphs`I1`. The substitution silently fails, and the solver produces wrong results or crashes downstream.

It took three emergency PRs in 24 hours to fix this across all five modules. The bugs had been latent for years — they only manifested when the full solver chain was run systematically for the first time.

**The lesson:** In Wolfram Language packages, always declare public symbols (with `::usage` messages) *before* `Begin["`Private`"]`. Test symbol visibility explicitly. And test the full pipeline, not just individual functions.

### 5. Code duplication across files

`ZAnd`, `ReZAnd`, and `ReplaceSolution` were implemented in at least three different files (`ZAnd.m`, `D2E2.m`, `D2E2-currents.m`) with subtle differences. When one was improved, the others weren't necessarily updated.

This is a natural consequence of research-style development: you copy a function to a new file, modify it, and forget to reconcile. The fix was simple but came late — consolidating everything into `DNFReduce.wl` in March 2026.

### 6. Vague commit messages

The git log contains at least 20 instances of "Update D2E2.m", 8 instances of "Association", and numerous others like "running examples", "running notebooks!", "minor changes", and "running." These tell you *when* something changed but not *what* or *why*.

Compare with the March 2026 messages: *"Performance: sequential And-Or early exit, NonLinear tolerance, TripleStep caching"* or *"Fix package context hygiene, solver bugs, and add infeasibility status."* The later messages are self-contained summaries that make the git history useful as documentation.

### 7. Global variable leaks

Variables like `style`, `A`, and `bool` in NonLinearSolver.m were defined at the top level rather than inside `Module` or `Block` scopes. This meant:
- Loading the package could silently overwrite user variables
- Multiple solver invocations could interfere with each other
- The behavior depended on evaluation order rather than explicit inputs

This is Mathematica's version of the global variable problem, and it's especially insidious because Mathematica's dynamic scoping means a global `A` can affect expressions evaluated inside unrelated functions.

---

## 8. By the Numbers

### Timeline

| Period | Commits | Key milestone |
|--------|---------|--------------|
| May–Nov 2020 | ~120 | First solver, ZAnd, switching costs |
| 2021 | ~150 | Braess paradox, Jamarat, first tests |
| 2022 | ~50 | D2E2 refinement, ZAnd investigation |
| 2023 | ~5 | Minimal activity |
| 2024 | ~15 | Monotone solver, D2E2 consolidation |
| Mar 2026 | ~180 | 18 PRs, full professionalization |
| **Total** | **~521** | |

### File evolution

| File (current name) | Born | Original name | Major rewrites |
|---------------------|------|---------------|----------------|
| DNFReduce.wl | Nov 2020 | ZAnd.m | Mar 2026 (memoization, rename) |
| DataToEquations.wl | May 2020 | DataToEquations.m → D2E2.m | Jun 2024 (consolidation), Mar 2026 (cleanup) |
| NonLinearSolver.wl | Aug 2020 | IterationFunction.m → NonLinearSolver.m | Mar 2026 (tolerance, caching) |
| Monotone.wl | Jun 2024 | Monotone.m | Mar 2026 (Hessian flow, context fixes) |
| ExamplesData.wl | May 2020 | ExamplesParameters.m → ExamplesData.m | Mar 2026 (34 cases) |

### Pull requests (March 2026)

| PR | Focus | Impact |
|----|-------|--------|
| #1–2 | Structure + README | Foundation |
| #4 | Benchmarking suite | Measurement infrastructure |
| #7 | DNFReduce + memoization | 47x speedup |
| #8–9 | Performance tracking + optimizations | And-Or early exit, TripleStep caching |
| #10–12 | Cleanup + naming + .wl extension | 574 KB removed, modern conventions |
| #13–17 | Context hygiene + solver fixes | 3 emergency fixes in 24 hours |
| #18 | Solver contracts + test runner | Standardized interfaces |

---

## 9. Epilogue: The Shape of Research Software

MFGraphs follows a pattern common to academic software projects:

1. **Rapid prototyping** — correctness is the only goal; code quality is irrelevant
2. **Long stabilization** — the code works well enough for papers; maintenance is minimal
3. **Professionalization crisis** — the code needs to be shared, tested, or extended; years of technical debt surface simultaneously

The transition from phase 2 to phase 3 is always painful. In MFGraphs' case, it was compressed into 16 days of intensive work. The context hygiene bugs that took three emergency PRs to fix had been present since 2020. The missing memoization that delivered 47x speedup had been *discussed* in commit messages since 2022 ("ZAnd is 'close' to optimal").

The lesson isn't "write perfect code from the start." Research code *should* prioritize mathematical correctness over engineering practices. The lesson is: **when you decide to professionalize, invest in measurement first.** The benchmarking suite (PR #4) made everything that followed possible — without it, the DNFReduce optimization would have been guesswork, the context fixes would have had no regression detection, and the cleanup would have been a leap of faith.

The package today is in a strong position: three solver implementations, 34 built-in test cases across four complexity tiers, comprehensive benchmarking infrastructure, performance history tracking, and standardized interfaces. It took six years of mathematical research and 16 days of engineering to get there.
