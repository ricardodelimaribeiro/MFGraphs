# Contributing to MFGraphs

This project uses Wolfram Language for Mean Field Games on networks. We follow specific automation and documentation workflows to ensure theoretical alignment and code quality.

## Pull Request Workflow

### Before starting
1. Create a feature branch: `git checkout -b feature/descriptive-name` (never commit to `master`)
2. Reference any related issues in your branch name or PR description

### Code style and conventions
- Avoid single-letter `System` symbols (e.g., `K`, `D`, `C`). Use descriptive names or suffixes like `KM`
- For CLI scripts in `Scripts/`, use inline lint suppression to allow `Print[]`:
  ```wolfram
  (* :!CodeAnalysis::BeginBlock:: *)
  (* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
  ...
  (* :!CodeAnalysis::EndBlock:: *)
  ```

### Testing before submission
Run the test suite locally before opening a PR:

```bash
# Fast regression suite (9 files, ~27 min) — run before PR
wolframscript -file Scripts/RunTests.wls fast

# Full test suite (slower, more comprehensive)
wolframscript -file Scripts/RunTests.wls all
```

### Commit message format
- **First line:** Imperative, present tense, 50 characters max (e.g., "Fix DNFReduce exponential branching")
- **Blank line**, then **body** (if needed): Explain why the change was made, not what it does
- **Example:**
  ```
  Improve CLAUDE.md with better organization

  Reorganize sections for better flow, add Quick Start and Debugging
  sections for new developers. Consolidate test/benchmark info into
  tables for easier reference.
  ```

### PR description format
When opening a PR:
1. **Summary:** 2-3 bullet points of what changed
2. **Test plan:** How to validate the changes
3. **Performance impact** (if applicable): Did benchmarks improve or regress?
4. **Related issues:** Link any issues this closes or relates to

## Benchmarking and Performance Tracking

### Running benchmarks
Before and after solver/optimization changes, benchmark to measure impact:

```bash
# Establish baseline before your changes
wolframscript -file Scripts/BenchmarkSuite.wls small --tag "baseline"

# Make your changes, then test
wolframscript -file Scripts/BenchmarkSuite.wls small --tag "my optimization"

# Results automatically appended to history files
```

### Interpreting results
- Performance history files: `DNF_PERFORMANCE_HISTORY.md`, `PARALLEL_PERFORMANCE_HISTORY.md`
- Speedup > 1 means improvement; < 1 means regression
- Always benchmark at least `core` tier for solver changes
- `stress` tier helps validate on edge cases but may timeout (see `BENCHMARKS.md`)

### Profiling bottlenecks
After changes, use the profiling script to identify where time is spent:

```bash
wolframscript -file Scripts/BottleneckReport.wls
# Generates Results/bottleneck_report.md with timing breakdowns
```

## Documentation

### Updating code documentation
Public API docs are auto-generated from `::usage` strings. After adding/modifying functions:

```bash
wolframscript -file Scripts/GenerateDocs.wls
# Regenerates API_REFERENCE.md
```

### Updating developer guides
- **CLAUDE.md** is the canonical AI/developer workflow guide
- **GEMINI.md** is a synchronized companion reference for the same workflows
- **AGENTS.md** must remain a thin pointer (do not duplicate long-form guidance there)
- **BENCHMARKS.md** documents the profiling infrastructure — update when adding new benchmark cases
- **CONTRIBUTING.md** (this file) — update when workflow changes

### What NOT to manually edit
- `API_REFERENCE.md` — auto-generated; edit `::usage` strings in source instead
- `DNF_PERFORMANCE_HISTORY.md` / `PARALLEL_PERFORMANCE_HISTORY.md` — auto-appended by benchmark scripts

## Repository Structure
- `MFGraphs/`: Core package source
- `Scripts/`: Test runners, benchmarking, documentation generation
- `MFGraphs/Tests/`: MUnit test files (`.mt`)
- `Results/`: Benchmark outputs (CSV, JSON, markdown reports)
- `History/`: Development history and decision records
- `CLAUDE.md`: Canonical AI/developer workflow guidance
- `GEMINI.md`: Companion AI guidance synchronized with CLAUDE.md
- `AGENTS.md`: Pointer file that redirects to canonical guidance
- `BENCHMARKS.md`: Performance profiling methodology
- `TROUBLESHOOTING.md`: FAQ and common issues
