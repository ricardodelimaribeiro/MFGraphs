# Contributing to MFGraphs

This project uses Wolfram Language for Mean Field Games on networks. We follow specific automation and documentation workflows to ensure theoretical alignment and code quality.

## Development Workflow

### 1. Code Style and Linting
- Avoid using single-letter `System` symbols (e.g., `K`, `D`, `C`). Use descriptive names or suffixes like `KM`.
- For CLI scripts in `Scripts/`, use inline lint suppression blocks to allow `Print[]` statements:
  ```wolfram
  (* :!CodeAnalysis::BeginBlock:: *)
  (* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
  ...
  (* :!CodeAnalysis::EndBlock:: *)
  ```

### 2. Testing
Run the test suite before submitting changes. Tests are located in `MFGraphs/Tests/`.

```bash
# Fast tests (regression)
wolframscript -file Scripts/RunTests.wls fast

# Full suite (includes slow nonlinear cases)
wolframscript -file Scripts/RunTests.wls all
```

GitHub Actions will automatically run the `fast` suite on every Pull Request.

### 3. Benchmarking
To track performance changes, use the benchmark suite:

```bash
wolframscript -file Scripts/BenchmarkSuite.wls --tag "FeatureName"
```

Performance history is recorded in `PARALLEL_PERFORMANCE_HISTORY.md` and `DNF_PERFORMANCE_HISTORY.md`.

### 4. Documentation
Documentation is derived from `::usage` strings in the package files. After adding or modifying public symbols, regenerate the API reference:

```bash
wolframscript -file Scripts/GenerateDocs.wls
```

## Repository Structure
- `MFGraphs/`: Core package source.
- `Scripts/`: Automation and benchmarking logic.
- `Tests/`: MUnit test files (`.mt`).
- `Results/`: Output from benchmark runs.
