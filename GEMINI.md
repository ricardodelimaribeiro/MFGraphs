# GEMINI.md - MFGraphs Project Guidance

For comprehensive guidance on working with MFGraphs, see **[CLAUDE.md](CLAUDE.md)**.

## Quick Reference

**MFGraphs** is a Wolfram Language package for solving Mean Field Games on networks with congestion and switching costs.

### Load the package
```mathematica
Needs["MFGraphs`"]
(* or: Get["/path/to/MFGraphs/MFGraphs/MFGraphs.wl"] *)
```

### Quick start
```mathematica
Data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
d2e = DataToEquations[Data];
result = CriticalCongestionSolver[d2e];
result["AssoCritical"]
```

### Run tests
```bash
wolframscript -file Scripts/RunTests.wls fast  # ~27 min
wolframscript -file Scripts/RunTests.wls slow  # longer
```

### Run benchmarks
```bash
wolframscript -file Scripts/BenchmarkSuite.wls small
wolframscript -file Scripts/CompareDNF.wls --tag "optimization description"
```

## For Detailed Documentation

→ See **[CLAUDE.md](CLAUDE.md)** for:
- Environment setup and prerequisites
- Complete Quick Start with code examples
- Architecture overview (pipeline, module load order, solver chain)
- Configuration (all global parameters and their defaults)
- Debugging & Profiling workflows
- Performance optimization notes

→ See **[CONTRIBUTING.md](CONTRIBUTING.md)** for:
- PR workflow and code style
- Testing and benchmarking procedures
- Commit message conventions
- Documentation update guidelines

→ See **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** for:
- Common issues and solutions
- Solver timeouts and infeasibility
- Parallel kernel debugging
- Data format validation

## Key Files

| File | Purpose |
|------|---------|
| `MFGraphs/MFGraphs.wl` | Package loader |
| `MFGraphs/DataToEquations.wl` | Network → equations converter |
| `MFGraphs/DNFReduce.wl` | Boolean algebra solver |
| `MFGraphs/NonLinearSolver.wl` | Iterative solver (Hamiltonian) |
| `MFGraphs/Monotone.wl` | ODE-based gradient flow solver |
| `Scripts/RunTests.wls` | Test suite runner |
| `Scripts/BenchmarkSuite.wls` | Performance benchmarking |
| `Results/` | Benchmark outputs and performance history |
