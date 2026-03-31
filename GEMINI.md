# GEMINI.md - MFGraphs Project

## Project Overview

This project contains `MFGraphs`, a Wolfram Language (Mathematica) package designed for solving Mean Field Games on networks. It specializes in scenarios involving congestion and switching costs. The core functionality of the package is to take a network's topological description (including vertices, edges, flows, and costs) and convert it into a system of equations. It then employs various solvers to find equilibrium flow distributions and value functions.

The project is structured as a Mathematica application, with source code primarily in `.wl` files. It includes several solvers: a "Critical Congestion" solver, a general "Non-linear" solver, and a "Monotone Operator" solver.

## Building and Running

The project doesn't have a traditional build process. It's used directly within a Mathematica environment.

### Running the code

1.  **Load the package:**
    Open a Mathematica session and load the main package file:
    ```mathematica
    Get["/path/to/MFGraphs/MFGraphs/MFGraphs.wl"]
    ```
    Or, if the package is in a location on Mathematica's `$Path`:
    ```mathematica
    Needs["MFGraphs`"]
    ```

2.  **Use the functions:**
    The `README.md` file provides extensive examples of how to define a network and use the solvers. A quick start example is:
    ```mathematica
    << MFGraphs`
    Data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
    d2e = DataToEquations[Data];
    result = CriticalCongestionSolver[d2e];
    result["AssoCritical"]
    ```

### Running Tests

The project includes a test suite that can be run from the command line using `wolframscript`. The test runner script is located in the `Scripts/` directory.

To run the tests, execute the following commands from the project root:

```bash
# Run the fast test suite
wolframscript -file Scripts/RunTests.wls fast

# Run the slow test suite
wolframscript -file Scripts/RunTests.wls slow

# Run all tests
wolframscript -file Scripts/RunTests.wls all
```

## Development Conventions

*   **Language:** The project is written in the Wolfram Language for Mathematica.
*   **File Structure:** The main source code is located in the `MFGraphs/` directory.
    *   `MFGraphs.wl`: Main package loader.
    *   `DataToEquations.wl`: Converts network data to equations.
    *   `NonLinearSolver.wl`: The main iterative solver.
    *   `Monotone.wl`: An ODE-based solver.
    *   `DNFReduce.wl`: Handles boolean algebra for Disjunctive Normal Form reduction.
    *   `Examples/ExamplesData.wl`: Contains built-in example networks.
    *   `Tests/`: Contains test files (`.mt`).
*   **Testing:** Tests are written in `.mt` files and are executed via the `RunTests.wls` script.
*   **Verbose Logging:** The package uses a global variable `$MFGraphsVerbose` to control the verbosity of the output. It is `True` by default.
