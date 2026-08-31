# lbm-simulation-lib Documentation

lbm-simulation-lib provides two connected components:

- [lbm-sim](../lbm-sim/html/index.html): core lattice Boltzmann simulation framework.
- [lbm-gen](../lbm-gen/html/index.html): parser, validator, and model builder for simulation configuration files.

## Documentation Scope

This documentation set is generated for the full project and includes:

- File reference for project headers and source files included by Doxygen.
- Class and struct reference.
- Method and function reference.

Namespace documentation is intentionally minimal.

## Main Project Areas

- Core simulation (`src/include/lbm-sim`):
  - velocity sets, grids, solvers, collision operators, and problem definitions.
- Simulation generator (`src/include/lbm-gen`):
  - lexer/parser/AST, schema validation, deferred checks, and codegen model visitors.
- Tool entrypoint:
  - `src/generate_simulation.cpp` for config parsing and generation pipeline execution.

## Build and Generate Docs

From repository root:

```bash
cmake -S . -B build
cmake --build build --target lbm-full-docs 
```
