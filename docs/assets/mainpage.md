# lbm-simulation-lib {#mainpage}

Documentation portal for **lbm-simulation-lib**, a header-only C++17 Lattice
Boltzmann library with an OpenMP and a CUDA backend.

## Components

- [lbm-sim](../lbm-sim/html/index.html): the simulation core -- velocity sets,
  lattice and grid, collision operators, boundary conditions, solvers, output
  listeners and error analysis. Its landing page carries the usage walkthrough,
  and its *Architecture* page the design and the extension points.

The portal itself contains no API reference: it only links the component
manuals, which are generated from the annotated headers under
`include/lbm-sim/`.

## Source layout

| Path | Contents |
|------|----------|
| `include/lbm-sim/` | the header-only library (everything is documented in the `lbm-sim` manual) |
| `include/lbm/`, `src/lbm/` | the Quill logging wrapper -- the project's only translation unit |
| `simulations/openmp/`, `simulations/cuda/` | example mains: Couette, Poiseuille, lid-driven cavity, one executable per file |
| `benchmarks/ghia/` | Ghia et al. (1982) reference tables for the lid-driven cavity |
| `scripts/` | Python post-processing: animations, velocity profiles, scaling plots |
| `docs/report/` | project report: theory, benchmarks, error results, profiling, flow animations |

## Build and generate the docs

From the repository root:

```bash
cmake -S . -B build
cmake --build build --target lbm-docs-full
```

| Target | Output | Contents |
|--------|--------|----------|
| `lbm-sim-docs` | `docs/doxygen/lbm-sim/` | the library manual: landing page, architecture page, and the header reference |
| `lbm-docs-full` | `docs/doxygen/` | this portal; depends on `lbm-sim-docs`, so it builds both |

Both targets are configured only when Doxygen is found at CMake configure time.
`*.cuh` headers -- the CUDA solver and its helpers -- are added to the input set
only when the project is configured with `-DLBM_ENABLE_CUDA=ON`.
