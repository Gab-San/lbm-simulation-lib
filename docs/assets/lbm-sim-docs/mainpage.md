# LBM SIM

`lbm-sim` is the core of **lbm-simulation-lib**: a header-only C++17
implementation of the Lattice Boltzmann Method (LBM) for incompressible flow,
with an OpenMP and a CUDA backend behind a single solver interface.

Everything under `include/lbm-sim/` is header-only and consumed through the
`lbm-sim` INTERFACE target; the only translation unit in the project is
`src/lbm/logging.cpp`, which backs the small Quill wrapper in
`include/lbm/logging.hpp`.

The design rests on three ideas, described in detail in @subpage architecture:

- **compile-time polymorphism** -- dimension, velocity set, collision operator
  and execution backend are template parameters, so the inner loop contains no
  virtual call and no run-time branch on the model;
- **one lattice, two population buffers** -- streaming and collision are fused
  into a single pass that reads one buffer and writes the other;
- **observer-based output** -- solvers push raw byte chunks to listeners and
  know nothing about file formats.

## What is implemented today

The template parameters are written for the general case, but the
implementation is currently 2D-only: the missing pieces are `static_assert`s,
not silent fallbacks, so an unsupported instantiation fails at compile time.

| Area | Available | Reserved / not implemented |
|------|-----------|----------------------------|
| Dimension | `dim == 2` | `dim == 3`: `Grid<3>::scalar_index()`, `Grid<3>::field_index()`, `CollisionParams`, `ErrorEvaluator` and `utils::ops::dot()` all `static_assert` |
| Velocity sets | `lbm::D2Q9` | `D3Q27` is present but commented out in `core/velocity-sets.hpp` |
| Collision operators | `CollisionModel::BGK`, `CollisionModel::TRT` | `CollisionModel::MRT`: the enumerator exists, `CollisionStrategy::apply()` `static_assert`s |
| Backends | `MPISolver2D` (OpenMP), `CUDASolver2D` (CUDA, opt-in) | -- |
| Boundary conditions | `BB_RIGID_WALL`, `BB_MOVING_WALL`, `PERIODIC`, `PRESSURE_PERIODIC_INLET`, `PRESSURE_PERIODIC_OUTLET` | conditions are resolved only on links leaving the domain, so obstacles **inside** the fluid are not handled yet |
| Geometry | `Segment`, `Circle`, `Parallelogram` rasterized into a per-node mask | `Circle` in 3D throws; no airfoil or shell shapes |
| Output | `AsyncBinaryWriter` (frames) and `LBMSimulation::output()` (centerline profile) | no VTK writer |
| Analysis | `L1`, `L2`, `L2^2`, `Linfty` against an `analysis::Function<2>`; Ghia et al. (1982) benchmark | `compute_ghia_error()` is 2D only by construction |

@note `MPISolver2D` is the **OpenMP** solver. The name is a leftover from an
earlier MPI prototype (the class carries a `FIXME: CHANGE NAME TO
OpenMPSolver2D`); it derives from `SolverBase2D<cm_t,
ExecutionBackend::OPEN_MP>` and contains no MPI code.

## A complete run

Condensed from `simulations/openmp/couette_flow_2d_bgk.cpp`. The six steps are
the same for every simulation in the project:

```cpp
using namespace lbm;
using types::Coordinate;
using utils::Vector;

constexpr unsigned short int DIM = 2;
constexpr auto CollisionType = CollisionModel::BGK;
using Simulation = LBMSimulation<DIM, D2Q9, CollisionType>;

const types::DimPoint<2> grid_size{129, 129};
const Vector<double, 2> init_vel{0.1, 0.0};

// --- 1. LOGGING ----------------------------------------------------------
logging::setup_quill();
quill::Logger *main_logger = logging::create_or_get_logger("main");

// --- 2. GEOMETRY ---------------------------------------------------------
// Each CollisionArea is one obstacle id; the id indexes obst_type_map,
// which says what boundary condition its nodes carry.
const Coordinate<2> A(0, 0), B(0, 128), C(128, 128), D(128, 0);

const std::vector<CollisionDetection::CollisionArea<DIM>> obstacles{
    CollisionDetection::CollisionArea(A, {CollisionDetection::Segment(A, D)}),
    CollisionDetection::CollisionArea(A, {CollisionDetection::Segment(B, C)}),
    CollisionDetection::CollisionArea(
        A, {CollisionDetection::Segment(A + Vector<int, DIM>(0, 1),
                                        B - Vector<int, DIM>(0, 1)),
            CollisionDetection::Segment(C - Vector<int, DIM>(0, 1),
                                        D + Vector<int, DIM>(0, 1))}),
};

const std::unordered_map<unsigned int, uint8_t> obst_type_map{
    {0, Solid::BB_RIGID_WALL},  // bottom wall
    {1, Solid::BB_MOVING_WALL}, // lid, driven at init_vel
    {2, Solid::PERIODIC},       // left and right sides
};

// --- 3. RASTERIZE THE GEOMETRY INTO A PER-NODE MASK ----------------------
types::boundary_mask_t boundary_mask =
    Solid::compute_boundary_mask<DIM>(obst_type_map, obstacles, grid_size);

// --- 4. WIRE THE OUTPUT --------------------------------------------------
// The writer must be attached to BOTH observables: the simulation emits the
// grid header, the solver emits the velocity-norm frames.
auto writer = std::make_shared<AsyncBinaryWriter>("out/norms.bin");

Simulation simulation(grid_size, boundary_mask,
                      CollisionParams<DIM, CollisionType>(
                          /*Re*/ 100.0, grid_size, init_vel));
simulation.attachListener(writer);

MPISolver2D<CollisionType> solver(/*iters*/ 100000, /*frames*/ 300);
solver.attachListener(writer);

// --- 5. SOLVE ------------------------------------------------------------
const LidCavity2D problem;
simulation.solve(solver, problem);

simulation.output("out/profile.bin",
                  functional::extract_dx_profile_along_y_center);

// --- 6. MEASURE THE ERROR ------------------------------------------------
const double H = static_cast<double>(grid_size.y - 1);
const analysis::CouetteSolution2D exact(H, init_vel.dx);
const double err_l2 =
    simulation.compute_error(analysis::NormType::L2, exact);

simulation.detachListener(writer);
solver.detachListener(writer);
```

Two things about this snippet are worth stating explicitly, because they are
easy to get wrong:

- **The lid velocity is not an initial condition.** `LidCavity2D::init()` is
  currently not called by `LBMSimulation::solve()` (a `FIXME` sits in its
  place), so every run starts from `u = 0`, `rho = 1`. The wall velocity enters
  the simulation only through the `BB_MOVING_WALL` momentum correction, which
  reads `CollisionParams::init_vel`.
- **`init_vel` has three jobs at once**: it is the reference velocity that sets
  the viscosity (`nu = init_vel.dx * ny / Re`), the velocity of every moving
  wall, and the value written into the profile header used to normalize the
  plots. Changing it changes the physics, not just a label.

## Header map

| Header | What it holds |
|--------|---------------|
| `lbm-simulation.hpp` | `lbm::LBMSimulation` -- owns the lattice, drives a solver, computes errors, writes profiles |
| `lattice.hpp` | `lbm::Lattice` -- grid, macroscopic fields (`u`, `rho`), boundary mask, inlet/outlet pressures |
| `backend.hpp` | `lbm::ExecutionBackend`, and the `detail::direction/weight/opposite` accessors that switch between host tables and CUDA `__constant__` memory |
| `boundaries.hpp` | `lbm::Solid` -- boundary type constants, `compute_boundary_mask()`, and one function per boundary condition |
| `functions.hpp` | `lbm::functional` -- centerline profile extractors passed to `LBMSimulation::output()` |
| `core/point.hpp`, `core/vector.hpp`, `core/operators.hpp` | immutable `Point`, mutable `Vector`, their mixed-type arithmetic and `dot`/`cross` |
| `core/types.hpp` | `dim_t`, `Coordinate`, `DimPoint`, `ValuePoint`, `boundary_mask_t` |
| `core/grid.hpp` | `lbm::Grid` -- extents plus the two index maps (`scalar_index`, `field_index`) |
| `core/velocity-sets.hpp` | `lbm::D2Q9` (weights, directions, opposites) and the CUDA `__constant__` mirrors |
| `collision-operators/metadata.hpp` | `CollisionModel` and the per-model `CollisionParams`, where `nu`, `tau`, `s_plus`, `s_minus` are derived |
| `collision-operators/collision-strategy.hpp` | `lbm::CollisionStrategy` -- BGK and TRT kernels, selected with `if constexpr` |
| `solver/solver-base.hpp` | `SolverBase`, `SolverBase2D` -- iteration/frame bookkeeping and the `solve()` contract |
| `solver/openmp-solver.hpp` | `MPISolver2D` -- the OpenMP fused stream-collide loop |
| `solver/cuda-solver.cuh` | `CUDASolver2D` and its two `__global__` kernels |
| `data/data-listener.hpp`, `data/data-observable.hpp` | the observer pair `IDataListener` / `DataObservable` |
| `data/async-binary-writer.hpp` | `AsyncBinaryWriter` -- queue plus writer thread |
| `analysis/types.hpp` | `analysis::Function`, `NormType`, `NormErrorResult` |
| `analysis/exact-solution.hpp` | `CouetteSolution2D`, `PoiseuilleSolution2D` |
| `analysis/error.hpp` | `ErrorEvaluator` -- per-cell difference and global norm reduction |
| `analysis/benchmarks.hpp` | Ghia et al. (1982) table parsing and `compute_ghia_error()` |
| `collision-detection/` | CRTP `Shape` hierarchy, `CollisionArea`, Bresenham rasterization |
| `cuda/annotations.hpp`, `omp/annotations.hpp` | the `LBM_HD_FUNC` and `UNROLL_FULL` portability macros |
| `cuda/utils.cuh` | `LBM_CUDA_CHECK` and `ceil_div` |

## Building

```bash
cmake -S . -B build && cmake --build build -j
```

Every `.cpp` under `simulations/` becomes an executable named after the file;
with `-DLBM_ENABLE_CUDA=ON` every `.cu` becomes an executable prefixed with
`cuda_`. The options that matter here:

| Option | Default | Meaning |
|--------|---------|---------|
| `LBM_ENABLE_CUDA` | `OFF` | enable the CUDA language and build `lbm-sim-cuda` plus the `cuda_*` executables |
| `CMAKE_CUDA_ARCHITECTURES` | `75` | target GPU architecture |
| `LBM_LOG_LEVEL` | `INFO` | compile-time Quill level; statements below it are compiled out entirely |

## Building this documentation

```bash
cmake --build build --target lbm-docs-full
```

`lbm-sim-docs` builds this manual -- these pages plus the annotated headers --
into `docs/doxygen/lbm-sim/`; `lbm-docs-full` additionally builds the portal in
`docs/doxygen/`. Only `*.hpp` headers are fed to Doxygen unless
`LBM_ENABLE_CUDA=ON`, which adds the `*.cuh` ones.
