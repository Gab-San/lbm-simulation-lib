# Configuration files {#configuration}

[TOC]

A configuration file is TOML. It holds a top-level array of tables named
`conf`, one entry per run, and is read by
`lbm::config::parse_config<dim>()`, which returns a
`std::vector<lbm::config::SimulationConfig<dim>>`.

Everything in this page is what `parse_table()` in `config/config-parser.hpp`
actually reads. Keys not listed here are accepted by the TOML parser and then
ignored.

## Schema

```toml
[[conf]]

[conf.lattice]
size = [129, 129]          # [nx, ny] in 2D, [nx, ny, nz] in 3D

[conf.physics]
reynolds = 100.0           # > 0
init_vel = [0.1, 0.0]      # reference velocity, one component per dimension

[conf.solver]
niters  = 100000           # > 0
nframes = 200              # must not exceed niters

[conf.output]
frames  = "out/frames_lid_cavity.bin"
profile = "out/profile_lid_cavity.dat"

[conf.backend]
n_threads = 8              # optional, defaults to 0 ("let OpenMP decide")
```

Repeat the `[[conf]]` block to describe several runs in one file; the main
loops over the vector it gets back. `configs/profiling.toml` is an example of
a thread-count sweep written this way.

## Field reference

| Key | Type | Required | Maps to |
|---|---|---|---|
| `conf.lattice.size` | array of `dim` integers | yes | `SimulationConfig::grid_size` |
| `conf.physics.reynolds` | float `> 0` | yes | `SimulationConfig::reynolds` |
| `conf.physics.init_vel` | array of `dim` floats | yes | `SimulationConfig::u0` |
| `conf.solver.niters` | integer `> 0` | yes | `SimulationConfig::niters` |
| `conf.solver.nframes` | integer `<= niters` | no (default `0`) | `SimulationConfig::nframes` |
| `conf.output.frames` | string | yes | `SimulationConfig::frames_out` |
| `conf.output.profile` | string | yes | `SimulationConfig::profile_out` |
| `conf.backend.n_threads` | integer | no (default `0`) | `SimulationConfig::n_threads` |

`SimulationConfig::name` is not a field either: the parser fills it with the
stem of the file name, so error messages and output basenames identify the run
without a redundant key.

### The dimension is not a field

`parse_config<dim>()` is instantiated with the dimension the binary was
compiled for. An array whose length does not match `dim` — a 3-element
`lattice.size` handed to a 2D binary — is rejected with a
`lbm::config::ConfigError` naming both. There is no key that changes `dim` at
run time.

### `physics.init_vel` is more than an initial condition

`u0` is copied into `lbm::CollisionParams`, which derives the kinematic
viscosity from it:

\f[
  \nu = \frac{u_0 \, N_y}{\mathrm{Re}}, \qquad \tau = \frac{1}{2} + 3\nu
\f]

so `init_vel[0]` (that is, `u0.dx`) has to be strictly positive **even in
problems where no wall moves**, such as Poiseuille or pipe flow. There it is
the characteristic velocity, and the value the imposed pressure drop is
expected to reproduce.

@warning The parser does **not** check this: it validates the array length and
that every element is a number, nothing more. A zero or negative `init_vel[0]`
gets through and produces `nu <= 0` downstream, where `CollisionParams` throws
on `tau <= 0.5`. The failure is reported, but by the collision parameters and
not by the file that caused it.

`CollisionParams` additionally throws if the resulting `tau <= 0.5`, and warns
outside the roughly `0.55 - 1.2` band where the scheme is comfortably stable.

### Pressures and geometry are not configuration

The inlet/outlet pressure drop of the Poiseuille and pipe cases is **derived**
by the main from the analytical solution, not read from the file — it is fixed
by `Re`, `u0` and the geometry, so exposing it as a free parameter would only
allow inconsistent inputs. The same goes for the pipe radius, which the main
derives from the grid and builds as a
`lbm::CollisionDetection::CylindricalShell`.

The formulas live in the mains, next to the code that applies them:
`simulations/openmp/poiseuille_d2q9_bgk.cpp` inverts the channel solution for
the pressure drop, and `simulations/cuda/pipe_poiseuille_d3q19_bgk.cu` inverts
Hagen-Poiseuille for the density jump, deriving the effective radius as
`r_eff = radius + 0.5` and printing it, together with the transient time
scale, in its start-up log.

## Errors

Every failure path raises `lbm::config::ConfigError`, a distinct type rather
than a plain `std::runtime_error` precisely so that a main can tell "the user
passed a bad configuration" (print the message, exit 1) from an error raised
during the solve:

```cpp
std::vector<config::SimulationConfig<DIM>> configs;
try {
  configs = config::parse_config<DIM>(argv[1]);
} catch (const config::ConfigError &err) {
  std::cerr << "Configuration error: " << err.what() << "\n";
  return 1;
}
```

`lbm::config::print_usage()` prints the one-line usage banner the mains show
when they are given no argument.

The messages name the offending key, and the file-level ones name the path.
The complete set, as written by `parse_config()` and `parse_table()`:

| Message | Raised when |
|---|---|
| `configuration file does not exist: <path>` | the path does not exist |
| `malformed TOML in <path>: <toml++ description>` | toml++ fails to parse the file |
| `field 'conf' is missing or is not an array in <path>` | no `[[conf]]` array at the top level |
| `invalid element in 'conf' in <path>` | an entry of the array is not a table |
| `config '<name>': [lattice].size is mandatory` | the key is missing or not an array |
| `config '<name>': [lattice].size has 3 entries, but this binary is 2D` | the array length does not match `dim` |
| `config '<name>': [lattice].size[i] is missing or is not an integer` | a non-integer entry |
| `config '<name>': [physics].reynolds is mandatory` / `must be > 0` | missing, or `<= 0` |
| `config '<name>': [physics].init_vel is mandatory` | the key is missing or not an array |
| `config '<name>': [physics].init_vel has N entries, but this binary is <dim>D` | the array length does not match `dim` |
| `config '<name>': [solver].niters is mandatory and must be > 0` | missing or zero |
| `config '<name>': [solver].nframes (N) cannot exceed [solver].niters (M)` | more frames than steps |
| `config '<name>': [output].frames is mandatory` / `[output].profile is mandatory` | either output path is missing |

`<name>` is the stem of the configuration file, so a message identifies the run
without a redundant key in the file.

## Where the files live

`configs/` is passed to every executable as an absolute path at configure
time, through the `LBM_CONFIGS_DIR` compile definition set by
`simulations/CMakeLists.txt`, and the Ghia tables likewise through
`LBM_BENCHMARKS_DIR`. The files are *not* copied into the build tree, so
editing a `.toml` takes effect on the next run without re-running CMake, and a
binary finds its inputs from any working directory.

Output paths, by contrast, are taken from the file as written and are relative
to the working directory.

## Known gaps

These are recorded here so the page matches the code rather than the
intention:

- `collision = "bgk"` appears as a per-entry key in
  `configs/lid_cavity_cuda.toml`, and is **not read** by `parse_table()`. The
  collision operator and the backend are compile-time parameters of the
  binary; the key is documentation only, and nothing checks that it agrees
  with the executable it is handed to.
- Nothing validates `init_vel[0] > 0`, even though the viscosity is derived
  from it (see the warning above).
- `[conf.backend].n_threads` is parsed into `SimulationConfig::n_threads`, but
  only `simulations/openmp/profiling_lid_cavity_d2q9.cpp` applies it, through
  `BackendProperties<OPEN_MP>::setNumThreads()`. Every other main ignores the
  field and runs with whatever `OMP_NUM_THREADS` gives it, so a sweep driven by
  this key only works with the profiling executable.
- All nine files under `configs/` follow the current `[[conf]]` schema; the
  older flat `[grid] nx/ny/nz` layout is gone from the tree. There is no
  configuration for the 2D Poiseuille channel: `poiseuille_d2q9_{bgk,trt}`
  have to be pointed at a hand-written file, or at `configs/example.toml`.
