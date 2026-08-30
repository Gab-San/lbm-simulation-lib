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
size     = [0.1, 0.0]      # reference velocity; size[0] must be > 0

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
| `conf.physics.size` | array of `dim` floats | yes | `SimulationConfig::u0` |
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

### `physics.size` is more than an initial condition

`u0` is copied into `lbm::CollisionParams`, which derives the kinematic
viscosity from it:

\f[
  \nu = \frac{u_0 \, N_y}{\mathrm{Re}}, \qquad \tau = \frac{1}{2} + 3\nu
\f]

so `size[0]` (that is, `u0.dx`) has to be strictly positive **even in problems
where no wall moves**, such as Poiseuille or pipe flow. There it is the
characteristic velocity, and the value the imposed pressure drop is expected to
reproduce. The parser enforces `size[0] > 0` and says why.

`CollisionParams` additionally throws if the resulting `tau <= 0.5`, and warns
outside the roughly `0.55 - 1.2` band where the scheme is comfortably stable.

### Pressures and geometry are not configuration

The inlet/outlet pressure drop of the Poiseuille and pipe cases is **derived**
by the main from the analytical solution, not read from the file — it is fixed
by `Re`, `u0` and the geometry, so exposing it as a free parameter would only
allow inconsistent inputs. The same goes for the pipe radius, which the main
builds as a `lbm::CollisionDetection::CylindricalShell`. The comments at the
top of `configs/pipe_3d_128x65x65_re100_bgk.toml` and
`configs/poiseuille_129x129_re100_bgk.toml` spell out the formulas.

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

The messages name the offending key, and the file-level ones name the path:

- `configuration file does not exist: <path>`
- `malformed TOML in <path>: <toml++ description>`
- `field 'conf' is missing or is not an array in <path>`
- `config '<name>': [lattice].size is mandatory`
- `config '<name>': [physics].init_vel_x must be > 0 (it is the characteristic velocity used for nu = u*Ny/Re)`

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

- `collision` and `backend`, which appear as per-entry keys in
  `configs/couette.toml` and in the older examples, are **not read** by
  `parse_table()`. The collision operator and the backend are compile-time
  parameters of the binary; the keys are currently documentation only.
- `lbm::config::ensure_compatible()` is written against
  `SimulationConfig::collision` and `SimulationConfig::backend`, which do not
  exist as members. No main calls it, so it never gets instantiated; it needs
  either those fields or removal before it can be used.
- Several files under `configs/` still use an older `[grid] nx/ny/nz` +
  `[physics] init_vel_x` layout that predates the `[[conf]]` array. They parse
  as valid TOML but produce a `ConfigError`. Only `couette.toml` and
  `profiling.toml` follow the current `[[conf]]` schema.
