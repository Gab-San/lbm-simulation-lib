# Output formats {#output_formats}

[TOC]

A run produces up to four artefacts: the **frame stream**, the **profile
file**, optionally a **VTK series**, and — in a profiling build — a
**profiling CSV**. This page gives their exact layout, so a reader can be
written without inspecting the sources.

All binary payloads use the **native byte order and the native `double` /
`float` representation**. They are meant to be read back on a matching
platform; nothing in the format records endianness.

## The frame stream

Emitted through the listener chain (@ref architecture), so its shape is
independent of which listener consumes it. The producers are two distinct
observables:

| Producer | When | Payload |
|---|---|---|
| `lbm::LBMSimulation::solve()` | once, before the solver runs | the grid extents |
| `lbm::OpenMPSolver` / `lbm::CUDASolver` | every `niters / nframes` steps | one frame of velocity norms |

### Header chunk

`dim` values of `int32_t`, native order:

```
2D:  int32 nx, int32 ny                 (8 bytes)
3D:  int32 nx, int32 ny, int32 nz      (12 bytes)
```

`LBMSimulation::write_header()` casts from the grid's `size_t` extents
deliberately, so the wire format does not follow the platform's `size_t`
width.

### Frame chunk

`nx * ny [* nz]` values of `float32`, one per node, in the order used by
`lbm::Grid::scalar_index()` (x fastest). Each value is the velocity
magnitude at that node. `lbm::OpenMPSolver::write_norms()` `memcpy`s the
scratch array straight into the buffer; no per-frame header, no separator.

A file written by `lbm::AsyncBinaryWriter` is therefore exactly

```
[header][frame 0][frame 1] ... [frame nframes-1]
```

and the frame count is recovered by division:

```python
n_frames = len(data) // (nx * ny)      # data read as float32 after the header
```

which is what `scripts/py/visualize_frames.py` does.

`AsyncBinaryWriter` writes on a dedicated thread and flushes every 32 chunks;
the destructor drains the queue and closes the file, so the writer must
outlive the run — hold it in a `std::shared_ptr` and detach it afterwards.

### VTK series

`lbm::VtkWriter` consumes the same stream and turns it into one `.vti`
(ImageData XML) per frame, named `<basename>_00000.vti`, `_00001.vti`, … plus
a `.pvd` series file. The scalar array is `velocity_magnitude`, `Float32`,
stored in `<AppendedData encoding="raw">` with a `UInt32` byte count in front,
so ParaView reads the payload without any numeric conversion on either side.

The `.pvd` is rewritten every 10 frames (and once more when the worker has
drained the queue), which is what lets ParaView open the series **while the
run is still going** without paying one rewrite per frame.

A frame whose size does not match `nx * ny * nz * sizeof(float)` is reported
on `stderr` and skipped rather than corrupting the series.

## The profile file

Written by `lbm::LBMSimulation::output()`, which creates the parent directory
chain if needed. One ASCII header line, then the payload:

```
%%profile <MODEL> <N> <U_REF>\n
<N × float64, native order>
```

- `MODEL` — `collision_model_to_string(cm_t)`: `BGK`, `TRT` or `MRT`;
- `N` — number of samples in the profile;
- `U_REF` — `params.init_vel.dx`, the reference velocity Ghia's tables are
  normalised by.

The samples themselves come from the extractor passed as the second argument.
The ready-made ones are in `lbm::functional`:

| Extractor | Returns |
|---|---|
| `extract_dx_profile_along_y_center()` | `u_x(y)` on the column `x = nx/2` (2D) |
| `extract_dy_profile_along_x_center()` | `u_y(x)` on the row `y = ny/2` (2D) |
| `extract_dx_profile_along_z_center()` | `u_x(z)` on the column `x = nx/2, y = ny/2` (3D) |

Any callable with signature `std::vector<double>(const Lattice<dim> &)` works.

@warning If the destination cannot be opened the failure is logged and
`output()` returns silently; the caller cannot detect it. This is tracked by a
`FIXME` in `lbm-simulation.hpp`.

## The profiling CSV

Only in a build configured with `-DLBM_ENABLE_PROFILING=ON`, which defines
`LBM_PROFILING`. `PROFILE_SCOPE(name)` accumulates wall time per label into a
process-wide registry; the solver dumps it at the end of `solve()` through
`lbm::formatting::CsvWriter`, using `lbm::ProfilingSchemaOpenMP`:

```
label,size,collision_model,backend,n_threads,total,avg,calls
```

with `total` in milliseconds formatted to two decimals, `avg = total / calls`.
The labels the OpenMP solver registers are `init_equilibrium`,
`stream_collide` and `solve_total` (the wall time of the whole loop).

`CsvWriter` is synchronous and formats on the calling thread, so no row can be
lost by skipping a shutdown. Opened with `append = true` it writes the header
only when the file was empty, which is how a sweep running one executable per
configuration accumulates into a single file.

`lbm::analysis::ErrorAnalysisSchema` uses the same mechanism for error tables:

```
profile,size,collision_model,niters,type,rel,abs
```

## Reading the outputs

```bash
python scripts/py/visualize_frames.py out/frames_lid_cavity.bin -o cavity.gif
```

```bash
python scripts/py/visualize_profile.py out/profile_lid_cavity.dat \
       benchmarks/ghia/data_y_100.txt --title "Re = 100" -o profile.png
```

```bash
python scripts/py/validate_outputs.py
```

`validate_outputs.py` is the structural check: headers, payload lengths,
finite values, and output paths declared twice by different sources.

Ghia's tabulated values are already normalised by the lid velocity, so
`visualize_profile.py` plots them as they are; rescaling them by their own
maximum would distort the profile shape and hide a magnitude mismatch.
