# Performance {#performance}

[TOC]

How a run is timed, how a scaling sweep is driven and plotted, and what the
sweeps measured on a 64-core node. The physics side of the same runs is in
@ref validation.

## Instrumenting a run

The instrumentation is compiled in only when the project is configured with
`-DLBM_ENABLE_PROFILING=ON`, which defines `LBM_PROFILING`. Without it
`PROFILE_SCOPE` expands to nothing and the `Profiler` never opens a file, so an
un-profiled build pays nothing.

`OpenMPSolver::solve()` carries three timers:

| Label | Scope | Reads as |
|---|---|---|
| `init_equilibrium` | the initial fill of the populations | one-off cost, outside the loop |
| `solve_total` | the whole iteration loop | the number a scaling plot is about |
| `stream_collide` | the fused kernel, accumulated over every iteration | `calls` is `niters`, so `avg` is the per-iteration cost |

At the end of `solve()` one CSV row per timer is appended through
`profiling::Profiler<ProfilingSchemaOpenMP>`:

```
label,size,collision_model,backend,n_threads,total,avg,calls
```

`total` and `avg` are milliseconds. The gap between `solve_total` and
`stream_collide` is everything the loop does besides the kernel -- essentially
the frame emission -- so the two nearly coincide when benchmark mode is on.

The file has to be opened by the main; nothing opens it implicitly:

```cpp
auto &profiler = profiling::Profiler<ProfilingSchemaOpenMP>::get();
profiler.open("out/prof/" + configs[0].name + job_id + ".csv");
```

Opening with `append = true` from several executables lets one CSV accumulate
a whole sweep: the header is written only into an empty file.

## Controlling the thread count

`profiling::BackendProperties<OPEN_MP>` is the process-wide handle on the
OpenMP knobs. It is filled once from the configuration and pushed into the
runtime around the solve:

```cpp
auto &prop = profiling::BackendProperties<ExecutionBackend::OPEN_MP>::get();
prop.setNumThreads(cfg.n_threads);   // [conf.backend].n_threads, 0 = leave it alone
prop.setBenchmarkMode(true);         // no frames: time the solver, not the writer
{
  const auto scope = prop.scopedApply();   // restores the previous ICVs on exit
  simulation.solve(solver);
}
```

Two properties of this design matter when reading a scaling table:

- `setDynamicThreads()` defaults to **off**. With dynamic adjustment on, the
  runtime may hand a parallel region fewer threads than requested, which makes
  a timing run unreproducible; the default exists so that a sweep measures what
  it says it measures.
- The `n_threads` column of the CSV is `props.getNumThreads()` -- the value
  that was **requested**, not one measured inside the parallel region with
  `omp_get_num_threads()`. A job that is granted fewer CPUs than the sweep asks
  for still writes the requested number into the file. When a scaling curve
  looks impossible, this is the first thing to rule out.

`simulations/openmp/profiling_lid_cavity_d2q9.cpp` is the main that does all of
the above; it is the one the sweeps below were run with.

## Running a sweep

The sweep itself is a list of `[[conf]]` blocks, one per point, differing only
in `[conf.backend].n_threads` (strong) or in that plus `[conf.lattice].size`
(weak):

- [`configs/profiling.toml`](../../configs/profiling.toml) -- 129x129,
  `Re = 100`, `niters = 100000`, threads 1 to 64;
- [`configs/weak_scaling.toml`](../../configs/weak_scaling.toml) -- same
  iteration count, grid grown alongside the thread count.

One repetition is one PBS job; `submit_simulations.sh` submits N of them and
each leaves its CSVs under `results/<JOBID>/prof/`:

```bash
./scripts/submit_simulations.sh build/simulations/profiling_lid_cavity_d2q9 configs/profiling.toml 5
```

```bash
python scripts/py/plot_strong.py results -o docs/report/assets/strong_scaling.png
```

```bash
python scripts/py/plot_weak.py results -o docs/report/assets/weak_scaling.png
```

Both plotters walk the directories recursively, average the repetitions, and
write two images: the figure and the summary table. Two rules they follow are
worth knowing before reading their output:

- **The baseline is the smallest thread count actually present**, per grid
  (strong) or per series (weak) -- never a hardcoded 1. A sweep that starts at
  2 threads is scored against its own first point.
- **`plot_weak.py` splits repetitions into series** by the work per thread at
  their baseline, so two sweeps with different targets landing in the same
  `results/` tree are averaged separately instead of being blended. The figure
  shows the series; the table image does not carry the series column.

## Strong scaling

Fixed 129x129 grid, `Re = 100`, 100 000 iterations, 5 repetitions.

![Strong scaling](../report/assets/strong_scaling.png)

| threads | mean [s] | std [s] | speedup | efficiency | MLUPS |
|--------:|---------:|--------:|--------:|-----------:|------:|
| 1  | 140.459 | 0.792 | 1.00  | 1.00 | 11.8 |
| 2  |  72.565 | 1.086 | 1.94  | 0.97 | 22.9 |
| 4  |  37.437 | 0.256 | 3.75  | 0.94 | 44.5 |
| 8  |  19.560 | 0.252 | 7.18  | 0.90 | 85.1 |
| 16 |  10.567 | 0.271 | 13.29 | 0.83 | 157.5 |
| 32 |   7.207 | 0.523 | 19.49 | 0.61 | 230.9 |
| 64 |  14.092 | 5.891 | 9.97  | 0.16 | 118.1 |

MLUPS is `nx * ny * niters / time`, in millions of lattice updates per second.
It is the metric to carry across grids: speedup compares a run only against
itself.

@note The wall times, speedups and efficiencies in this table and in the weak
one below are transcribed from the summary images the plotters write,
`docs/report/assets/strong_scaling_table.png` and `weak_scaling_table.png`.
Those images are the record of the sweep; the MLUPS columns are derived from
them. Re-running a sweep replaces the images, and both tables have to be
checked against them.

The node loop holds up to 16 threads (efficiency 0.83) and still pays at 32
(0.61). At 64 it *regresses*: the time doubles back to 14.1 s and the standard
deviation reaches 5.9 s, 40% of the mean. A slower mean with an unstable spread
is contention, not saturation -- with 16641 / 64 = 260 nodes per thread the
loop body is shorter than the barrier that closes it, and 64 logical CPUs over
fewer physical cores puts two threads on one set of execution units. Neither is
something the solver can fix: bind the threads (`OMP_PROC_BIND=close`,
`OMP_PLACES=cores`) and give 64 threads a grid worth their while.

@warning The 300x300 series in the figure is flat -- 78.0 s on 1 thread against
75.5 s on 64, a speedup of 1.03 -- and should not be quoted as a property of
the solver. A bandwidth-bound loop still improves several-fold before
flattening; a curve flat *from the first doubling* means the extra threads
never ran, and the 64-thread run has the same throughput as its own serial one.
Because the CSV records the requested thread count (see above), a job granted
fewer CPUs than it asked for is indistinguishable from one that scaled badly.
Re-run it with `select=...:ncpus=` matching the sweep.

## Weak scaling

Grid grown with the thread count, same 100 000 iterations per point.

![Weak scaling](../report/assets/weak_scaling.png)

| threads | grid | cells/thread | mean [s] | MLUPS | MLUPS per thread |
|--------:|------|-------------:|---------:|------:|-----------------:|
| 1  | 129x129   | 16641 | 141.774 | 11.7  | 11.7 |
| 2  | 183x183   | 16744 | 146.129 | 22.9  | 11.5 |
| 4  | 257x257   | 16512 | 146.945 | 45.0  | 11.2 |
| 8  | 365x365   | 16653 | 151.584 | 87.9  | 11.0 |
| 16 | 515x515   | 16577 | 151.768 | 174.8 | 10.9 |
| 32 | 1000x1000 | 31250 | 325.416 | 307.3 |  9.6 |
| 64 | 2000x2000 | 62500 | 1447.020| 276.4 |  4.3 |

Weak scaling means constant work per thread, so the ideal is a **flat** time
curve rather than a falling one. The first five points are exactly that: about
16 600 nodes per thread throughout, 141.8 s against 151.8 s from 1 to 16
threads, efficiency 0.934.

@note The last two rows are a second series, not a continuation of the first.
`weak_scaling.toml` gives 32 threads 31 250 nodes each and 64 threads 62 500 --
the work per thread doubles, then doubles again -- so their efficiency column
restarts from its own 32-thread baseline and reads 1.000 by construction.
`plot_weak.py` draws them as a separate series for this reason; the table image
does not show that split.

Per-thread throughput sidesteps the ambiguity, because it is comparable across
every row: 11.7 to 10.9 MLUPS/thread from 1 to 16 threads (a 7% loss), 9.6 at
32 (82% of the serial rate) and 4.3 at 64 (37%).

The last drop is the memory hierarchy, and the footprint says so. Two
population buffers of `ndir = 9` doubles are 144 bytes per node, plus roughly
28 bytes for `u`, `rho` and the norm buffer:

| grid | populations | working set |
|------|------------:|------------:|
| 515x515   |  38 MB |  46 MB |
| 1000x1000 | 144 MB | 172 MB |
| 2000x2000 | 576 MB | 688 MB |

Up to 515x515 the whole simulation fits in the aggregate last-level cache of
the node; at 2000x2000 every iteration streams 688 MB from DRAM twice, once to
read `ffrom` and once to write `fto`. A stream-collide step does a few dozen
flops against 144 bytes of traffic per node, so once the working set leaves
cache the memory controllers set the pace -- which is why the 1-to-16 sweep
loses 7% and the 64-thread point loses a factor of 2.5.

## Before the next sweep

- Keep the work per thread constant across the whole weak series: 16 641 nodes
  per thread means 730x730 at 32 threads and 1032x1032 at 64. The 2000x2000 run
  is worth keeping as its own series -- it is the only point that measures the
  solver against DRAM rather than against cache.
- Record the thread count *measured* inside the parallel region, not the one
  requested, so a mis-sized job is visible in the CSV.
- Raise the walltime in `scripts/cpu_job.pbs` before adding repetitions at
  2000x2000: `solve_total` alone is 1447 s of the 1800 s budget, leaving no
  margin for initialization and for copying results back.
- Bind threads (`OMP_PROC_BIND`, `OMP_PLACES`) for the 32- and 64-thread points
  and check whether the 64-thread regression survives.
