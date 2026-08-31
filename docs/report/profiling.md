# Profiling
## Scaling results

Strong scaling -- 129x129 grid, `Re = 100`, 10 000 steps (red: speedup, blue:
wall time in seconds):

![strong scaling](assets/strong_scaling.png)

Weak scaling -- problem size grown with the thread count (blue: wall time in
seconds), over 50^2 and 100^2 at `Re = 100`, 150^2 and 200^2 at `Re = 500`,
250^2 and 300^2 at `Re = 1000`:

![weak scaling](assets/weak_scaling.png)

> These plots were produced with `scripts/strong_scaling.py` and
> `scripts/weak_scaling.py` for the 2D OpenMP solver of the original hands-on.
> Both scripts still emit the old flat configuration format and need updating
> before they can drive the current binaries.

