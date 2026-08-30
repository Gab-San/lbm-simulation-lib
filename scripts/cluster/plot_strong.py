#!/usr/bin/env python3
"""Plot the strong scaling curve from the results left behind by the PBS jobs.

Reads every results/<JOBID>/prof/profiling.csv produced by run_properties.pbs --
one file per repetition, one row per thread count -- keeps the `solve_total`
row, which is the wall time of the whole iteration loop, and averages across
repetitions.

    python3 scripts/cluster/plot_strong.py
    python3 scripts/cluster/plot_strong.py --results results -o strong.png
    python3 scripts/cluster/plot_strong.py --show

Strong scaling means the problem stays the same size while threads are added,
so the ideal is time halving each time the thread count doubles.
"""

import argparse
import csv
import os
import sys

# Columns written by lbm::ProfilingSchemaOpenMP.
FIELDS = ("label", "size", "collision_model", "backend",
          "n_threads", "total", "avg", "calls")

# The solver registers three timers; this is the one a scaling plot is about.
LABEL = "solve_total"


def read_csv(path):
    """(n_threads, total_ms, size) for the solve loop, skipping bad rows.

    A job killed by the walltime can leave a half-written line: it is dropped
    with a warning rather than failing the whole plot, which would throw away
    every repetition that did complete.
    """
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or tuple(reader.fieldnames) != FIELDS:
            print("attenzione: intestazione inattesa in %s, file saltato" % path,
                  file=sys.stderr)
            return rows
        for row in reader:
            if row.get("label") != LABEL:
                continue
            try:
                rows.append((int(row["n_threads"]), float(row["total"]),
                             row["size"]))
            except (TypeError, ValueError):
                print("attenzione: riga scartata in %s: %r" % (path, row),
                      file=sys.stderr)
    return rows


def load(results_dir, csv_name):
    """{n_threads: (times_ms, size)} gathered over every repetition."""
    by_threads = {}
    n_files = 0
    for entry in sorted(os.listdir(results_dir)):
        path = os.path.join(results_dir, entry, "prof", csv_name)
        if not os.path.isfile(path):
            continue
        n_files += 1
        for n, total_ms, size in read_csv(path):
            times, _ = by_threads.setdefault(n, ([], size))
            times.append(total_ms)
    return by_threads, n_files


def summarize(by_threads):
    """Per thread count: mean, sample std, speedup and efficiency.

    Speedup is measured against the smallest thread count actually present, not
    against a hardcoded 1: a sweep starting at 2 threads would otherwise produce
    silently wrong numbers.
    """
    threads = sorted(by_threads)
    base = threads[0]
    base_mean = sum(by_threads[base][0]) / len(by_threads[base][0])

    stats = []
    for n in threads:
        vals = by_threads[n][0]
        mean = sum(vals) / len(vals)
        if len(vals) > 1:
            var = sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)
            std = var ** 0.5
        else:
            # One repetition measures no spread. 0 draws no error bar, which is
            # honest; the run count is printed in the table below.
            std = 0.0
        speedup = base_mean / mean
        stats.append({
            "n_threads": n,
            "runs": len(vals),
            "size": by_threads[n][1],
            "mean_s": mean / 1000.0,
            "std_s": std / 1000.0,
            "speedup": speedup,
            "efficiency": speedup / (n / base),
        })
    return stats, base


def print_table(stats):
    print("%-9s %5s %-11s %10s %9s %9s %11s"
          % ("n_threads", "runs", "size", "mean [s]", "std [s]",
             "speedup", "efficiency"))
    for s in stats:
        print("%-9d %5d %-11s %10.3f %9.3f %9.3f %11.3f"
              % (s["n_threads"], s["runs"], s["size"], s["mean_s"],
                 s["std_s"], s["speedup"], s["efficiency"]))


def plot(stats, base, out_path, show):
    import matplotlib
    if not show:
        # No DISPLAY on a compute or login node: pick the file-only backend
        # before pyplot is imported, or the import itself fails there.
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    threads = [s["n_threads"] for s in stats]
    mean = [s["mean_s"] for s in stats]
    std = [s["std_s"] for s in stats]
    speedup = [s["speedup"] for s in stats]
    efficiency = [s["efficiency"] for s in stats]
    ideal = [n / base for n in threads]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    ax = axes[0]
    ax.errorbar(threads, mean, yerr=std, fmt="o-", capsize=4, label="misurato")
    # Perfect scaling halves the time at each doubling: a straight line of slope
    # -1 through the first point, which log-log axes make easy to read off.
    ax.plot(threads, [mean[0] * base / n for n in threads], "k--",
            alpha=.6, label="ideale")
    ax.set_xlabel("thread")
    ax.set_ylabel("tempo del solve [s]")
    ax.set_title("Tempo")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")

    ax = axes[1]
    ax.errorbar(threads, speedup, fmt="o-", capsize=4, label="misurato")
    ax.plot(threads, ideal, "k--", alpha=.6, label="ideale")
    ax.set_xlabel("thread")
    ax.set_ylabel("speedup")
    ax.set_title("Speedup")
    ax.set_xscale("log", base=2)

    ax = axes[2]
    ax.plot(threads, efficiency, "o-")
    ax.axhline(1.0, color="k", ls="--", alpha=.6)
    ax.set_xlabel("thread")
    ax.set_ylabel("efficienza")
    ax.set_title("Efficienza")
    ax.set_xscale("log", base=2)
    # Superlinear speedup is real (the smaller per-thread working set starts
    # fitting in cache), so a hardcoded 1.1 ceiling would clip those points off
    # the plot without saying so.
    ax.set_ylim(0, max(1.1, max(efficiency) * 1.05))

    for ax in axes:
        # Powers of two on a log axis default to 2^k labels; the plain thread
        # counts are what the reader is looking for.
        ax.set_xticks(threads)
        ax.set_xticklabels([str(t) for t in threads])
        ax.grid(True, which="both", alpha=.2)
        if ax is not axes[2]:
            ax.legend()

    size = stats[0]["size"]
    runs = max(s["runs"] for s in stats)
    fig.suptitle("Strong scaling -- griglia %s, %d ripetizioni" % (size, runs))
    fig.tight_layout()

    fig.savefig(out_path, dpi=150)
    print("scritto %s" % out_path)
    if show:
        plt.show()


def main():
    ap = argparse.ArgumentParser(
        description="Grafico di strong scaling dai risultati dei job PBS.")
    ap.add_argument("--results", default="results",
                    help="cartella con una sottodirectory per job (default: results)")
    ap.add_argument("--csv-name", default="profiling.csv",
                    help="nome del CSV del profiler dentro <JOBID>/prof "
                         "(default: profiling.csv)")
    ap.add_argument("-o", "--output", default="results/strong_scaling.png",
                    help="immagine da scrivere (default: results/strong_scaling.png)")
    ap.add_argument("--show", action="store_true",
                    help="apre la finestra invece di limitarsi a salvare")
    args = ap.parse_args()

    if not os.path.isdir(args.results):
        sys.exit("cartella dei risultati inesistente: %s" % args.results)

    by_threads, n_files = load(args.results, args.csv_name)
    if not by_threads:
        sys.exit("nessun dato in %s/*/prof/%s "
                 "(il binario e' stato compilato con LBM_ENABLE_PROFILING=ON?)"
                 % (args.results, args.csv_name))

    stats, base = summarize(by_threads)
    print("%d ripetizioni lette da %s" % (n_files, args.results))
    print_table(stats)
    plot(stats, base, args.output, args.show)


if __name__ == "__main__":
    main()
