#!/usr/bin/env python3
"""Plot the weak scaling curve from the results left behind by the PBS jobs.

Reads every results/<JOBID>/prof/weak_scaling.csv produced by run_weak.pbs --
one file per repetition, one row per configuration -- keeps the `solve_total`
row, which is the wall time of the whole iteration loop, and averages across
repetitions.

    python3 scripts/cluster/plot_weak.py
    python3 scripts/cluster/plot_weak.py --results results -o weak.png
    python3 scripts/cluster/plot_weak.py --show

Weak scaling means the lattice grows with the thread count so that the work per
thread stays constant. The ideal is therefore a FLAT time curve, not a falling
one: any slope is parallel overhead plus the cost of the larger footprint.
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


def cells(size):
    """Cell count from a "129x129" tag, 0 if it cannot be parsed."""
    try:
        n = 1
        for part in size.split("x"):
            n *= int(part)
        return n
    except ValueError:
        return 0


def summarize(by_threads):
    """Per thread count: mean, sample std and weak efficiency.

    Weak efficiency is T(base)/T(N): with constant work per thread a perfect
    run takes the same time everywhere, so this sits at 1. The baseline is the
    smallest thread count actually present rather than a hardcoded 1.

    cells/thread is reported too, because the whole premise of the plot is that
    it stays constant: if the config drifts, the curve stops meaning anything
    and this column is where it shows.
    """
    threads = sorted(by_threads)
    base = threads[0]
    base_mean = sum(by_threads[base][0]) / len(by_threads[base][0])

    stats = []
    for n in threads:
        vals = by_threads[n][0]
        size = by_threads[n][1]
        mean = sum(vals) / len(vals)
        if len(vals) > 1:
            var = sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)
            std = var ** 0.5
        else:
            # One repetition measures no spread. 0 draws no error bar, which is
            # honest; the run count is printed in the table below.
            std = 0.0
        stats.append({
            "n_threads": n,
            "runs": len(vals),
            "size": size,
            "cells_per_thread": cells(size) / n,
            "mean_s": mean / 1000.0,
            "std_s": std / 1000.0,
            "efficiency": base_mean / mean,
        })
    return stats, base


def print_table(stats):
    print("%-9s %5s %-11s %14s %10s %9s %11s"
          % ("n_threads", "runs", "size", "celle/thread", "mean [s]",
             "std [s]", "efficienza"))
    for s in stats:
        print("%-9d %5d %-11s %14.0f %10.3f %9.3f %11.3f"
              % (s["n_threads"], s["runs"], s["size"], s["cells_per_thread"],
                 s["mean_s"], s["std_s"], s["efficiency"]))


def plot(stats, out_path, show):
    import matplotlib
    if not show:
        # No DISPLAY on a compute or login node: pick the file-only backend
        # before pyplot is imported, or the import itself fails there.
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    threads = [s["n_threads"] for s in stats]
    mean = [s["mean_s"] for s in stats]
    std = [s["std_s"] for s in stats]
    efficiency = [s["efficiency"] for s in stats]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axes[0]
    ax.errorbar(threads, mean, yerr=std, fmt="o-", capsize=4, label="misurato")
    # Constant work per thread means a perfect run never gets slower: the ideal
    # is the baseline time held flat, not a descending line as in strong scaling.
    ax.axhline(mean[0], color="k", ls="--", alpha=.6, label="ideale (piatto)")
    ax.set_xlabel("thread")
    ax.set_ylabel("tempo del solve [s]")
    ax.set_title("Tempo")
    ax.set_xscale("log", base=2)
    # Headroom for the grid labels annotated above each point, which would
    # otherwise be clipped by the axes box.
    ax.set_ylim(0, max(mean) * 1.25)
    # The grid is what changes between points here, so label each with its own.
    for s, t, m in zip(stats, threads, mean):
        ax.annotate(s["size"], (t, m), textcoords="offset points",
                    xytext=(0, 9), ha="center", fontsize=8)

    ax = axes[1]
    ax.plot(threads, efficiency, "o-")
    ax.axhline(1.0, color="k", ls="--", alpha=.6)
    ax.set_xlabel("thread")
    ax.set_ylabel("efficienza  T(%d)/T(N)" % threads[0])
    ax.set_title("Efficienza")
    ax.set_xscale("log", base=2)
    # Weak efficiency can exceed 1 when the larger lattice happens to vectorise
    # better; a hardcoded ceiling would clip those points off without saying so.
    ax.set_ylim(0, max(1.1, max(efficiency) * 1.05))

    for ax in axes:
        # Powers of two on a log axis default to 2^k labels; the plain thread
        # counts are what the reader is looking for.
        ax.set_xticks(threads)
        ax.set_xticklabels([str(t) for t in threads])
        ax.grid(True, which="both", alpha=.2)
    axes[0].legend()

    cpt = stats[0]["cells_per_thread"]
    runs = max(s["runs"] for s in stats)
    fig.suptitle("Weak scaling -- ~%.0f celle/thread, %d ripetizioni"
                 % (cpt, runs))
    fig.tight_layout()

    fig.savefig(out_path, dpi=150)
    print("scritto %s" % out_path)
    if show:
        plt.show()


def main():
    ap = argparse.ArgumentParser(
        description="Grafico di weak scaling dai risultati dei job PBS.")
    ap.add_argument("--results", default="results",
                    help="cartella con una sottodirectory per job (default: results)")
    ap.add_argument("--csv-name", default="weak_scaling.csv",
                    help="nome del CSV del profiler dentro <JOBID>/prof "
                         "(default: weak_scaling.csv)")
    ap.add_argument("-o", "--output", default="results/weak_scaling.png",
                    help="immagine da scrivere (default: results/weak_scaling.png)")
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

    stats, _ = summarize(by_threads)
    print("%d ripetizioni lette da %s" % (n_files, args.results))
    print_table(stats)
    plot(stats, args.output, args.show)


if __name__ == "__main__":
    main()
