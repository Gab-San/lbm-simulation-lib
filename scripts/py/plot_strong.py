#!/usr/bin/env python3
"""Plot the strong scaling curve from the CSVs written by lbm's profiler.

Accepts any mix of files, globs and directories:

    python3 scripts/cluster/plot_strong.py                    # scans ./results
    python3 scripts/cluster/plot_strong.py results/1234/prof/profiling.csv
    python3 scripts/cluster/plot_strong.py run_a.csv run_b.csv -o strong.png
    python3 scripts/cluster/plot_strong.py 'results/*/prof/*.csv' --show
    python3 scripts/cluster/plot_strong.py results other_results --show

Directories are walked recursively, so the results/<JOBID>/prof/ layout the PBS
jobs leave behind is found without depending on it. Each file is one
repetition; rows are kept for the `solve_total` timer, which is the wall time
of the whole iteration loop, and averaged across repetitions.

Strong scaling means the problem stays the same size while threads are added,
so the ideal is time halving each time the thread count doubles.
"""

import argparse
import glob
import os
import sys

import pandas as pd

# Columns written by lbm::ProfilingSchemaOpenMP.
FIELDS = [
    "label",
    "size",
    "collision_model",
    "backend",
    "n_threads",
    "total",
    "avg",
    "calls",
]

# The solver registers three timers; this is the one a scaling plot is about.
LABEL = "solve_total"

def expand(paths, csv_name):
    """Every CSV named by the arguments: files, globs, directories, in order.

    Duplicates are dropped (a directory and one of the files inside it can both
    be passed) while the order the user gave is kept, so the reading warnings
    come out in a predictable sequence.
    """
    found, seen = [], set()

    def add(path):
        key = os.path.realpath(path)
        if key not in seen:
            seen.add(key)
            found.append(path)

    for arg in paths:
        if os.path.isdir(arg):
            for path in sorted(glob.glob(os.path.join(arg, "**", csv_name),
                                         recursive=True)):
                add(path)
        elif os.path.isfile(arg):
            add(arg)
        else:
            matches = sorted(glob.glob(arg, recursive=True))
            if not matches:
                print("attenzione: nessun file per %s" % arg, file=sys.stderr)
            for path in matches:
                if os.path.isdir(path):
                    for sub in sorted(glob.glob(os.path.join(path, "**", csv_name),
                                                recursive=True)):
                        add(sub)
                else:
                    add(path)
    return found


def read_one(path):
    """One repetition as a DataFrame, or None if the file is unusable.

    A job killed by the walltime can leave a half-written line: it is dropped
    with a warning rather than failing the whole plot, which would throw away
    every repetition that did complete.
    """
    try:
        df = pd.read_csv(path, index_col=False, on_bad_lines="skip")
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as err:
        print("attenzione: %s non interpretabile (%s)" % (path, err), file=sys.stderr)
        return None

    if list(df.columns) != FIELDS:
        print("attenzione: intestazione inattesa in %s, file saltato" % path,
              file=sys.stderr)
        return None

    df = df[df["label"] == LABEL].copy()
    df["n_threads"] = pd.to_numeric(df["n_threads"], errors="coerce")
    df["total"] = pd.to_numeric(df["total"], errors="coerce")

    # A line cut in half by the walltime keeps its leading fields, so
    # n_threads and total can both look fine while the tail is missing: an
    # incomplete row is not a measurement and every column has to be there.
    bad = df[FIELDS].isna().any(axis=1)
    for _, row in df[bad].iterrows():
        print("attenzione: riga scartata in %s: %s" % (path, row.to_dict()),
              file=sys.stderr)

    df = df[~bad]
    if df.empty:
        return None

    df["n_threads"] = df["n_threads"].astype(int)
    df["source"] = path
    return df


def load(paths):
    """One long frame with every row of every repetition, plus the file count."""
    frames = [df for df in (read_one(p) for p in paths) if df is not None]
    if not frames:
        return pd.DataFrame(columns=FIELDS + ["source"]), 0
    return pd.concat(frames, ignore_index=True), len(frames)


def summarize(df):
    """Per thread count: mean, sample std, speedup and efficiency.

    Speedup is measured against the smallest thread count actually present, not
    against a hardcoded 1: a sweep starting at 2 threads would otherwise produce
    silently wrong numbers.
    """
    dup = df.duplicated(["source", "n_threads"], keep=False)
    if dup.any():
        # Two solve_total rows for the same thread count in one file means the
        # binary was run twice into the same CSV; averaging them is fine, but a
        # silent average would hide a botched job script.
        print("attenzione: righe solve_total ripetute per lo stesso n_threads in %s"
              % ", ".join(sorted(df[dup]["source"].unique())), file=sys.stderr)

    stats = df.groupby("n_threads").agg(
        runs=("total", "size"),
        size=("size", "first"),
        mean_s=("total", "mean"),
        std_s=("total", "std"),
    )
    stats[["mean_s", "std_s"]] /= 1000.0
    # std() is NaN on a single repetition, which measures no spread. 0 draws no
    # error bar, which is honest; the run count is printed in the table below.
    stats["std_s"] = stats["std_s"].fillna(0.0)

    base = int(stats.index[0])
    stats["speedup"] = stats.loc[base, "mean_s"] / stats["mean_s"]
    stats["efficiency"] = stats["speedup"] / (stats.index / base)

    sizes = df["size"].unique()
    if len(sizes) > 1:
        print("attenzione: griglie diverse nello stesso grafico: %s"
              % ", ".join(map(str, sizes)), file=sys.stderr)
    return stats, base


def print_table(stats):
    out = stats.reset_index()[
        ["n_threads", "runs", "size", "mean_s", "std_s", "speedup", "efficiency"]
    ]
    print(out.to_string(
        index=False,
        formatters={
            "mean_s": "{:.3f}".format,
            "std_s": "{:.3f}".format,
            "speedup": "{:.3f}".format,
            "efficiency": "{:.3f}".format,
        },
    ))


def plot(stats, base, out_path, show):
    import matplotlib

    if not show:
        # No DISPLAY on a compute or login node: pick the file-only backend
        # before pyplot is imported, or the import itself fails there.
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    threads = stats.index.to_list()
    mean = stats["mean_s"].to_list()
    ideal = [n / base for n in threads]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    ax = axes[0]
    ax.errorbar(threads, mean, yerr=stats["std_s"], fmt="o-", capsize=4,
                label="misurato")
    # Perfect scaling halves the time at each doubling: a straight line of slope
    # -1 through the first point, which log-log axes make easy to read off.
    ax.plot(threads, [mean[0] * base / n for n in threads], "k--", alpha=0.6,
            label="ideale")
    ax.set_xlabel("thread")
    ax.set_ylabel("tempo del solve [s]")
    ax.set_title("Tempo")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")

    ax = axes[1]
    ax.plot(threads, stats["speedup"], "o-", label="misurato")
    ax.plot(threads, ideal, "k--", alpha=0.6, label="ideale")
    ax.set_xlabel("thread")
    ax.set_ylabel("speedup")
    ax.set_title("Speedup")
    ax.set_xscale("log", base=2)

    ax = axes[2]
    ax.plot(threads, stats["efficiency"], "o-")
    ax.axhline(1.0, color="k", ls="--", alpha=0.6)
    ax.set_xlabel("thread")
    ax.set_ylabel("efficienza")
    ax.set_title("Efficienza")
    ax.set_xscale("log", base=2)
    # Superlinear speedup is real (the smaller per-thread working set starts
    # fitting in cache), so a hardcoded 1.1 ceiling would clip those points off
    # the plot without saying so.
    ax.set_ylim(0, max(1.1, stats["efficiency"].max() * 1.05))

    for ax in axes:
        # Powers of two on a log axis default to 2^k labels; the plain thread
        # counts are what the reader is looking for.
        ax.set_xticks(threads)
        ax.set_xticklabels([str(t) for t in threads])
        ax.grid(True, which="both", alpha=0.2)
        if ax is not axes[2]:
            ax.legend()

    fig.suptitle("Strong scaling -- griglia %s, %d ripetizioni"
                 % (stats["size"].iloc[0], int(stats["runs"].max())))
    fig.tight_layout()

    fig.savefig(out_path, dpi=150)
    print("scritto %s" % out_path)
    if show:
        plt.show()


def main():
    ap = argparse.ArgumentParser(
        description="Grafico di strong scaling dai CSV del profiler.")
    ap.add_argument(
        "paths",
        nargs="*",
        help="file CSV, glob o cartelle da leggere (default: results). "
             "Le cartelle sono esplorate ricorsivamente.",
    )
    ap.add_argument(
        "--csv-name",
        default="profiling.csv",
        help="nome (o glob) dei CSV da cercare dentro le cartelle "
             "(default: profiling.csv)",
    )
    ap.add_argument(
        "-o",
        "--output",
        default="results/strong_scaling.png",
        help="immagine da scrivere (default: results/strong_scaling.png)",
    )
    ap.add_argument(
        "--csv-out",
        help="salva anche la tabella riassuntiva in questo CSV",
    )
    ap.add_argument(
        "--show",
        action="store_true",
        help="apre la finestra invece di limitarsi a salvare",
    )
    args = ap.parse_args()

    paths = expand(args.paths or ["results"], args.csv_name)
    if not paths:
        sys.exit("nessun file trovato in: %s" % ", ".join(args.paths or ["results"]))

    df, n_files = load(paths)
    if df.empty:
        sys.exit("nessuna riga %s leggibile in %d file "
                 "(il binario e' stato compilato con LBM_ENABLE_PROFILING=ON?)"
                 % (LABEL, len(paths)))

    stats, base = summarize(df)
    print("%d ripetizioni lette (%d file esaminati)" % (n_files, len(paths)))
    print_table(stats)
    if args.csv_out:
        stats.to_csv(args.csv_out)
        print("scritto %s" % args.csv_out)
    plot(stats, base, args.output, args.show)


if __name__ == "__main__":
    main()
