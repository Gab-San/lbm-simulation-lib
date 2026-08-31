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
                print("warning: no file for %s" % arg, file=sys.stderr)
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
        print("warning: %s could not be parsed (%s)" % (path, err), file=sys.stderr)
        return None

    if list(df.columns) != FIELDS:
        print("warning: unexpected header in %s, file skipped" % path,
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
        print("warning: row dropped in %s: %s" % (path, row.to_dict()),
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
    dup = df.duplicated(["source", "size", "n_threads"], keep=False)
    if dup.any():
        # Two solve_total rows for the same size/thread count in one file means
        # the binary was run twice into the same CSV; averaging them is fine,
        # but a silent average would hide a botched job script.
        print("warning: repeated solve_total rows for the same grid/n_threads in %s"
              % ", ".join(sorted(df[dup]["source"].unique())), file=sys.stderr)

    stats = df.groupby(["size", "n_threads"]).agg(
        runs=("total", "size"),
        mean_s=("total", "mean"),
        std_s=("total", "std"),
    ).reset_index()
    stats[["mean_s", "std_s"]] /= 1000.0
    # std() is NaN on a single repetition, which measures no spread. 0 draws no
    # error bar, which is honest; the run count is printed in the table below.
    stats["std_s"] = stats["std_s"].fillna(0.0)

    # Each grid size scales on its own: speedup/efficiency are measured
    # against the smallest thread count present for THAT size, not against a
    # hardcoded 1 or another size's base.
    stats["base"] = stats.groupby("size")["n_threads"].transform("min")
    base_mean = stats[stats["n_threads"] == stats["base"]].set_index("size")["mean_s"]
    stats["speedup"] = base_mean.reindex(stats["size"]).values / stats["mean_s"]
    stats["efficiency"] = stats["speedup"] / (stats["n_threads"] / stats["base"])
    return stats.sort_values(["size", "n_threads"])


TABLE_COLUMNS = ["size", "n_threads", "runs", "mean_s", "std_s", "speedup", "efficiency"]
TABLE_FORMATTERS = {
    "mean_s": "{:.3f}".format,
    "std_s": "{:.3f}".format,
    "speedup": "{:.3f}".format,
    "efficiency": "{:.3f}".format,
}


def format_table(stats):
    """The summary table with numeric columns pre-formatted as strings, so the
    text printout and the image export show exactly the same numbers."""
    out = stats[TABLE_COLUMNS].copy()
    for col, fmt in TABLE_FORMATTERS.items():
        out[col] = out[col].map(fmt)
    return out


def print_table(stats):
    print(format_table(stats).to_string(index=False))


def table_image_path(out_path):
    """<output>_table.<ext>, keeping png/pdf and falling back to png otherwise."""
    base, ext = os.path.splitext(out_path)
    if ext.lower() not in (".png", ".pdf"):
        ext = ".png"
    return base + "_table" + ext


def save_table_image(stats, path):
    # Reuses whichever backend plot() already selected (Agg unless --show),
    # so this must run after plot() rather than choosing its own.
    import matplotlib.pyplot as plt

    df = format_table(stats)
    fig, ax = plt.subplots(figsize=(max(6, 1.1 * len(df.columns)),
                                    0.35 * (len(df) + 1) + 0.6))
    ax.axis("off")
    tbl = ax.table(cellText=df.values, colLabels=df.columns, loc="center",
                   cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.4)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote %s" % path)


def plot(stats, out_path, show):
    import matplotlib

    if not show:
        # No DISPLAY on a compute or login node: pick the file-only backend
        # before pyplot is imported, or the import itself fails there.
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sizes = stats["size"].unique()
    colors = plt.cm.tab10.colors
    all_threads = sorted(stats["n_threads"].unique())

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    for i, size in enumerate(sizes):
        g = stats[stats["size"] == size]
        threads = g["n_threads"].to_list()
        mean = g["mean_s"].to_list()
        base = int(g["base"].iloc[0])
        c = colors[i % len(colors)]

        ax = axes[0]
        ax.errorbar(threads, mean, yerr=g["std_s"], fmt="o-", capsize=4,
                    color=c, label=str(size))
        # Perfect scaling halves the time at each doubling: a straight line of
        # slope -1 through the first point, which log-log axes make easy to
        # read off.
        ax.plot(threads, [mean[0] * base / n for n in threads], "--", color=c,
                alpha=0.5)

        ax = axes[1]
        ax.plot(threads, g["speedup"], "o-", color=c, label=str(size))
        ax.plot(threads, [n / base for n in threads], "--", color=c, alpha=0.5)

        ax = axes[2]
        ax.plot(threads, g["efficiency"], "o-", color=c, label=str(size))

    axes[0].set_xlabel("threads")
    axes[0].set_ylabel("solve time [s]")
    axes[0].set_title("Time")
    axes[0].set_yscale("log")

    axes[1].set_xlabel("threads")
    axes[1].set_ylabel("speedup")
    axes[1].set_title("Speedup")

    axes[2].set_xlabel("threads")
    axes[2].set_ylabel("efficiency")
    axes[2].set_title("Efficiency")
    axes[2].axhline(1.0, color="k", ls="--", alpha=0.6)
    # Superlinear speedup is real (the smaller per-thread working set starts
    # fitting in cache), so a hardcoded 1.1 ceiling would clip those points off
    # the plot without saying so.
    axes[2].set_ylim(0, max(1.1, stats["efficiency"].max() * 1.05))

    for ax in axes:
        # Powers of two on a log axis default to 2^k labels; the plain thread
        # counts are what the reader is looking for.
        ax.set_xscale("log", base=2)
        ax.set_xticks(all_threads)
        ax.set_xticklabels([str(t) for t in all_threads])
        ax.set_xlim(all_threads[0] / 1.4, all_threads[-1] * 1.4)
        ax.grid(True, which="both", alpha=0.2)
        ax.legend(title="grid" if len(sizes) > 1 else None)

    n_reps = int(stats["runs"].max())
    if len(sizes) > 1:
        fig.suptitle("Strong scaling -- %d grids, up to %d repetitions"
                     % (len(sizes), n_reps))
    else:
        fig.suptitle("Strong scaling -- grid %s, %d repetitions"
                     % (sizes[0], n_reps))
    fig.tight_layout()

    fig.savefig(out_path, dpi=150)
    print("wrote %s" % out_path)
    if show:
        plt.show()


def main():
    ap = argparse.ArgumentParser(
        description="Plot the strong scaling curve from the profiler CSVs.")
    ap.add_argument(
        "paths",
        nargs="*",
        help="CSV files, globs or directories to read (default: results). "
             "Directories are walked recursively.",
    )
    ap.add_argument(
        "--csv-name",
        default="profiling.csv",
        help="name (or glob) of the CSVs to look for inside directories "
             "(default: profiling.csv)",
    )
    ap.add_argument(
        "-o",
        "--output",
        default="results/strong_scaling.png",
        help="image file to write (default: results/strong_scaling.png)",
    )
    ap.add_argument(
        "--csv-out",
        help="also save the summary table to this CSV",
    )
    ap.add_argument(
        "--show",
        action="store_true",
        help="open a window instead of only saving",
    )
    args = ap.parse_args()

    paths = expand(args.paths or ["results"], args.csv_name)
    if not paths:
        sys.exit("no files found in: %s" % ", ".join(args.paths or ["results"]))

    df, n_files = load(paths)
    if df.empty:
        sys.exit("no readable %s rows in %d files "
                 "(was the binary compiled with LBM_ENABLE_PROFILING=ON?)"
                 % (LABEL, len(paths)))

    stats = summarize(df)
    print("%d repetitions read (%d files scanned)" % (n_files, len(paths)))
    print_table(stats)
    if args.csv_out:
        stats.to_csv(args.csv_out, index=False)
        print("wrote %s" % args.csv_out)
    plot(stats, args.output, args.show)
    save_table_image(stats, table_image_path(args.output))


if __name__ == "__main__":
    main()
