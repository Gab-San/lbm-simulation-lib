import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_header(f, path):
    """Read one '%%kind ...' line, return (kind, [args]). Leaves cursor after the line."""
    line = f.readline().decode("ascii", errors="ignore").rstrip("\n")
    if not line.startswith("%%"):
        raise ValueError(f"File {path} has no header!")
    parts = line[2:].split()
    kind, args = parts[0], parts[1:]
    return kind, args


def read_ghia_data(f):
    # Ghia's tables are already normalized by the lid velocity (u/U_lid),
    # so we return them as-is. Rescaling by their own max (as before)
    # distorts the true profile shape and silently hides magnitude
    # mismatches between the simulation and the benchmark.
    ghia_data = []

    for line in f:
        line_str = line.strip()
        if not line_str:
            continue
        x, y = line_str.split()
        ghia_data.append((float(x), float(y)))

    return np.array(ghia_data)


def read_benchmark(f, args):
    """%%benchmark <name>\\n"""
    for arg in args:
        if arg == "ghia":
            return read_ghia_data(f)
        else:
            raise ValueError(f"Invalid benchmark {args}")


def read_profile(f, args):
    """%%profile <model> <n> <u_ref>\\n then raw float64 binary data.
    Payload is a 1D centerline (length n).

    Normalized by u_ref (same convention Ghia's tables and
    compute_ghia_error() use) rather than by its own max, so the plotted
    curve is directly comparable in magnitude, not just shape."""
    data = np.fromfile(f, dtype=np.float64)

    expected = int(args[1])
    if data.size != expected:
        raise ValueError("Number of data inputs not matching data size.")

    if len(args) < 3:
        raise ValueError(
            "Profile header is missing u_ref (old file format) -- "
            "regenerate this file with the updated output(); the plotted "
            "magnitude cannot be trusted otherwise."
        )

    u_ref = float(args[2])
    if u_ref == 0.0:
        raise ValueError("u_ref in header is zero, cannot normalize.")

    return data / u_ref


def read_exact(f, args):
    """%%exact <name> <npoints>\\n then interleaved (x, u) float64 pairs.

    Coordinates are written by the simulation already normalized to [0, 1]
    with the same cell-centered convention as the profile output, so they
    collocate exactly with the simulation nodes and the vertical gap at
    each marker is the pointwise error."""
    data = np.fromfile(f, dtype=np.float64)

    expected = 2 * int(args[1])
    if data.size != expected:
        raise ValueError("Number of data inputs not matching data size.")

    return data.reshape(-1, 2)


def read_data(files):
    """Parse each file's header and dispatch to the right reader.
    Returns a dict: {path: (kind, args, data)}."""
    results = {}
    for file in files:
        with open(file, "rb") as f:
            kind, args = parse_header(f, file)

            if kind == "benchmark":
                text_f = (line.decode("ascii", errors="ignore") for line in f)
                data = read_benchmark(text_f, args)
            elif kind == "profile":
                data = read_profile(f, args)
            elif kind == "exact":
                data = read_exact(f, args)
            else:
                raise ValueError(f"Unknown header kind '{kind}' in {file}")

        results[file] = (kind, args, data)

    return results


def plot_ghia_data(plt, ghia_data):
    plt.plot(ghia_data[:, 0], ghia_data[:, 1], label="Ghia et al. (1982)", marker="s")


def create_graph(parsed_results, title, xlabel, ylabel, output_file=None):
    # creo un grafico per i valori di data e lo mostro
    plt.figure()

    for file, (kind, args, data) in parsed_results.items():
        if kind == "benchmark":
            match args[0]:
                case "ghia":
                    plot_ghia_data(plt, data)
                case _:
                    raise ValueError(
                        f"Plotting unknown benchmark data ({kind}, {args})"
                    )
        elif kind == "exact":
            plt.plot(
                data[:, 0],
                data[:, 1],
                linestyle="none",
                marker="o",
                markersize=4,
                fillstyle="none",
                markevery=max(1, len(data) // 40),
                label=f"{args[0]} (exact)",
            )
        elif kind == "profile":
            # data parte da x=0 a x=1 con passo 1/(len(data)-1)
            N = len(data)
            x = (np.arange(N) + 0.5) / N
            plt.plot(x, data, label=args[0])
        else:
            raise ValueError(f"Plotting unknown data kind '{kind}' ({args})")

    # Aggiungo i dati di Ghia et al. 1982 per confronto
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(fontsize="small", loc="best")
    plt.grid(True)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Plot saved in {output_file}")

    plt.show()


def existing_file(path_str: str) -> Path:
    """Validate that the path exists and is a regular file."""
    path = Path(path_str)
    if not path.exists():
        raise argparse.ArgumentTypeError(f"File '{path_str}' does not exist.")
    if not path.is_file():
        raise argparse.ArgumentTypeError(
            f"Path '{path_str}' is a directory, not a file."
        )
    return path


class VerboseHelpFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter
):
    pass


def parse_args():
    """Configure and parse command-line arguments"""

    description = """
    ===================================================================
                      FLUID PROFILE VISUALIZATION TOOL
    ===================================================================
    Visualize a fluid profile or compare multiple fluid profiles against one another.
        """

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=VerboseHelpFormatter,
    )

    parser.add_argument(
        "files",
        type=existing_file,
        metavar="PROFILE_DATA",
        nargs="+",
        help="Path to fluid profile.",
    )

    parser.add_argument(
        "--title",
        nargs="?",
        type=str,
        default="Velocity Profile",
        help="Custom title",
    )

    parser.add_argument(
        "--xlabel",
        nargs="?",
        type=str,
        default="x-coordinate evaluated at y/2",
        help="Custom xlabel",
    )

    parser.add_argument(
        "--ylabel",
        nargs="?",
        type=str,
        default="uy(x, y/2)",
        help="Custom xlabel",
    )

    parser.add_argument("-o", "--output", nargs="?", type=str, help="Output file path")

    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    results = read_data(args.files)
    create_graph(
        results,
        title=args.title,
        xlabel=args.xlabel,
        ylabel=args.ylabel,
        output_file=args.output,
    )


if __name__ == "__main__":
    main()
