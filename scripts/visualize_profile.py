import argparse
import sys
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
    ghia_data = []

    for line in f:
        line_str = line.strip()
        if not line_str:
            continue
        x, y = line_str.split()
        ghia_data.append((float(x), float(y)))

    # Normalizzo i dati di Ghia rispetto al massimo valore in modulo
    max_value_ghia = max([abs(y) for _, y in ghia_data])
    if max_value_ghia > 0:
        ghia_data = [(x, y / max_value_ghia) for x, y in ghia_data]

    return np.array(ghia_data)


def read_benchmark(f, args):
    """%%benchmark <name>\\n"""
    for arg in args:
        if arg == "ghia":
            return read_ghia_data(f)
        else:
            raise ValueError(f"Invalid benchmark {args}")


def read_profile(f, args):
    """%%profile <model> <nx|ny>\\n then raw float64 binary data.
    Payload is a 1D centerline (length nx or ny)."""
    data = np.fromfile(f, dtype=np.float64)

    expected = int(args[1])
    if data.size != expected:
        raise ValueError("Number of data inputs not matching data size.")

    # Normalizzo i dati LBM
    max_value = np.abs(data).max() if len(data) > 0 else 0
    if max_value > 0:
        data = data / max_value

    return data  # fallback: unknown shape, return flat


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
            else:
                raise ValueError(f"Unknown header kind '{kind}' in {file}")

        results[file] = (kind, args, data)

    return results


def plot_ghia_data(plt, ghia_data):
    plt.plot(ghia_data[:, 0], ghia_data[:, 1], label="Ghia et al. (1982)", marker="s")


def create_graph(parsed_results, title, xlabel, output_file=None):
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
        else:
            # data parte da x=0 a x=1 con passo 1/(len(data)-1)
            N = len(data)
            x = (np.arange(N) + 0.5) / N
            plt.plot(x, data, label=args[0])

    # Aggiungo i dati di Ghia et al. 1982 per confronto
    plt.xlabel(xlabel)
    plt.ylabel("u velocity")
    plt.title(title)
    plt.legend(fontsize="small", loc="best")
    plt.grid(True)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Grafico salvato in: {output_file}")

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
        default="Comparison of u-velocity profiles",
        help="Custom title",
    )

    parser.add_argument(
        "--xlabel",
        nargs="?",
        type=str,
        default="y coordinate",
        help="Custom xlabel",
    )

    parser.add_argument("-o", "--output", nargs="?", type=str, help="Output file path")

    args = parser.parse_args()

    return args


def main():
    args = parse_args()

    results = read_data(args.files)
    create_graph(results, title=args.title, xlabel=args.xlabel, output_file=args.output)


if __name__ == "__main__":
    main()
