#!/usr/bin/env python3
"""Validate LBM output files structurally and numerically.

Checks:
- every output path declared by simulations/cuda/*.cu exists (unless --allow-missing)
- no two CUDA sources declare the same output path
- profile files have a valid %%profile header, expected payload length and finite float64 data
- norms files contain nx, ny int32 header, an integral number of float32 frames and finite data
"""

from __future__ import annotations

import argparse
import re
import struct
from collections import defaultdict
from pathlib import Path

import numpy as np


def declared_cuda_outputs(repo: Path):
    occurrences: dict[str, list[str]] = defaultdict(list)
    for src in sorted((repo / "simulations" / "cuda").glob("*.cu")):
        text = src.read_text(errors="ignore")
        for rel in re.findall(r'"(out/[^"\\]+\.bin)"', text):
            occurrences[rel].append(src.name)
    return occurrences


def validate_profile(path: Path):
    with path.open("rb") as f:
        line = f.readline().decode("ascii", errors="strict").strip()
        parts = line.split()
        if len(parts) != 3 or parts[0] != "%%profile":
            raise ValueError(f"invalid profile header: {line!r}")
        model = parts[1]
        expected = int(parts[2])
        data = np.fromfile(f, dtype=np.float64)
    if data.size != expected:
        raise ValueError(f"profile payload has {data.size} values, expected {expected}")
    if not np.isfinite(data).all():
        raise ValueError("profile contains NaN/Inf")
    return f"profile {model}, n={data.size}, min={data.min():.6g}, max={data.max():.6g}"


def validate_norms(path: Path):
    with path.open("rb") as f:
        raw = f.read(8)
        if len(raw) != 8:
            raise ValueError("norm file shorter than 8-byte grid header")
        nx, ny = struct.unpack("=ii", raw)
        if nx <= 0 or ny <= 0:
            raise ValueError(f"invalid grid dimensions {nx}x{ny}")
        data = np.fromfile(f, dtype=np.float32)
    area = nx * ny
    if data.size == 0 or data.size % area != 0:
        raise ValueError(
            f"norm payload has {data.size} floats, not a positive multiple of grid area {area}"
        )
    frames = data.size // area
    if not np.isfinite(data).all():
        raise ValueError("norm data contains NaN/Inf")
    return f"norms {nx}x{ny}, frames={frames}, min={data.min():.6g}, max={data.max():.6g}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--allow-missing", action="store_true")
    args = ap.parse_args()
    repo = args.repo.resolve()

    declared = declared_cuda_outputs(repo)
    duplicates = {p: srcs for p, srcs in declared.items() if len(srcs) > 1}
    failed = False

    print(f"Declared CUDA output files: {len(declared)}")
    if duplicates:
        failed = True
        print("\n[ERROR] Duplicate output paths in CUDA sources:")
        for rel, srcs in sorted(duplicates.items()):
            print(f"  {rel}: {', '.join(srcs)}")

    print("\nOutput validation:")
    for rel in sorted(declared):
        path = repo / rel
        if not path.exists():
            tag = "WARN" if args.allow_missing else "ERROR"
            print(f"[{tag}] {rel}: missing")
            failed |= not args.allow_missing
            continue
        try:
            if path.name.startswith("data_"):
                details = validate_profile(path)
            elif path.name.startswith("norms_"):
                details = validate_norms(path)
            else:
                raise ValueError("unknown output file naming convention")
            print(f"[ OK ] {rel}: {details}")
        except Exception as exc:
            failed = True
            print(f"[ERROR] {rel}: {exc}")

    if failed:
        raise SystemExit(1)
    print("\nAll declared CUDA outputs passed validation.")


if __name__ == "__main__":
    main()
