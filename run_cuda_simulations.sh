#!/usr/bin/env bash
#
# run_cuda_simulations.sh - Run CUDA simulation executables sequentially
#
# Usage:
#   ./run_cuda_simulations.sh sim1 sim2 sim3
#   ./run_cuda_simulations.sh          # runs all cuda-* executables found in BIN_DIR

set -uo pipefail

# Directory where your cuda-<simulation_name> executables live
BIN_DIR="./bin"

# Directory to store logs (one per simulation)
LOG_DIR="./logs"
mkdir -p "$LOG_DIR"

# If simulation names are passed as args, use those; otherwise auto-discover
if [ "$#" -gt 0 ]; then
    SIMULATIONS=("$@")
else
    SIMULATIONS=()
    for exe in "$BIN_DIR"/cuda-*; do
        [ -x "$exe" ] || continue
        name="$(basename "$exe")"
        SIMULATIONS+=("${name#cuda-}")
    done
fi

if [ "${#SIMULATIONS[@]}" -eq 0 ]; then
    echo "No simulations found or specified." >&2
    exit 1
fi

echo "Running ${#SIMULATIONS[@]} simulation(s): ${SIMULATIONS[*]}"
echo "----------------------------------------"

FAILED=()

for sim in "${SIMULATIONS[@]}"; do
    exe="$BIN_DIR/cuda-$sim"
    log_file="$LOG_DIR/${sim}.log"

    if [ ! -x "$exe" ]; then
        echo "[SKIP] $exe not found or not executable"
        FAILED+=("$sim (missing)")
        continue
    fi

    echo "[START] cuda-$sim"
    start_time=$(date +%s)

    if "$exe" > "$log_file" 2>&1; then
        elapsed=$(( $(date +%s) - start_time ))
        echo "[DONE]  cuda-$sim (${elapsed}s) -> $log_file"
    else
        elapsed=$(( $(date +%s) - start_time ))
        echo "[FAIL]  cuda-$sim (${elapsed}s) -> see $log_file"
        FAILED+=("$sim")
    fi
done

echo "----------------------------------------"
if [ "${#FAILED[@]}" -eq 0 ]; then
    echo "All simulations completed successfully."
else
    echo "Completed with failures: ${FAILED[*]}"
    exit 1
fi
