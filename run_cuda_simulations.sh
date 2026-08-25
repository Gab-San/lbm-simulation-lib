#!/usr/bin/env bash
#
# run_cuda_simulations.sh - Run CUDA simulation executables sequentially
#
# Usage:
#   ./run_cuda_simulations.sh sim1 sim2 sim3
#   ./run_cuda_simulations.sh          # runs all cuda-* executables found in BIN_DIR

# Optional environment variables:
#   BIN_DIR=/path/to/executables   (default: ./build/simulations)
#   LOG_DIR=/path/to/logs          (default: ./logs)

set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${BIN_DIR:-$ROOT_DIR/build/simulations}"
LOG_DIR="${LOG_DIR:-$ROOT_DIR/logs}"
mkdir -p "$LOG_DIR" "$ROOT_DIR/out"

# If simulation names are passed as args, use those; otherwise auto-discover
if [ "$#" -gt 0 ]; then
    SIMULATIONS=("$@")
else
    SIMULATIONS=()
    while IFS= read -r exe; do
        name="$(basename "$exe")"
        SIMULATIONS+=("${name#cuda_}")
    done < <(find "$BIN_DIR" -maxdepth 1 -type f -name 'cuda_*' -perm -u+x | sort)
fi

if [ "${#SIMULATIONS[@]}" -eq 0 ]; then
    echo "No CUDA simulations found or specified in $BIN_DIR" >&2
    echo "Build first with: cmake -S . -B build -DLBM_ENABLE_CUDA=ON && cmake --build build -j2" >&2
    exit 1
fi

echo "Running ${#SIMULATIONS[@]} CUDA simulation(s): ${SIMULATIONS[*]}"
echo "Executables: $BIN_DIR"
echo "Logs:       $LOG_DIR"
echo "----------------------------------------"

FAILED=()

for sim in "${SIMULATIONS[@]}"; do
    # Accept either the short name or the complete cuda_* target name.
    if [[ "$sim" == cuda_* ]]; then
        exe_name="$sim"
        short_name="${sim#cuda_}"
    else
        exe_name="cuda_$sim"
        short_name="$sim"
    fi

    exe="$BIN_DIR/$exe_name"
    log_file="$LOG_DIR/${exe_name}.log"

    if [ ! -x "$exe" ]; then
        echo "[SKIP] $exe not found or not executable"
        FAILED+=("$short_name (missing)")
        continue
    fi

    echo "[START] $exe_name"
    start_time=$(date +%s)

    if (cd "$ROOT_DIR" && "$exe") >"$log_file" 2>&1; then
        elapsed=$(( $(date +%s) - start_time ))
        echo "[DONE]  $exe_name (${elapsed}s) -> $log_file"
    else
        status=$?
        elapsed=$(( $(date +%s) - start_time ))
        echo "[FAIL]  $exe_name (${elapsed}s, exit $status) -> $log_file"
        FAILED+=("$short_name")
    fi
done

echo "----------------------------------------"
if [ "${#FAILED[@]}" -eq 0 ]; then
    echo "All simulations completed successfully."
else
    echo "Completed with failures: ${FAILED[*]}" >&2
    exit 1
fi
