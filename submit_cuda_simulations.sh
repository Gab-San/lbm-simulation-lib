#!/bin/bash

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${BIN_DIR:-$ROOT_DIR/build/cuda/simulations}"
LOG_DIR="${LOG_DIR:-$ROOT_DIR/logs}"
mkdir -p "$LOG_DIR" "$ROOT_DIR/out"

echo "${ROOT_DIR}"
echo "${BIN_DIR}"
echo "${LOG_DIR}"

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

echo "Submitting ${#SIMULATIONS[@]} CUDA simulation(s): ${SIMULATIONS[*]}"
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

    # Actually submit to the scheduler with qsub -F. Running the .pbs file
    # directly with "sh" (as before) just executes it inline on the login
    # node - the #PBS lines are inert comments to sh, so nothing is queued
    # and no GPU is ever reserved.
    if JOBID=$(qsub -F "$exe_name" submit_cuda_job.pbs); then
        echo "  $short_name -> ${JOBID}" | tee -a "$log_file"
    else
        echo "  $short_name -> qsub failed" | tee -a "$log_file"
        FAILED+=("$short_name (qsub failed)")
    fi
done

echo "----------------------------------------"
if [ "${#FAILED[@]}" -eq 0 ]; then
    echo "All simulations submitted successfully."
else
    echo "Completed with failures: ${FAILED[*]}" >&2
    exit 1
fi
