#!/bin/bash

EXE=$1
CONFIG=$2
N_JOBS=${3:-1}
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

mkdir -p logs

if [ ! -x ${EXE} ]; then
  echo "Binary not found: building before submitting the jobs..."
  (cd build && cmake --build . --target properties_test -j"$(nproc)")
  exit 1;
fi

filename="${EXE##*/}"   # Removes path -> my.document.txt
stem="${filename%.*}"    # Removes extension -> my.document

if [[ $stem == cuda_* ]]; then
  for i in $(seq 1 "$N_JOBS"); do
    JOBID=$(qsub -v "EXE=${EXE},CONFIG=${CONFIG}" "$SCRIPT_DIR/cuda_job.pbs")
    echo "  run $i -> ${JOBID}"
  done
else
  for i in $(seq 1 "$N_JOBS"); do
    JOBID=$(qsub -v "EXE=${EXE},CONFIG=${CONFIG}" "$SCRIPT_DIR/cpu_job.pbs")
    echo "  run $i -> ${JOBID}"
  done
fi

echo "Done. Check with: qstat -u \$(whoami)"
echo "Check where it is running with: pbsnoodes -aSj"
echo "Logs will appear in: $(pwd)/logs/lbm_properties.o<jobid>"
