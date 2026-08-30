#!/bin/bash
# Run it FROM THE LOGIN NODE with: bash submit_weak.sh [N]
N_JOBS=${1:-10}

for i in $(seq 1 "$N_JOBS"); do
  JOBID=$(qsub run_weak.pbs)
  echo "  run $i -> ${JOBID}"
done

echo "Fatto. Controlla con: qstat -u \$(whoami)"
