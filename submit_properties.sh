#!/bin/bash
# Lanciare DAL LOGIN NODE con: bash submit_properties.sh [N]
N_JOBS=${1:-10}

for i in $(seq 1 "$N_JOBS"); do
  JOBID=$(qsub run_properties.pbs)
  echo "  run $i -> ${JOBID}"
done

echo "Fatto. Controlla con: qstat -u \$(whoami)"
