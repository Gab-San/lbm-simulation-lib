#!/bin/bash
# Run FROM THE LOGIN NODE with: bash submit_properties.sh [N]
set -euo pipefail

N_JOBS=${1:-10}

# Build once, here, BEFORE submitting the jobs. If this is left to
# run_properties.pbs instead, the N_JOBS jobs start almost simultaneously,
# all see the binary missing, and all launch "cmake --build" in parallel
# against the same shared build/ (NFS, not scratch) -> a corrupted/failed
# build for most of the jobs, with no visible error other than in the job
# logs.
if [ ! -x build/simulations/properties_test ]; then
  echo "Binary not found: building before submitting the jobs..."
  (cd build && cmake --build . --target properties_test -j"$(nproc)")
fi

for i in $(seq 1 "$N_JOBS"); do
  JOBID=$(qsub run_properties.pbs)
  echo "  run $i -> ${JOBID}"
done

echo "Done. Check with: qstat -u \$(whoami)"
