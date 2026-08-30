#!/bin/bash
# Run FROM THE LOGIN NODE with: bash submit_properties.sh [N]
set -euo pipefail

N_JOBS=${1:-10}

# run_properties.pbs writes its merged stdout/stderr log to logs/ (see
# "#PBS -o logs/" in that file). PBS does not create the output directory
# for you - if it's missing, some sites silently drop the log instead of
# erroring, which is indistinguishable from "job ran and did nothing".
mkdir -p logs

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
echo "Logs will appear in: $(pwd)/logs/lbm_properties.o<jobid>"
