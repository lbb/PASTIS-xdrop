
#!/bin/bash
set -euo pipefail
pushd build/ && make -j$(nproc) && popd

set +euo pipefail
scancel $(squeue -n runpastis -h | awk '{print $1}')
rm pastis-*.log ./build/sim_mat.mtx 
rm slurm.log
touch slurm.log

set -euo pipefail
sbatch runpastis

sleep 1

squeue --me

tail -f slurm.log