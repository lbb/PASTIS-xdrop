
#!/bin/bash
rm pastis-*.log ./build/sim_mat.mtx slurm.ipu* 
cd build/ && make -j$(nproc) && cd .. && sbatch runpastis