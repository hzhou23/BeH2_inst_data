#!/bin/bash
#SBATCH -N 8
#SBATCH -C knl
#SBATCH -q regular
#SBATCH -t 15:00:00
#SBATCH -o dmc001

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
srun -n 256 -c 1 --cpu_bind=cores python atomic_nexus.py

