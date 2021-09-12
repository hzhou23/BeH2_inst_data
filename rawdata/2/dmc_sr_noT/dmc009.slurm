#!/bin/bash
#SBATCH -N 64
#SBATCH -C knl
#SBATCH -q overrun
#SBATCH -t 01:00:00
#SBATCH -o dmc0001

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
srun -n 17408 -c 1 --cpu_bind=cores qwalk-cori BeH2009.ldmc
