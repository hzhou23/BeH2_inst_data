#!/bin/bash
#SBATCH -N 1
#SBATCH -C knl
#SBATCH -q low
#SBATCH -t 00:10:00
#SBATCH -o hf.o

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

module load gamess

#run the application:
rungms cisd_4.inp 00 8 &> cisd_4.out

