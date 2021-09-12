#!/bin/bash
#SBATCH -A m2113
#SBATCH -q low
#SBATCH -C knl
#SBATCH -J vmcnoj_multi
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=256
#SBATCH -o vmc.out
#SBATCH -e vmc.err
#SBATCH --export=ALL

echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR


export HDF5_USE_FILE_LOCKING=FALSE
module unload cray-libsci
module load cray-hdf5
#module load python/3.7-anaconda-2019.10
source activate myenv
module list

export OMP_NUM_THREADS=256
srun qmcpack cisd_4.qmc.in-wfnoj.xml

