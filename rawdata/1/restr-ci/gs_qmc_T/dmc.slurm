#!/bin/bash
#SBATCH -A m2113
#SBATCH -q flex --time-min=02:00:00 --time=04:00:00 
#SBATCH -C knl
#SBATCH -J BeH2_2_ref_T
#SBATCH -N 16
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=272
#SBATCH -o 001.multi.out
#SBATCH -e 001.multi.err
#SBATCH --export=ALL

echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR


export HDF5_USE_FILE_LOCKING=FALSE
module unload cray-libsci
module load cray-hdf5
#module load python/3.7-anaconda-2019.10
source activate myenv
module list

export OMP_NUM_THREADS=1
 srun qmcpack gs_qmc.xml
~                            
