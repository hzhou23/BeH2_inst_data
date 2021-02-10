#!/usr/bin/env bash

source ~/anaconda3/etc/profile.d/conda.sh
conda activate py2
source /home/haihan/qp2/quantum_package.rc
export OMP_NUM_THREADS=4
filename='scf'


qp_convert_output_to_ezfio ${filename}.out -o ${filename}.ezfio > convert.out

qp set_mo_class --core="[]" --act="[1-31]" --inact="[32-150]" --del="[151-250]"  ${filename}.ezfio > active_space.out

qp_run fci ${filename}.ezfio &> ${filename}.fci

echo '1e-8' > ${filename}.ezfio/qmcpack/ci_threshold
qp_run truncate_wf_spin ${filename}.ezfio &> ${filename}.trunc
QP_STATE=1 qp_run save_for_qmcpack ${filename}.ezfio > ${filename}.dump
mv ${filename}.dump qmc.dump
mv QP2QMCACK.h5 qmc.h5

convert4qmc -QP qmc.dump -threshold 1e-20 -add3BodyJ -hdf5 -production
