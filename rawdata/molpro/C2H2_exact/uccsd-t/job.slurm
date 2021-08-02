#!/bin/bash -l
#SBATCH -J exact
#SBATCH -q debug
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -t 00:05:00
#SBATCH -o slurm.out
#!SBATCH -o slurm.o%j

module load molpro

#molpro -n 12 -s 2.com
#molpro -n 12 -s 3.com
molpro -n 12 -s 4.com
#molpro -n 12 -s 5.com
#molpro -n 12 -s 6.com

