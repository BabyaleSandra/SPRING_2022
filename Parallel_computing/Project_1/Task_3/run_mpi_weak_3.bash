#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH -N 3
#SBATCH -n 64
#SBATCH --partition=shortq
#SBATCH --job-name=std_dev
#SBATCH --output=std_dev.o%j

module load gcc
module load mpich
module load slurm
#source ~/.bashrc

p=(1 2 4 8 16 32 64)
echo "nproc,N,std_dev,elapsed_time" >> task3_weak.csv
for i in 0 1 2 3 4 5 6 
do 
mpirun -np ${p[i]} ./task_3 $((2**(22+$i))) >> task3_weak.csv
done