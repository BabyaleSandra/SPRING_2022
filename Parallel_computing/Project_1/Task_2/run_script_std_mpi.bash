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


echo "nproc,std_dev,end_timer_broadcast,end_timer_allreduce,end_timer_reduce,communication_time,computational_time,elapsed_timer" >> t1strong.csv
for p in 1 2 4 8 16 32 64 
do 
mpirun -np $p ./task_1 $((2**28)) >> t1strong.csv
done