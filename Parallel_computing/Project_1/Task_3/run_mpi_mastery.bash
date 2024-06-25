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


echo "nproc,std_dev,gather_comm_time,reduce_comm_time,communication_time,comput_timer1,comput_timer2,computational_time,elapsed_time" >> task3_mastery.csv
for p in 1 2 4 8 16 32 64 
do 
mpirun -np $p ./task_3 $((2**28)) >> task3_mastery.csv
done