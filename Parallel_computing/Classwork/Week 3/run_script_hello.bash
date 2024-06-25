#!/bin/bash
###
###
#SBATCH --time=01:00:00
#SBATCH --tasks=28
#SBATCH --partition=shortq
#SBATCH --job-name=hello
#SBATCH --output=hello.o%j

module load gcc
module load openmpi/gcc

mpirun -np 28 ./hello_mpi
