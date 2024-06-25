#!/bin/bash
###
###
#SBATCH --time=00:10:00
#SBATCH --tasks=16
#SBATCH --partition=shortq
#SBATCH --job-name=trapez
#SBATCH --output=trapez.o%j

module load gcc
module load mpich
module load slurm

mpirun -np 1 ./trapez_integration_mpi 1000000
mpirun -np 2 ./trapez_integration_mpi 1000000
mpirun -np 4 ./trapez_integration_mpi 1000000
mpirun -np 8 ./trapez_integration_mpi  1000000
mpirun -np 16 ./trapez_integration_mpi  1000000
