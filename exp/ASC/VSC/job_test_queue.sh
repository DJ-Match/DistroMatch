#!/bin/bash

#SBATCH -J scaling
#SBATCH -A p72884
#SBATCH -t 00:10:00
#SBATCH --partition=zen3_0512
#SBATCH --qos=zen3_0512_devel
#SBATCH --nodes=4               # Request nodes
#SBATCH --ntasks=4              # MPI tasks total
#SBATCH --ntasks-per-node=1     # 1 task per node
#SBATCH --cpus-per-task=1       # 1 core per task
#SBATCH --gres=gpu:0 

# module purge
# module load ...

#execute your program
./run.sh run_s24_p1 1

#uname -a
#hostnamectl

