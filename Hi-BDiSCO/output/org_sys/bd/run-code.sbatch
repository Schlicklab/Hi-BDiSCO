#!/bin/bash
#SBATCH --job-name=cuda
#SBATCH --gres=gpu
#SBATCH --partition=v100
#SBATCH --nodes=1
#SBATCH --time=24:00:00

module purge
module load openmpi/intel/4.1.1
module load cuda/11.3.1

./code 
