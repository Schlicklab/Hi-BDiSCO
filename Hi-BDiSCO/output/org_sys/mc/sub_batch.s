#!/bin/bash
#
#SBATCH --job-name=gene_lh
#SBATCH --nodes=4
#SBATCH --time=48:00:00
#SBATCH --tasks-per-node=24
#SBATCH --mem=50GB
 
module purge
module load openmpi/intel/4.1.1

srun ./chrom_ap1.x dim.in input.1 LH_N0G6C22.in LH_N0G6C22_equil.in constraints.dat >> log_ch
