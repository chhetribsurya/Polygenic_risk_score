#!/bin/bash -l
#SBATCH
#SBATCH --job-name=simulation
#SBATCH --partition=gpup100
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
#SBATCH --mem=30GB
#SBATCH --time=12:00:00
#SBATCH --mail-user=marvani1@jhu.edu
#SBATCH --mail-type=ALL
#SBATCH --account=abattle4

cd /work-zfs/abattle4/marios/BP_MR/750k_DBP_new_prs

module load gcc/5.5.0
module load R/3.5.1

Rscript ../run_simulation_fracpoly.r


