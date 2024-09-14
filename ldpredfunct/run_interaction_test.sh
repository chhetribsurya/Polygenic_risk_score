#!/bin/bash
#SBATCH --job-name=BP_interaction
#SBATCH --time=72:0:0
#SBATCH --partition=lrgmem
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --mem=50GB

Rscript run_interaction_test.r
