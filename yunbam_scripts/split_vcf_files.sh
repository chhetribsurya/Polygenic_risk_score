#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=3:00:00

chr="$1"
vcftools --gzvcf batch1.vcf.bgz --chr chr${chr} --maf 0.01 --out batch1.chr${chr} --recode

