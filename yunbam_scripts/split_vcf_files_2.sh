#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00

chr="$1"
vcftools --gzvcf Hapmap.recode.vcf.gz --chr chr${chr} --maf 0.01 --out Hapmap.chr${chr} --recode

