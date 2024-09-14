#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00

ml vcftools
chr="$1"

vcftools --gzvcf merged.chr${chr}.recode.vcf.gz --keep samples_final.txt | bgzip -c > merged.chr${chr}.recode.filterSamples.vcf.gz 
vcftools --gzvcf batch1.chr${chr}.recode.vcf.gz  --keep samples_final.txt | bgzip -c > batch1.chr${chr}.recode.filterSamples.vcf.gz
