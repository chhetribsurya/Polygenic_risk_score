#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00

ml vcftools
chr="$1"

filename=out.chr${chr}.maf0.05.recode.vcf.gz
vcftools --gzvcf ${filename} --max-alleles 2 --out final.biallelic.maf005.chr${chr} --recode

output=final.biallelic.maf005.chr${chr}.recode.vcf
bgzip ${output}
tabix -p vcf ${output}.gz

#rm out.chr${chr}.maf0.05.recode.vcf
