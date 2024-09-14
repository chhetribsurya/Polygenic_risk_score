#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH -p skylake

sample="$1"
root_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/PRJNA484801
cd ${root_dir}/ATAC_seq/alignment_bowtie/VCF_files/

VCF=${sample}.filtered.recode.vcf.gz
 
source ~/.bash_profile
ml gcc/5.5.0

save_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Imputation/${sample}
mkdir -p ${save_dir}
ml vcftools

for chromosome in {1..22}
do
	refPanel=/work-zfs/abattle4/lab_data/imputation_reference_panel/${chromosome}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz
	vcfFn=GRCh37/${sample}/${sample}_chr${chromosome}.vcf.gz

	outFn=GRCh37/${sample}/${sample}_chr${chromosome}.imputed.GRCh37
	minimac4 --refHaps ${refPanel} \
         --haps ${vcfFn} \
         --prefix ${outFn} \
 	 --ChunkLengthMb 100

done

cd GRCh37/${sample}
for chromosome in {1..22}
do
	outFn=${sample}_chr${chromosome}.imputed.GRCh37
	gunzip ${outFn}.dose.vcf.gz
	vcf-to-tab < ${outFn}.dose.vcf > ${save_dir}/chr${chromosome}.imputed.GRCh37.genotype.txt
done
