#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --time=24:00:00
#SBATCH -p skylake


root_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/ENA/ATAC_seq/alignment_bowtie
ml fastqtl
ml bcftools
ml htslib

ncovs="$1"
chromosome="$2"
sample_file=/work-zfs/abattle4/heyuan/Variant_calling/QTL_calling/samples_for_final_peaks_exist.txt

peakDir=${root_dir}/Peaks_MACS2_clean
peakFn=${peakDir}/peaks_chr${chromosome}.txt
#rm -f  ${peakFn}.gz
#bgzip ${peakFn}&& tabix -p bed ${peakFn}.gz

cov=${peakDir}/Derived_PCs_uniqueIDs_v3.txt
VCFfn=/work-zfs/abattle4/heyuan/Variant_calling/datasets/ENA/Gencove/chr${chromosome}.maf005.biallelic.recode.peakSamples.vcf.gz

outdir=${root_dir}/fastQTL
mkdir -p ${outdir}

for cisDis in 100000
do
	fastqtl --vcf ${VCFfn} --bed ${peakFn}.gz --cov ${cov} --include-samples ${sample_file} --include-covariates covs/covs_${ncovs}.txt --window ${cisDis} --region chr${chromosome}:1-300000000 --out ${outdir}/chr${chromosome}.window${cisDis}.cov${ncovs}.fastq.results
	fastqtl --vcf ${VCFfn} --bed ${peakFn}.gz --cov ${cov} --include-samples ${sample_file} --include-covariates covs/covs_${ncovs}.txt --window ${cisDis}  --region chr${chromosome}:1-300000000 --permute 1000 --out ${outdir}/chr${chromosome}.window${cisDis}.cov${ncovs}.fastq.permutation.results
done

