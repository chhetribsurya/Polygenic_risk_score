#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=12:00:00
#SBATCH -p skylake

sample="$1"
root_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/PRJNA484801
cd ${root_dir}/ATAC_seq/alignment_bowtie/VCF_files/
mkdir -p GRCh37/

VCF=${sample}.filtered.recode.vcf
outputFn=GRCh37/${sample}.filtered.recode.GRCh37.vcf

VCF=${VCF/.gz/}
cat ${VCF} |  grep "^#" > ${VCF}_temp
cat ${VCF} |  grep -v "^#" | awk '{print "chr"$0}' >> ${VCF}_temp

ml picard
picard LiftoverVcf I=${VCF}_temp O=${outputFn} CHAIN=/work-zfs/abattle4/heyuan/tools/liftOver/hg38ToHg19.over.chain.gz REJECT=${outputFn/.vcf/.rejected.vcf} R=/work-zfs/abattle4/heyuan/database/GRCh37_reference/hg19.fa
rm ${VCF}_temp


ml htslib
bgzip ${VCF}
