#!/bin/bash

sample="$1"

root_dir=/work-zfs/abattle4/heyuan/Variant_calling/datasets/GEO/PRJNA484801
cd ${root_dir}/ATAC_seq/alignment_bowtie/VCF_files/GRCh37/
VCF=${sample}.filtered.recode.GRCh37.vcf

ml htslib
ml vcftools

bgzip $VCF #compress vcf

VCF=${sample}.filtered.recode.GRCh37.vcf.gz
tabix -p vcf ${VCF} # index compressed vcf


mkdir -p ${sample}
tabix -H ${VCF} > ${sample}/${VCF}_header.txt
rm -f ${sample}/${sample}_chr*
for chromosome in {1..22}
do
	outfn=${sample}/${sample}_chr${chromosome}.vcf
	cat ${sample}/${VCF}_header.txt > ${outfn}
        tabix ${VCF} chr${chromosome} | sed 's/chr//g' | sed '/	<NON_REF>/d' | sed 's/,<NON_REF>//g' >> ${outfn}
done

cd ${sample}
ml bcftools
for chromosome in {1..22}
do
	outfn=${sample}_chr${chromosome}.vcf
	bgzip ${outfn}
        tabix -p vcf ${outfn}.gz
	bcftools annotate -x INFO ${outfn}.gz |  bcftools annotate -x FORMAT > temp
	rm ${outfn}.gz*
	mv temp ${outfn}
	bgzip ${outfn}
	tabix -p vcf ${outfn}.gz
done
