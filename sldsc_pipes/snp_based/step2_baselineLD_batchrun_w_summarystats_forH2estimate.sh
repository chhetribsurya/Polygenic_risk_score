#!/usr/bin/bash

####################
##run s-LDSC program:
####################

#source activate ldsc-new
eval "$(conda shell.bash hook)"
conda activate ldsc-new
 
prefix_argv=$1
annotpath_argv=$2

#scriptpath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc"
scriptpath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc"

#bfile_path="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_EUR_Phase3_plink"
bfile_path="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_EUR_Phase3_plink"

#hapmap3_snpspath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/hapmap3_snps"
hapmap3_snpspath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/hapmap3_snps"

#freqpath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_Phase3_frq"
freqpath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_Phase3_frq"

#regweights_path="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_Phase3_weights_hm3_no_MHC"
regweights_path="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_Phase3_weights_hm3_no_MHC"

##################
#parameter changes
annotpath=$annotpath_argv
prefix=$prefix_argv

annotation_dir="${annotpath}/${prefix}"

genome_1kg_prefix="1000G.EUR.QC"

hm3snp_prefix="weights.hm3_noMHC"


#compute ldscore estimates
for chr in {1..22}; do \
    python ${scriptpath}/ldsc.py \
      --l2 \
      --bfile ${bfile_path}/1000G.EUR.QC.${chr} \
      --ld-wind-cm 1 \
      --print-snps ${hapmap3_snpspath}/hm.${chr}.snp \
      --annot ${annotation_dir}/${prefix}.${chr}.annot.gz \
      --out ${annotation_dir}/${prefix}.${chr}
done

#compute snp-heritability estimates for traits
sumstats_path="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/summary_stats_ldsc_format"
#sumstats_path="/work-zfs/abattle4/surya/datasets/prs_project/UKBB/highly_heritable_traits/ldsc_sumstats"
echo -e "\n\nStarting H2estimate batchrun for summarystats files ...\n\n"

#batch run with sumstats file list:
#for sumstatsfile in ${sumstats_path}/*cancer*.bgz; do
#for sumstatsfile in ${sumstats_path}/*.bgz; do
for sumstatsfile in ${sumstats_path}/*blood*; do

    filename=$(basename ${sumstatsfile} | cut -f1 -d "."); #without pheno_id
    filename=$(basename ${sumstatsfile} | cut -f1-2 -d "."); #with pheno_id
    outputdir="${annotation_dir}/results_${filename}";
    if [[ ! -d ${outputdir} ]];then mkdir -p ${outputdir}; fi

    #compute proportion of heritability explained by snps
    python ${scriptpath}/ldsc.py \
        --h2 ${sumstatsfile} \
        --ref-ld-chr ${annotation_dir}/${prefix}. \
        --frqfile-chr ${freqpath}/${genome_1kg_prefix}. \
        --w-ld-chr ${regweights_path}/${hm3snp_prefix}. \
        --overlap-annot \
        --print-coefficients \
        --print-delete-vals \
        --out ${outputdir}/${prefix}.ldsc;

done

#compute snp-heritability estimates for highly heritable traits
#sumstats_path="/work-zfs/abattle4/surya/datasets/prs_project/UKBB/highly_heritable_traits/ldsc_sumstats"
#echo -e "\n\nStarting H2estimate batchrun for summarystats files ...\n\n"
#for sumstatsfile in ${sumstats_path}/*.bgz; do

#    filename=$(basename ${sumstatsfile} | cut -f1 -d "."); #without pheno_id
#    filename=$(basename ${sumstatsfile} | cut -f1-2 -d "."); #with pheno_id
#    #outputdir="${annotation_dir}/results_${filename}";
#    outputdir="${annotation_dir}/highly_heritableTrait_results_${filename}";
#    if [[ ! -d ${outputdir} ]];then mkdir -p ${outputdir}; fi
#
#    #compute proportion of heritability explained by snps
#    python ${scriptpath}/ldsc.py \
#        --h2 ${sumstatsfile} \
#        --ref-ld-chr ${annotation_dir}/${prefix}. \
#        --frqfile-chr ${freqpath}/${genome_1kg_prefix}. \
#        --w-ld-chr ${regweights_path}/${hm3snp_prefix}. \
#        --overlap-annot \
#        --print-coefficients \
#        --print-delete-vals \
#        --out ${outputdir}/${prefix}.ldsc;
#
#done

echo -e "\nTask completed ..."
#echo -e "\nCheck result dir: ${outputdir}\n"

#run step3 script
bash /work-zfs/abattle4/surya/datasets/prs_project/scripts/sldsc_pipes/snp_based/run_step3.sh $prefix $annotpath

