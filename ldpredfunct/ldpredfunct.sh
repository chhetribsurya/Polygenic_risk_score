#!/usr/bin/bash
cohort="EUR"
workdir="/work-zfs/abattle4/surya/datasets/WHIMS"
outputdir=$workdir/ldpredfunct_output
mkdir -p $outputdir; 

#PHENO_EUR="/work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/WHIMS/race5/unrelated/ph    eno_mendrand.txt"
PHENO_AA="/work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/PAGE/race3/unrelated/pheno_mendrand.txt"

plinkfile="$workdir/$cohort/whimsUnrelated_${cohort}_[1:22]"; 
statsfile="$workdir/UKBB_sumstats/SBP_summarystats_baseline_ldpredFormat.txt"; 
phenotype="$workdir/phenotype_${cohort}.txt"; 
#functfile="$workdir/functfile_allsnp_baseline.txt"; 
functfile="$workdir/functfile_eQTL_baseline.txt"; 
functfile="/work-zfs/abattle4/surya/tools/LDpred-funct/test/test_functfile.txt"; 

outCoord="$outputdir/Coord_Final"; 
outLdpredfunct="$outputdir/ldpredfunct_posterior_means"; 
outValidate="$outputdir/ldpredfunct_prs"; 

N=373000; 
h2=0.1512; 
K=100

#source activate ldsc-new
eval "$(conda shell.bash hook)"
conda activate ldpred_funct_new

cd /work-zfs/abattle4/surya/tools/LDpred-funct
python /work-zfs/abattle4/surya/tools/LDpred-funct/ldpredfunct.py \
    --gf=${plinkfile} \
    --pf=${phenotype} \ 
    --ssf=${statsfile} \
    --FUNCT_FILE=${functfile} \
    --N=${N} \
    --H2=${h2} \ 
    --K=${K} \
    --coord=${outCoord} \ 
    --posterior_means=${outLdpredfunct} \  
     --out=${outValidate} > ${outValidate}.log
