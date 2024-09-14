#!/usr/bin/bash
export PARTITION="skylake"
export CORE=14
export MEMORY="45G"
export TIME=30:00:00
export MAILTYPE="END,FAIL"
export MAILUSER="chhetribsurya@gmail.com"

#input params
export COHORT="AA"
export WHIMDIR="/work-zfs/abattle4/surya/datasets/WHIMS"
N=340159
h2=0.1197

#fileinfo
export WORKDIR="/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"
export OUTDIR="$WORKDIR/HighMEM45_LDpredfunct_Output_${COHORT}"
export LOGDIR="$OUTDIR/Logfiles"

mkdir -p $OUTDIR
mkdir -p $LOGDIR

#scriptinfo
export SCRIPTPATH="/work-zfs/abattle4/surya/datasets/prs_project/scripts/ldpredfunct"

#activate ldpred-funct
eval "$(conda shell.bash hook)"
conda activate ldpred_funct_new

#setup args variable
#functfile=$1
#N=$2
#h2=$3
#outdir=$4 #$WORKDIR
#outfilename=$5 #FILENAME+TRAIT

FUNCTFILE_LIST="/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based/EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_*/funct_file/*"
#FUNCTFILE_LIST="/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based/EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_WITHbaseONLYwithPPH4/funct_file/systolic_bloodPressure.4080_irnt.ldpredfunct_file.txt"

#PHENO_EUR="/work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/WHIMS/race5/unrelated/pheno_mendrand.txt"
#PHENO_AA="/work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/PAGE/race3/unrelated/pheno_mendrand.txt"

for FUNCTFILE in ${FUNCTFILE_LIST}; do

    #echo -e "\n\nProcessing $(basename $FUNCTFILE) ...\n\n"
    PREFIX=$(basename $(dirname $(dirname $FUNCTFILE)))
    TRAIT_ID=$(basename $FUNCTFILE | awk -F ".ldpredfunct" '{print $1}')
    FUNCT_ID=$(echo $(basename $PREFIX) | awk -F "colocthresh075_" '{print $2}')
    FILETYPE=${TRAIT_ID}_${FUNCT_ID}
    echo -e "\nProcessing ${FILETYPE}\n" 

    JOBNAME="LDpred45_${COHORT}_${FUNCT_ID}"
    OUTFILE=$LOGDIR/${FILETYPE}_${COHORT}_%N_%J.out
    ERRORFILE=$LOGDIR/${FILETYPE}_${COHORT}_%N_%J.err

    #$SCRIPTPATH/new_ldpredfunct_AA.sh ${FUNCTFILE} ${N} ${h2} ${COHORT} ${WHIMDIR} ${OUTDIR} ${FILETYPE};
    sbatch -p $PARTITION -c $CORE -t $TIME --mem=$MEMORY -J "${JOBNAME}" -o "${OUTFILE}" -e "${ERRORFILE}" --mail-type "$MAILTYPE" --mail-user "$MAILUSER" $SCRIPTPATH/new_ldpredfunct_AA.sh ${FUNCTFILE} ${N} ${h2} ${COHORT} ${WHIMDIR} ${OUTDIR} ${FILETYPE};

done
