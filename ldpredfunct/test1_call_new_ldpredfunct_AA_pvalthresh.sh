#!/usr/bin/bash
export PARTITION="skylake"
export CORE=22
export MEMORY="70G"
export TIME=20:00:00
export MAILTYPE="END,FAIL"
export MAILUSER="chhetribsurya@gmail.com"

#input params
export COHORT="AA"
export WHIMDIR="/work-zfs/abattle4/surya/datasets/WHIMS"
N=340159
h2=0.1197

#fileinfo
export WORKDIR="/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"
export OUTDIR="$WORKDIR/HighMEM45_LDpredfunct_Output_${COHORT}/corrected_pvalthresh_prs_1"
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
FUNCTFILE_LIST="/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based/EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_WITHbaseONLYwithPPH4/funct_file/systolic_bloodPressure.4080_irnt.ldpredfunct_file.txt"
#FUNCTFILE_LIST="/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based/EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_baselineLDwithPPH4/funct_file/systolic_bloodPressure.4080_irnt.ldpredfunct_file.txt"

#PHENO_EUR="/work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/WHIMS/race5/unrelated/pheno_mendrand.txt"
#PHENO_AA="/work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/PAGE/race3/unrelated/pheno_mendrand.txt"

PVAL_THRESH=(0.5 1)
for FUNCTFILE in ${FUNCTFILE_LIST}; do

    #echo -e "\n\nProcessing $(basename $FUNCTFILE) ...\n\n"
    PREFIX=$(basename $(dirname $(dirname $FUNCTFILE)))
    TRAIT_ID=$(basename $FUNCTFILE | awk -F ".ldpredfunct" '{print $1}')
    FUNCT_ID=$(echo $(basename $PREFIX) | awk -F "colocthresh075_" '{print $2}')
    FILETYPE=${TRAIT_ID}_${FUNCT_ID}
    echo -e "\nProcessing ${FILETYPE}\n"
    
    #FILENAME=${FUNCTFILE}.pvalthesh.txt
    #filter funct files base on their pvals
    for PVAL in ${PVAL_THRESH[@]};do 
        echo -e "\nProcessing pval thresh: $PVAL"; 
        awk -v var=$PVAL '$3 <= var {print}' ${FUNCTFILE} > ${FUNCTFILE}.pvalthresh${PVAL}.txt; 
        PVAL_FUNCTFILE=${FUNCTFILE}.pvalthresh${PVAL}.txt 
        PVAL_FILETYPE=${TRAIT_ID}_${FUNCT_ID}.pvalthresh${PVAL}.txt

        JOBNAME="LDpred45_${COHORT}_${FUNCT_ID}.pvalthresh${PVAL}"
        OUTFILE=$LOGDIR/${FILETYPE}_${COHORT}.pvalthresh${PVAL}_%N_%J.out
        ERRORFILE=$LOGDIR/${FILETYPE}_${COHORT}.pvalthresh${PVAL}_%N_%J.err
        #echo -e "outfile:$OUTFILE"
        #echo -e "errorfile:$ERRORFILE"

        #$SCRIPTPATH/new_ldpredfunct_AA.sh ${PVAL_FUNCTFILE} ${N} ${h2} ${COHORT} ${WHIMDIR} ${OUTDIR} ${PVAL_FILETYPE};
        sbatch -p $PARTITION -c $CORE -t $TIME --mem=$MEMORY -J "${JOBNAME}" -o "${OUTFILE}" -e "${ERRORFILE}" --mail-type "$MAILTYPE" --mail-user "$MAILUSER" $SCRIPTPATH/new_ldpredfunct_${COHORT}.sh ${PVAL_FUNCTFILE} ${N} ${h2} ${COHORT} ${WHIMDIR} ${OUTDIR} ${PVAL_FILETYPE};
    done


done
