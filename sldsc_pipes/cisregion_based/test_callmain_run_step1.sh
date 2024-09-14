#!/usr/bin/nash
export PARTITION="skylake"
export CORE=6
export MEMORY="20G"
export TIME=24:00:00

#fileinfo
export OUTDIR="/work-zfs/abattle4/surya/datasets/prs_project/scripts/sldsc_pipes/cisregion_based"
export LOGDIR="$OUTDIR/logfiles"
mkdir -p $LOGDIR

export JOBNAME="sLDSC"
export OUTFILE=$LOGDIR/${JOBNAME}_annotat_h2est_combineRes_%N_%J.out
export ERRORFILE=$LOGDIR/${JOBNAME}_annotat_h2est_combineRes_%N_%J.err

#notification
export MAILTYPE="END,FAIL"
export MAILUSER="chhetribsurya@gmail.com"

#scriptinfo
export SCRIPTPATH="/work-zfs/abattle4/surya/datasets/prs_project/scripts/sldsc_pipes/cisregion_based"
#export SCRIPTFILE1="step1.2_baselineLDannot_generate_genebased_ldscAnnotation.py"
#export SCRIPTFILE2="step2.1_baselineLD_batchrun_w_summarystats_forH2estimate.sh"
#export SCRIPTFILE3="step3.1_baselineLD_concat_all_ldsc_resultfiles.R"

#source activate ldsc-new
eval "$(conda shell.bash hook)"
conda activate ldsc-new

#FILELIST=""
#for FILE_PREFIX in $FILELIST; do 
#    echo "processing: job $each ..."; 
#    #sbatch -p $PARTITION -c $CORE -t $TIME --mem=$MEMORY -J "$JOBNAME" -o "$OUTFILE" -e "$ERRORFILE" slurm-test/slurm_test1.sh ${FILE_PREFIX};
#    sbatch -p $PARTITION -c $CORE -t $TIME --mem=$MEMORY -J "$JOBNAME" -o "$OUTFILE" -e "$ERRORFILE" --mail-type "$MAILTYPE" --mail-user "$MAILUSER" $SCRIPTPATH/$SCRIPTFILE1 ${FILE_PREFIX};
#    #sbatch -p $PARTITION -c $CORE -t $TIME --mem=$MEMORY -J "$JOBNAME" -o "$OUTFILE" -e "$ERRORFILE" --mail-type "$MAILTYPE" --mail-user "$MAILUSER" $SCRIPTPATH/$SCRIPTFILE2 ${FILE_PREFIX};
#    #sbatch -p $PARTITION -c $CORE -t $TIME --mem=$MEMORY -J "$JOBNAME" -o "$OUTFILE" -e "$ERRORFILE" --mail-type "$MAILTYPE" --mail-user "$MAILUSER" $SCRIPTPATH/$SCRIPTFILE3 ${FILE_PREFIX};
#done

export WORKDIR="/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps"
export ANNOTPATH="$WORKDIR/ldsc_annot/cisregion_based"

export PREFIX0="EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_WITHbaseONLY"
export PREFIX1="EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_baselineLDONLY"
export PREFIX2="EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_baselineLDwithPPH0to4"
export PREFIX3="EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_baselineLDwithPPH3plus4"
export PREFIX4="EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_baselineLDwithPPH4"

export TESTING="TRUEVAL"

#heritability enrichment with base only
sbatch -p $PARTITION -c $CORE -t $TIME --mem=$MEMORY -J "${JOBNAME}_0" -o "${OUTFILE}.0" -e "${ERRORFILE}.0" --mail-type "$MAILTYPE" --mail-user "$MAILUSER" $SCRIPTPATH/step1_WITHbaseONLY.sh ${PREFIX0} ${ANNOTPATH};

#heritability enrichment with baseline LD only
#sbatch -p $PARTITION -c $CORE -t $TIME --mem=$MEMORY -J "${JOBNAME}_1" -o "${OUTFILE}.1" -e "${ERRORFILE}.1" --mail-type "$MAILTYPE" --mail-user "$MAILUSER" $SCRIPTPATH/step1_baselineLDONLY.sh ${PREFIX1} ${ANNOTPATH};

#heritability enrichment with PPH0, PPH1, PPH2, PPH3 and PPH4
#sbatch -p $PARTITION -c $CORE -t $TIME --mem=$MEMORY -J "${JOBNAME}_2" -o "${OUTFILE}.2" -e "${ERRORFILE}.2" --mail-type "$MAILTYPE" --mail-user "$MAILUSER" $SCRIPTPATH/step1_baselineLDwithPPH0to4.py ${PREFIX2} ${ANNOTPATH} ${WORKDIR};

#heritability enrichment with aggregated PPH3 and PPH4
#sbatch -p $PARTITION -c $CORE -t $TIME --mem=$MEMORY -J "${JOBNAME}_3" -o "${OUTFILE}.3" -e "${ERRORFILE}.3" --mail-type "$MAILTYPE" --mail-user "$MAILUSER" $SCRIPTPATH/step1_baselineLDwithPPH3plus4.py ${PREFIX3} ${ANNOTPATH} ${WORKDIR};

#heritability enrichment with PPH4 only
#sbatch -p $PARTITION -c $CORE -t $TIME --mem=$MEMORY -J "${JOBNAME}_4" -o "${OUTFILE}.4" -e "${ERRORFILE}.4" --mail-type "$MAILTYPE" --mail-user "$MAILUSER" $SCRIPTPATH/step1_baselineLDwithPPH4.py ${PREFIX4} ${ANNOTPATH} ${WORKDIR};



