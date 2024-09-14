#!/usr/bin/bash

#compute instance
export PARTITION="skylake"
export CORE=4
export MEMORY="10G"
export TIME=6:00:00

#fileinfo
export OUTDIR="/work-zfs/abattle4/surya/datasets/prs_project/scripts/sldsc_pipes/snp_based"
export JOBNAME="sLDSC"
export OUTFILE=$OUTDIR/${JOBNAME}_step3combineResults_%N_%J.out
export ERRORFILE=$OUTDIR/${JOBNAME}_step3combineResults_%N_%J.out

#notification
export MAILTYPE="END,FAIL"
export MAILUSER="chhetribsurya@gmail.com"

#script info
export SCRIPTPATH="/work-zfs/abattle4/surya/datasets/prs_project/scripts/sldsc_pipes/snp_based"

#args variable
export PREFIX=$1
export ANNOTPATH=$2

#source activate ldsc-new
eval "$(conda shell.bash hook)"
conda activate r4base

#sbatch run
#sbatch -p $PARTITION -c $CORE -t $TIME --mem=$MEMORY -J "$JOBNAME" -o "$OUTFILE" -e "$ERRORFILE" --mail-type "$MAILTYPE" --mail-user "$MAILUSER" Rscript $SCRIPTPATH/step3_baselineLD_concat_all_ldsc_resultfiles.R ${PREFIX} ${ANNOTPATH};

Rscript $SCRIPTPATH/step3_concat_all_ldsc_resultfiles.R ${PREFIX} ${ANNOTPATH};
