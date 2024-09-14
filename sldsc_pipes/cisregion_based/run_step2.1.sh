#!/usr/bin/bash

#compute instance
export PARTITION="skylake"
export CORE=8
export MEMORY="25G"
export TIME=24:00:00

#fileinfo
export OUTDIR="/work-zfs/abattle4/surya/datasets/prs_project/scripts/sldsc_pipes/cisregion_based"
export JOBNAME="sLDSC"
export OUTFILE=$OUTDIR/${JOBNAME}_step2heritability_%N_%J.out
export ERRORFILE=$OUTDIR/${JOBNAME}_step2heritability_%N_%J.out

#notification
export MAILTYPE="END,FAIL"
export MAILUSER="chhetribsurya@gmail.com"

#script info
export SCRIPTPATH="/work-zfs/abattle4/surya/datasets/prs_project/scripts/sldsc_pipes/cisregion_based"

#args variable
export PREFIX=$1
export ANNOTPATH=$2

#sbatch run
#sbatch -p $PARTITION -c $CORE -t $TIME --mem=$MEMORY -J "$JOBNAME" -o "$OUTFILE" -e "$ERRORFILE" --mail-type "$MAILTYPE" --mail-user "$MAILUSER" $SCRIPTPATH/step2_baselineLD_batchrun_w_summarystats_forH2estimate.sh ${PREFIX} ${ANNOTPATH};

$SCRIPTPATH/step2.1_baselineLD_batchrun_w_summarystats_forH2estimate.sh ${PREFIX} ${ANNOTPATH};

echo "partition again: $PARTITION"
echo "time again: $TIME"
echo -e "\tcore used again: $CORE"
echo -e "\ttesting export again: $TESTING"
