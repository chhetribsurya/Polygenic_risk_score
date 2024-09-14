#!/usr/bin/bash

#activate ldpred-funct
eval "$(conda shell.bash hook)"
conda activate ldpred_funct_new

#args variables
FUNCTFILE=$1
N=$2
h2=$3

COHORT=$4
WHIMDIR=$5
OUTDIR=$6
FILETYPE=$7

#sanity check
echo "Following parameters provided..."
echo -e "\n check params functfile: $FUNCTFILE"
echo -e "\n check params N: $N"
echo -e "\n check params h2: $h2"
echo -e "\n check params cohort: $COHORT"
echo -e "\n check params whimdir: $WHIMDIR"
echo -e "\n check params outdir: $OUTDIR"
echo -e "\n check params filetype: $FILETYPE"

#input params for ldpredfunct script
plinkfile="$WHIMDIR/${COHORT}_dedup/whimsUnrelated_${COHORT}_[1:22]"; 
statsfile="$WHIMDIR/UKBB_sumstats/SBP_Unique_summarystats_baseline_ldpredFormat.txt"; 
phenotype="$WHIMDIR/phenotype_${COHORT}.txt"; 

functfile="$FUNCTFILE"; 

outCoord="$OUTDIR/Coord_Final${FILETYPE}"; 
outLdpredfunct="$OUTDIR/ldpredfunct_posterior_means_${FILETYPE}"; 
outValidate="$OUTDIR/ldpredfunct_PRS_${FILETYPE}"; 

outLdpredfunct="$OUTDIR/ldpredfunct_posterior_means_${FILETYPE}";
echo -e "\nplinkfile: $plinkfile"
echo -e "\nphenotype file: $phenotype"
echo -e "\noutfile: $outLdpredfunct" 
echo -e "\nprocessing file debug: ${outValidate}.log\n"

#N=340159
#h2=0.1197
K=10

# Note: Make sure that the file ${outCoord}  does not exists already.
if [[ -f $outCoord ]];then rm $outCoord; fi
#cd /work-zfs/abattle4/surya/tools/LDpred-funct
python /home/schhetr1/scratch16-abattle4/surya/tools/LDpred-funct/ldpredfunct.py \
    --gf=${plinkfile} \
    --pf=${phenotype} \
    --ssf=${statsfile} --FUNCT_FILE=${functfile} --N=${N} \
    --H2=${h2} \
    --K=${K} \
    --coord=${outCoord} \
    --posterior_means=${outLdpredfunct} \
    --out=${outValidate} > ${outValidate}.log

    #python /work-zfs/abattle4/surya/tools/LDpred-funct/ldpredfunct.py --gf=${plinkfile} --pf=${phenotype} --FUNCT_FILE=${functfile}  --coord=${outCoord} --ssf=${statsfile} --N=${N} --posterior_means=${outLdpredfunct}  --H2=${h2} --out=${outValidate} > ${outValidate}.log

