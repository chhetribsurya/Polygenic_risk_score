#!/usr/bin/bash

prefix=$1
annotpath=$2

#prefix="EBV_LCL_Binary_100kb_colocthresh075_baselineLD_minusPPH0to4"
#mkdir -p "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/$prefix" 

mkdir -p ${annotpath}/${prefix}

#cp previous original baseline ldscore and annot files to new dir:
for each in /work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/baselineLD_annot/*;do
    filename="$prefix".$(echo $(basename $each) | cut -f2- -d"."); 
    echo cp ${each} ${annotpath}/${prefix}/${filename}; 
    echo processed $filename;
done

echo "partition: $PARTITION"
echo "time: $TIME"
echo -e "\tcore used: $CORE"
echo -e "\ttesting export: $TESTING"

#call step2 script
echo -e "\nRunning step2 shell script\n"
bash /work-zfs/abattle4/surya/datasets/prs_project/scripts/sldsc_pipes/snp_based/run_step2.1.sh $prefix $annotpath
