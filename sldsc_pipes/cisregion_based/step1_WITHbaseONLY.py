#!/usr/bin/env python

#############################################
##prepare snp based annot file for ldsc run:
#############################################
from __future__ import print_function
import pandas as pd
import numpy as np
from pybedtools import BedTool
import gzip, os, sys, argparse
from os.path import join

prefix_argv = sys.argv[1]
annotpath_argv = sys.argv[2]

#input filenames:
tissue = "Cells_EBV-transformed_lymphocytes"
#work_dir = "/home-4/schhetr1@jhu.edu/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps"
#work_dir = "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/high_ldsnps"

#annotpath = join(work_dir, "ldsc_annot")
annotpath = annotpath_argv

#prefix="EBV_LCL_Binary_100kb_colocthresh075_WITHbaseONLY"
prefix = prefix_argv

annotation_dir=join(annotpath, prefix)
if not os.path.exists(annotation_dir):
    os.makedirs(annotation_dir)

#prepare annotation files for gene based input file:
chrom_list = []
chromrange=22
for chrom_no in range(chromrange):
    chrom_list.append(str(chrom_no + 1))

#chrom_list = ["20", "21", "22"]
for chrom in chrom_list:
    #chrom=11
    #args_bimfile="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_EUR_Phase3_plink/1000G.EUR.QC.{}.bim".format(chrom)
    args_bimfile="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_EUR_Phase3_plink/1000G.EUR.QC.{}.bim".format(chrom)

    outfilename= prefix + ".{}.annot.gz".format(chrom)
    args_annot_file=join(annotation_dir, outfilename)

    print("\nprocessing chrom:{} ...".format(chrom))
    print('generating annot file...')
 
    df_bim = pd.read_csv(args_bimfile,
                delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    df_bim["base"] = 1
    df_bim.to_csv(args_annot_file, sep = "\t", index = False, compression="gzip")

print("Task completed...")

#call step2 script
import subprocess
print("\nRunning step2 shell script\n")
CMD_RUN="/work-zfs/abattle4/surya/datasets/prs_project/scripts/sldsc_pipes/cisregion_based/run_step2.sh {} {}".format(prefix, annotpath)
subprocess.check_call(CMD_RUN, shell=True)

#sbatch -p $PARTITION -c $CORE -t $TIME --mem=$MEMORY -J "$JOBNAME" -o "$OUTFILE" -e "$ERRORFILE" --mail-type "$MAILTYPE" --mail-user "$MAILUSER" $SCRIPTPATH/$SCRIPTFILE1 ${PREFIX} ${ANNOTPATH};
