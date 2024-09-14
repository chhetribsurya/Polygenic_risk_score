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
workdir_argv = sys.argv[3]

#input filenames:
tissue = "Cells_EBV-transformed_lymphocytes"
#work_dir = "/home-4/schhetr1@jhu.edu/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps"
#work_dir = "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/high_ldsnps"
work_dir = workdir_argv

#annotpath = join(work_dir, "ldsc_annot")
annotpath = annotpath_argv

#prefix="EBV_LCL_Binary_100kb_colocthresh075_baselineLD"
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

    df_list = []
    annotation_list = ["PP.H4.75"]
    #annotation_list = ["PP.H1.75", "PP.H2.75", "PP.H3.75", "PP.H4.75"]
    
    for annotation in annotation_list:
        coloc_geneqtl_file = join(work_dir, "{}.coloc.abf.topsnpList_forLDSC.txt".format(annotation))
        #coloc_geneqtl_file = join(work_dir, "{}.coloc.abf.topsnpList_forLDSC.txt.ldSnps.txt".format(annotation))

        print("\nprocessing chrom:{} annotation {} ...".format(chrom, annotation))
        print('generating annot file...')

        #load bimfile and snp file for merge
        df_bim = pd.read_csv(args_bimfile,
                delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])

        df_int = pd.read_csv(coloc_geneqtl_file, sep="\t", header=None)
        df_int.columns = ["SNP"]
        df_int.drop_duplicates(["SNP"], inplace=True)
        df_int["ANNOT"] =1

        df_annot = pd.merge(df_bim, df_int, how='left', on='SNP')
        df_annot.fillna(0, inplace=True)
        df_annot["ANNOT"] = df_annot[['ANNOT']].astype(int)
        arrange_cols = ["CHR", "BP", "SNP", "CM", "ANNOT"]
        df_annot = df_annot.loc[:,arrange_cols]
        df_annot.rename(columns={"ANNOT":annotation}, inplace=True)
        df_list.append(df_annot)

    df_merge_idx = df_list[0]
    for df in df_list[1:]:
        df_merge_idx = pd.merge(df_merge_idx, df, how='left', on=["CHR", "BP", "SNP", "CM"])

    df_merged = df_merge_idx.fillna(0)
    df_merged["base"] = 1
    select_cols = ["CHR", "BP", "SNP", "CM"] + ["base"] + annotation_list
    df_merged_f = df_merged.loc[:,select_cols]
    df_merged_f.to_csv(args_annot_file, sep = "\t", index = False, compression="gzip")
    #print("\n\nmerging custom annot with baselineLD annotation ...\n\n")  
    #df_merged["base"] = 1
    #select_cols = ["CHR", "BP", "SNP", "CM"] + ["base"] + annotation_list
    #df_merged_f = df_merged.loc[:,select_cols]
    #baselineLD_dir = "/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/baselineLD_annot"
    #df_baselineLD = pd.read_csv(join(baselineLD_dir, "baselineLD.{}.annot.gz".format(chrom)), sep="\t") 
    #df_merged_baselineLD = pd.merge(df_baselineLD, df_merged, how='left', on=["CHR", "BP", "SNP", "CM"])
    #df_merged_baselineLD.to_csv(args_annot_file, sep = "\t", index = False, compression="gzip")

print("Task completed...")

#call step2 script
import subprocess
print("\nRunning step2 shell script\n")
CMD_RUN="/work-zfs/abattle4/surya/datasets/prs_project/scripts/sldsc_pipes/snp_based/run_step2.sh {} {}".format(prefix, annotpath)
subprocess.check_call(CMD_RUN, shell=True)

