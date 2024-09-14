#!/usr/bin/env python

#################################
##prepare gene based annot file for ldsc run:
#################################
from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool
import gzip
import os, sys
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

def gene_set_to_bed(geneset_file, genecoord_file, windowsize):
    print('generating gene set bed file...')
    GeneSet = pd.read_csv(geneset_file, header = None, names = ['GENE'])
    all_genes = pd.read_csv(genecoord_file, delim_whitespace = True)
    df = pd.merge(GeneSet, all_genes, on = 'GENE', how = 'inner')
    df['START'] = np.maximum(1, df['START'] - args_windowsize)
    df['END'] = df['END'] + args_windowsize
    iter_df = [['chr'+(str(x1).lstrip('chr')), x2 - 1, x3] for (x1,x2,x3) in np.array(df[['CHR', 'START', 'END']])]
    bed_for_annot = BedTool(iter_df).sort().merge()
    #bed_for_annot = BedTool(iter_df).sort()
    return(bed_for_annot)


"""Input from args annotation_list"""
args_gene_coord_file = "/home-4/schhetr1@jhu.edu/surya/datasets/prs_project/ENSG_coord.txt"
#args_gene_coord_file = "/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/make_annot_sample_files/ENSG_coord.txt"

args_windowsize = 500000

####################################################
#prepare annotation files for gene based input file:
####################################################

chrom_list = []
chromrange=22
for chrom_no in range(chromrange):
    chrom_list.append(str(chrom_no + 1))

#chrom_list = ["20", "21", "22"]
for chrom in chrom_list:
    #chrom=22
    #args_bimfile="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_EUR_Phase3_plink/1000G.EUR.QC.{}.bim".format(chrom)
    args_bimfile="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_EUR_Phase3_plink/1000G.EUR.QC.{}.bim".format(chrom)

    outfilename= prefix + ".{}.annot.gz".format(chrom)
    args_annot_file=join(annotation_dir, outfilename)

    df_list = []
    
    #annotation_list = ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]
    annotation_list = ["PP.H0.75", "PP.H1.75", "PP.H2.75", "PP.H3.75", "PP.H4.75"]
    # annotation_list = ["PP.H0.75", "PP.H1.75", "PP.H2.75", "PP.H3.75", "PP.H4.75", "PP.H4.8", "PP.H4.9", "PP.H4.95"]
    
    for annotation in annotation_list:
        #annotation = "Null"
        args_gene_set_file = join(work_dir, "{}.coloc.abf.txt".format(annotation))
        #args_gene_set_file = "/home-4/schhetr1@jhu.edu/surya/datasets/prs_project/coloc_output/Coloc_Geneset_EBV_LCC_{}.abf.txt".format(fileidx)

        print("\nprocessing chrom:{} annotation {} ...".format(chrom, annotation))
        bed_for_annot = gene_set_to_bed(args_gene_set_file, args_gene_coord_file, args_windowsize)

        print('generating annot file...')
        df_bim = pd.read_csv(args_bimfile,
                delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
        iter_bim = [['chr'+str(x1), int(x2) - 1, int(x2)] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
        bimbed = BedTool(iter_bim)
        annotbed = bimbed.intersect(bed_for_annot)

        bp = [x.start + 1 for x in annotbed]
        #df_int = pd.DataFrame({'BP': bp, 'ANNOT':1})
        df_int = pd.DataFrame({'BP': bp})
        df_int["ANNOT"] = 1

        df_annot = pd.merge(df_bim, df_int, how='left', on='BP')
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
CMD_RUN="/work-zfs/abattle4/surya/datasets/prs_project/scripts/sldsc_pipes/cisregion_based/run_step2.sh {} {}".format(prefix, annotpath)
subprocess.check_call(CMD_RUN, shell=True)
