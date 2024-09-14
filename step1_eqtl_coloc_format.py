#!/usr/bin/env python
import numpy as np, pandas as pd
import seaborn as sns
import pybedtools
import os, re, pickle

from glob import glob
from os.path import basename, join, splitext
from collections import defaultdict
from functools import reduce 

import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy import stats

"""supress warnings"""
import time, warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None

""" Input files """
aa_filepath = "/Users/suryachhetri/datasets/prs_project/gtex/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_AFR_eQTL_all_associations"
ea_filepath = "/Users/suryachhetri/datasets/prs_project/gtex/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations"
gwas_stats = "/Users/suryachhetri/datasets/prs_project/gwas_stats/bmi/21001_irnt.gwas.imputed_v3.both_sexes.tsv"

#lifted_gwas_stats = "/Users/suryachhetri/datasets/prs_project/gwas_stats/bmi/lifted/updated_output.bed"

output_path = "/Users/suryachhetri/datasets/prs_project/final_hg19"
plot_path = join(output_path, "plots")
plot_path = join(output_path, "plots_chrX")
plot_path = join(output_path, "redo_chr1")

# Requires `liftOver`from UCSC to be on the path and a `chainfile`
chainfile = "/Users/suryachhetri/tools/chain_files/hg38ToHg19.over.chain"

""" create plot output dir"""
if not os.path.exists(plot_path):
    os.makedirs(plot_path)


def pop_shared_variants_colocformat(aa_filepath, ea_filepath, chromrange, pval, tissue, chrom_X=False):
    
    """ 
    Finds shared variants and specfic variants between 
    population cohorts without similar effects/eQTL info

    """

    ## Enter the chromosome no. range that you would like to analyse the data on. By default, it would take the autosomal 
    ## gene model from (chr1:chr22, excluding chrX & chrY), however, set chrom_XY=True for normal gene model analysis.
    #chromrange = 1

    chrom_list = []
    for chrom_no in range(chromrange):
        chrom_list.append("chr" + str(chrom_no + 1))
    if chrom_X:
        chrom_list.append("chrX")
        #chrom_list.append("chrY")

    # read chromwise for AFR and EUR:
    aa_filelist = glob(join(aa_filepath, str(tissue) + "*.parquet"))
    aa_dict = {basename(file).split(".")[4] :file for file in aa_filelist}
    ea_filelist = glob(join(ea_filepath, str(tissue) + "*.parquet"))
    ea_dict = {basename(file).split(".")[4] :file for file in ea_filelist}

    #eQTL_filt_stats = []
    pop_shared_betas_ses = [] 
    shared_cis_genes = []
    #genecount = []
    #variantcount = []

    for chrom in chrom_list:
        print("\nProcessing : {} ...\n".format(chrom))

        """ parse afr cohorts """
        read_aa = pd.read_parquet(aa_dict.get(chrom), engine="fastparquet")
        read_aa.reset_index(inplace=True, drop=True)
        aa_thresh = read_aa.groupby(["phenotype_id"]).apply(lambda X:  X[X["pval_nominal"] <= pval])
        #aa_df = read_aa[read_aa["phenotype_id"] =="ENSG00000001461.16"]

        """ parse eur cohorts"""
        read_ea = pd.read_parquet(ea_dict.get(chrom), engine="fastparquet")
        read_ea.reset_index(inplace=True, drop=True)
        ea_thresh = read_ea.groupby(["phenotype_id"]).apply(lambda X:  X[X["pval_nominal"] <= pval])
        #ea_df = read_ea[read_ea["phenotype_id"] =="ENSG00000001461.16"]
        
        # List of gene indexes with at least one variant having significant assoc:
        # Union gives list of genes where at least one cohort has one variant per gene with sig assoc:
        aa_genes = pd.Series(aa_thresh["phenotype_id"].unique(), name="phenotype_id")
        aa_genes = aa_genes.str.replace("\\..*", "") #remove transcript iso
        ea_genes = pd.Series(ea_thresh["phenotype_id"].unique(), name="phenotype_id")
        ea_genes = ea_genes.str.replace("\\..*", "") #remove transcript iso
        ea_aa_union = pd.merge(aa_genes, ea_genes, how="outer")
        shared_cis_genes.append(ea_aa_union.shape[0])

        # Filter gene-variant based on gene idx
        aa_filt = pd.merge(read_aa, ea_aa_union, how="inner") 
        ea_filt = pd.merge(read_ea, ea_aa_union, how="inner")

        """Find shared gene-variants between EA and AA pops on same gene"""
        """ merge afr eur cohorts based on selected gene idx"""
        shared_vars = pd.merge(ea_filt, aa_filt, on=["phenotype_id", "variant_id"], how = "outer", suffixes=["_EA", "_AA"], indicator=True)
        shared_variants = shared_vars[shared_vars["_merge"]=="both"]        
        ea_variants = shared_vars[shared_vars["_merge"]=="left_only"]
        aa_variants = shared_vars[shared_vars["_merge"]=="right_only"]

        pop_betas_ses = shared_variants.loc[:, ["variant_id", "slope_EA", "slope_AA", "slope_se_EA", "slope_se_AA", "maf_EA", "maf_AA", "phenotype_id"]]
        pop_shared_betas_ses.append(pop_betas_ses)
        """Basic summary stats"""
        print('''\tTotal variants across cohorts ::
        European cohorts : {0}
        African cohorts : {1}'''.format(read_ea.shape[0], read_aa.shape[0]))
        print('''\tSignificant gene count with at least 1 significant hit ::
        European cohorts : {0}
        African cohorts : {1}
        EUR-AFR Gene Union : {2}'''.format(ea_genes.shape[0], aa_genes.shape[0], ea_aa_union.shape[0]))
        print('''\t\nGene-Variants with pval {0} thresh ::
        Shared (at least one pop thresh) : {1}
        EUR specific : {2}
        AFR specific : {3}'''.format(pval, shared_variants.shape[0],ea_variants.shape[0],aa_variants.shape[0]))

    betas_ses_concat = pd.concat(pop_shared_betas_ses, ignore_index=True)
    print("\ntotal shared eQTL counts :  {}".format(betas_ses_concat.shape[0]))
    print("\ntotal shared cis_genes counts :  {}".format(sum(shared_cis_genes)))
    #betas_ses_concat.to_csv(join(output_path, tissue + "_coloc_betas_ses.txt"), sep="\t", index=False, header=True)
    return(betas_ses_concat)


# liftover chrom wise:
def liftover(dataframe):
    dataframe = cis_eqtl_coloc_df.copy()
    dataframe.rename(columns={dataframe.columns[0]: "newsnp"}, inplace=True)
    dataframe["newsnp"] = dataframe["newsnp"].str.replace("_b38", "").str.replace("_", ":")
    dataframe["idx"] =  np.arange(len(dataframe))
    dfsubset = dataframe["newsnp"].str.split(":", expand=True)

    #df_test[["chrom", "end", "ref", "alt"]] = pd.DataFrame([ x.split(':') for x in list(df_test['newsnp']) ])
    #df_test[["chrom", "end", "ref", "alt"]] = df_test["newsnp"].str.split(":", return_type='frame')
    dfsubset.columns = ["chrom", "end", "ref", "alt"]
    dataframe["ref"] = dfsubset["ref"]
    dataframe["alt"] = dfsubset["alt"]
    dfsubset["start"] = dfsubset["end"].astype(int) - 1
    dfsubset["idx"] = dataframe["idx"] #.str.cat(dataframe[dataframe.columns[1:].values], sep="_")
    select_cols = ["chrom", "start", "end", "idx"]
    subset_df = dfsubset.loc[:,select_cols]

    pybed_df = pybedtools.BedTool.from_dataframe(subset_df)
    hg19_pybed = pybed_df.liftover(chainfile)
    print("\nLiftover intermediate")
    hg19_df = pd.read_csv(hg19_pybed.fn, sep="\t", header=None)
    hg19_df.columns = ["chrom", "start", "end", "idx"]
    df_merged = pd.merge(hg19_df, dataframe, on="idx", how="inner")
    print("Liftover dataframe merged")

    #df_merged[["chr", "pos", "allele"]] = df_merged["newsnp"].str.split(":", 2, expand=True)
    df_merged["hg19_snp"] = df_merged["chrom"] + ":" +  df_merged["end"].astype(str) + ":" +  \
                            df_merged["ref"] + ":" + df_merged["alt"]
    select_cols = ['hg19_snp', 'newsnp']+ dataframe.columns[1:-3].values.tolist()
    hg19_final_df = df_merged.loc[:,select_cols]
    print("\nLiftover completed")
    return(hg19_final_df)

#hg19_eqtl_coloc_df = liftover(cis_eqtl_coloc_df)


# Common GWAS with hg19 assembly:
def eQTL_gwas_colocformat(eqtl_filt_loci, gwas_summary_stats_file):
    gwas_df = pd.read_csv(gwas_summary_stats_file, sep="\t")
    select_cols = ["variant", "beta", "se"]
    hg19_gwas = gwas_df.loc[:,select_cols]
    hg19_gwas["variant"] =  "chr" + hg19_gwas["variant"].astype(str)

    # merge eQTL GWAS:
    merged_final_df = pd.merge(eqtl_filt_loci, hg19_gwas, left_on=["hg19_snp"], right_on=["variant"], how="inner")
    
    # compute percent overlap:
    unique_eqtl_df = eqtl_filt_loci.drop_duplicates(subset=["hg19_snp"])
    unique_merged_df = merged_final_df.drop_duplicates(subset=["hg19_snp"])
    total_unique_eQTL_count = unique_eqtl_df.shape[0]
    total_eQTL_GWAS_count = unique_merged_df.shape[0]
    percent = round((total_eQTL_GWAS_count/float(total_unique_eQTL_count))*100, 2)
    print("Total shared unique-eQTL counts :  {}".format(total_unique_eQTL_count))
    print("Total eQTL-GWAS counts :  {}".format(total_eQTL_GWAS_count))
    print("eQTL Percent Overlap with GWAS  :  {}%%".format(percent)) 
    return(merged_final_df)

#all_eQTL_gwas_loci = eQTL_gwas_colocformat(hg19_eqtl_coloc_df, gwas_stats)


"""function call"""
if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    start_time = time.time()

    """Input parameters"""
    tissue = "Adipose_Subcutaneous" #tissue = "Heart_Atrial"
    chromrange = 22
    chrom_X = False
    pval = 1e-5
    pval = 1e-4

    # Shared eQTLs and storage in compressed formats: read/write with pickle
    cis_eqtl_coloc_df = pop_shared_variants_colocformat(aa_filepath, ea_filepath, chromrange, pval, tissue, chrom_X=True)    
    cis_eqtl_coloc_df.to_pickle(join(output_path, tissue + "_coloc_betas_ses_1e4.pkl"))
    #cis_eqtl_coloc_df = pd.read_pickle(join(output_path, tissue + "_coloc_betas_ses_1e4.pkl"))

    # Dframe based genomelift of hg38 GTEx tissues to hg19:
    hg19_eqtl_coloc_df = liftover(cis_eqtl_coloc_df)
    hg19_eqtl_coloc_df.to_pickle(join(output_path, tissue + "_coloc_betas_ses_hg19_1e4.pkl"))
    #hg19_eqtl_coloc_df = pd.read_pickle(join(output_path, tissue + "_coloc_betas_ses_hg19_1e4.pkl"))

    # eQTL loci associated to GWAS loci:
    eQTL_gwas_coloci  = eQTL_gwas_colocformat(hg19_eqtl_coloc_df, gwas_stats)
    eQTL_gwas_coloci.drop(["newsnp", "variant"], inplace=True, axis=1)
    eQTL_gwas_coloci.to_csv(join(output_path, tissue + "_coloc_betas_ses_eqtlgwas_1e4.txt"), sep="\t", index=False, header=True )
    eQTL_gwas_coloci.to_pickle(join(output_path, tissue + "_coloc_betas_ses_eqtlgwas_1e4.pkl"))
    #eQTL_gwas_coloci = pd.read_pickle(join(output_path, tissue + "_coloc_betas_ses_eqtlgwas_1e4.pkl"))

    end_time = time.time()
    total_time = end_time - start_time
    print("Time for analysis : {}".format(total_time))
    print("Task completed!")


#group by genes to calculate correlation:
eQTL_gwas_coloci["EA_zscore"] = eQTL_gwas_coloci["slope_EA"].astype(float)/eQTL_gwas_coloci["slope_se_EA"].astype(float)
eQTL_gwas_coloci["AA_zscore"] = eQTL_gwas_coloci["slope_AA"].astype(float)/eQTL_gwas_coloci["slope_se_AA"].astype(float)
grouped_df = eQTL_gwas_coloci.groupby("phenotype_id")
corr = grouped_df[["EA_zscore", "AA_zscore"]].corr()
final_corr = corr.iloc[0::2,-1]
final_corr = final_corr.reset_index().iloc[:, [0,2]]

#####
correlation  = corr.reset_index().iloc[:, [1,2,3]]
correlation.to_csv("ea_aa_cis_effects_correlation.txt", sep="\t", header=True)



