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

ldgenome_path = ""

# Requires `liftOver`from UCSC to be on the path and a `chainfile`
chainfile = "/Users/suryachhetri/tools/chain_files/hg38ToHg19.over.chain"

""" create plot output dir"""
if not os.path.exists(plot_path):
    os.makedirs(plot_path)


"""liftover files"""
inputfile = "/Users/suryachhetri/datasets/encode_screen/GRCh38-ccREs.bed"
filesuffix = ".hg19.bed"

def genomelift(input_bed, file_suffix):
    output_bed = splitext(input_bed)[0] + file_suffix
    unmapped_bed = splitext(input_bed)[0] + file_suffix + ".unmapped.txt"
    cmds = ['liftOver', input_bed, chainfile, output_bed, unmapped_bed]
    os.system(' '.join(cmds))
    return(output_bed)

#liftover_filepath = genomelift(inputfile, filesuffix)


"""calculate variant-gene pairwise slope corr(r) 
    and coeff of determination (rsquare) for each gene"""

def stats_linregress(df): 
    #compute slope,intercept,r_value,p_value,std_err
    result_df = stats.linregress(df["slope_AA"], df["slope_EA"])[2]
    return(result_df)

def stats_zlinregress(df): 
    #compute slope,intercept,r_value,p_value,std_err
    result_df = stats.linregress(df["ZScore_AA"], df["ZScore_EA"])[2]
    return(result_df)

def stats_pearsonr(df):
    result_df = stats.pearsonr(df["ZScore_AA"], df["ZScore_EA"])[0]
    return(result_df)

def stats_spearmanr(df):
    result_df = stats.spearmanr(df["ZScore_AA"], df["ZScore_EA"])[0]
    return(result_df)

"""Distribution differentiation test-statistics"""

def stats_KS_2samp(df):
    result_df = stats.ks_2samp(df["ZScore_AA"], df["ZScore_EA"])[1]
    return(result_df)

def stats_wilcoxon(df):
    result_df = stats.wilcoxon(df["ZScore_AA"], df["ZScore_EA"])[1]
    return(result_df)

def stats_KS_2samp_adj(df):
    result_df = stats.ks_2samp(df["scaled_Z_AA"], df["scaled_Z_EA"])[1]
    return(result_df)

def stats_wilcoxon_adj(df):
    result_df = stats.wilcoxon(df["scaled_Z_AA"], df["scaled_Z_EA"])[1]
    return(result_df)

"""Mean-centered scaling - adjusting marked distribution 
    diff with mean 0 and std dev 1"""

def zscore_scaling(s):
    result_val = (s - s.mean())/s.std()
    return(result_val)

"""Population shared eQTLs compute based on Beta estimates, where 
    Beta estimates are zscored normalized and similarity metrics 
    are calculated with rsquare of Betabases Zscore across acncestry"""

def pop_shared_variants_betabased(aa_filepath, ea_filepath, chromrange, pval, tissue, zscore_thresh, chrom_X=False):
    
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
    eQTL_unfilt_stats = []
    aa_specific_df = []
    ea_specific_df = []
    shared_eQTL_count = []
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

        # Filter gene-variant based on gene idx
        aa_filt = pd.merge(read_aa, ea_aa_union, how="inner") 
        ea_filt = pd.merge(read_ea, ea_aa_union, how="inner")

        """Find shared gene-variants between EA and AA pops on same gene"""
        """ merge afr eur cohorts based on selected gene idx"""
        shared_vars = pd.merge(ea_filt, aa_filt, on=["phenotype_id", "variant_id"], how = "outer", suffixes=["_EA", "_AA"], indicator=True)
        shared_variants = shared_vars[shared_vars["_merge"]=="both"]        
        ea_variants = shared_vars[shared_vars["_merge"]=="left_only"]
        aa_variants = shared_vars[shared_vars["_merge"]=="right_only"]

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

        """calculate correlation"""
        shared_variants["logpval_nominal_EA"] = -(np.log10(shared_variants["pval_nominal_EA"]))
        shared_variants["logpval_nominal_AA"] = -(np.log10(shared_variants["pval_nominal_AA"]))
        shared_variants["ZScore_EA"] = shared_variants["slope_EA"].astype(float)/shared_variants["slope_se_EA"].astype(float)
        shared_variants["ZScore_AA"] = shared_variants["slope_AA"].astype(float)/shared_variants["slope_se_AA"].astype(float)
        shared_variants_df = shared_variants.loc[:,["phenotype_id", "variant_id", "pval_nominal_EA", "pval_nominal_AA", "slope_EA", 
                                                    "slope_AA", "ZScore_EA", "ZScore_AA", "tss_distance_EA", "maf_EA", "maf_AA",
                                                    "logpval_nominal_EA", "logpval_nominal_AA"]]
        shared_variants_df.rename(columns={"tss_distance_EA": "tss_dist"}, inplace=True) #given: [tss_distance_EA==tss_distance_AA]
        shared_variants_df["tss_dist_mb"] = shared_variants_df["tss_dist"].astype(float)/1000000
        shared_variants_df["pval_ratio"] = shared_variants["logpval_nominal_EA"].astype(float)/shared_variants["logpval_nominal_AA"].astype(float)
        shared_variants_df["pval_ratio"] = shared_variants_df["pval_ratio"].astype(float)/1000 #standardize by the factor of 1000

        shared_variants_df.reset_index(drop=True, inplace=True)
        shared_variants_df["scaled_Z_EA"] = shared_variants_df.groupby(["phenotype_id"])["ZScore_EA"].transform(zscore_scaling)
        shared_variants_df["scaled_Z_AA"] = shared_variants_df.groupby(["phenotype_id"])["ZScore_AA"].transform(zscore_scaling)

        """correlation test-statistics"""
        stats_beta = shared_variants_df.groupby(["phenotype_id"]).apply(stats_linregress)
        stats_beta = stats_beta.reset_index(name="Beta_corr")
        stats_pearson = shared_variants_df.groupby(["phenotype_id"]).apply(stats_pearsonr)
        stats_pearson = stats_pearson.reset_index(name="Pearson_corr")
        stats_spearman = shared_variants_df.groupby(["phenotype_id"]).apply(stats_spearmanr)
        stats_spearman = stats_spearman.reset_index(name="Spearman_corr")
        stats_zscore = shared_variants_df.groupby(["phenotype_id"]).apply(stats_zlinregress)
        stats_zscore = stats_zscore.reset_index(name="Zscore_corr")

        """KS test-statistics"""
        stats_kstest = shared_variants_df.groupby(["phenotype_id"]).apply(stats_KS_2samp)
        stats_kstest = stats_kstest.reset_index(name="KS_score")
        stats_kstest_adj = shared_variants_df.groupby(["phenotype_id"]).apply(stats_KS_2samp_adj)
        stats_kstest_adj = stats_kstest_adj.reset_index(name="KS_score_adj")

        """Wilcoxon rank test-statistics"""
        stats_wilcoxon_df = shared_variants_df.groupby(["phenotype_id"]).apply(stats_wilcoxon)
        stats_wilcoxontest = stats_wilcoxon_df.reset_index(name="wilcoxon_score")
        stats_wilcoxontest_adj = shared_variants_df.groupby(["phenotype_id"]).apply(stats_wilcoxon_adj)
        stats_wilcoxontest_adj = stats_wilcoxontest_adj.reset_index(name="wilcoxon_score_adj")

        """for python3 moved to functools"""
        dfs = [stats_beta, stats_pearson, stats_spearman, stats_zscore, stats_kstest, stats_kstest_adj, stats_wilcoxontest, stats_wilcoxontest_adj]
        stats_filt_df = reduce(lambda left, right: pd.merge(left,right,on=["phenotype_id"], how="outer"), dfs)
        stats_filt_df["KS_score_log"] = -(np.log10(stats_filt_df["KS_score"]))
        stats_filt_df["KS_score_adj_log"] = -(np.log10(stats_filt_df["KS_score_adj"]))
        stats_filt_df["wilcoxon_score_log"] = -(np.log10(stats_filt_df["wilcoxon_score"]))
        stats_filt_df["wilcoxon_score_adj_log"] = -(np.log10(stats_filt_df["wilcoxon_score_adj"]))
        stats_filt_df["Beta_rsquare"] = stats_filt_df["Beta_corr"]**2
        stats_filt_df["Zscore_rsquare"] = stats_filt_df["Zscore_corr"]**2

        stats_filt_df["chrom"] = chrom
        cols = stats_filt_df.columns.tolist()
        select_cols = cols[-1:] + cols[:-1] #Flip chrom to first pos
        stats_filt_df = stats_filt_df.loc[:,select_cols]
        #stats_filt_all.append(stats_filt_df)

        """variant, genecount and median stats"""
        corr_median = stats_filt_df["Beta_corr"].median()
        rsquare_median = stats_filt_df["Beta_rsquare"].median()
        Z_rsquare_median = stats_filt_df["Zscore_corr"].median()
        Z_corr_median = stats_filt_df["Zscore_rsquare"].median()
        print("Total filtered genes :  {}".format(stats_filt_df.shape[0]))

        # Merge dataframes after genewise computation of test-stats:
        select_cols = ["phenotype_id", "Zscore_rsquare","wilcoxon_score_adj_log"]
        genewise_select = stats_filt_df.loc[:, select_cols]

        select_cols = ["phenotype_id", "variant_id", "pval_nominal_EA", "pval_nominal_EA", "scaled_Z_EA", "scaled_Z_AA"]
        shared_vars_select = shared_variants_df.loc[:, select_cols]
        shared_vars_select["newsnp"] = shared_vars_select["variant_id"].str.replace("_b38", "").str.replace("_", ":")
        merged_eqtl_stats = pd.merge(shared_vars_select, genewise_select, on="phenotype_id", how="inner")
        merged_eqtl_statsfilt = merged_eqtl_stats.loc[:,["newsnp", "phenotype_id", "Zscore_rsquare", "wilcoxon_score_adj_log"]]
        shared_eQTL_count.append(shared_vars_select.shape[0])

        """prioritizing multigene associated snps based on top-zscore to find unique snps"""
        sorted_df = merged_eqtl_statsfilt.sort_values("Zscore_rsquare", ascending=False)
        grouped = sorted_df.groupby("newsnp")
        zscore_sorted_df = grouped.head(1).reset_index(drop=True)
        # grouped.get_group("chr10:72787275:C:T") #example

        """prioritize multigene associated snps based on wilcoxon_score to find unique snps"""
        # sorted_df1 = merged_eqtl_statsfilt.sort_values("wilcoxon_score_adj_log", ascending=True)
        # grouped1 = sorted_df1.groupby("newsnp")

        # wilcox_sorted_df = grouped1.head(1).reset_index(drop=True)
        # eQTL_merged_df = pd.merge(zscore_sorted_df, wilcox_sorted_df, on=["newsnp"], how="inner")

        """find number of genes associated per snp"""
        genecount_merged_df = grouped.size().reset_index(name="gene_count")
        eQTL_gene_merged_df = pd.merge(zscore_sorted_df, genecount_merged_df, on="newsnp", how="inner")
        
        #select_cols = ['chr', 'end', 'ref', 'alt', 'phenotype_id', 'Zscore_rsquare', 'wilcoxon_score_adj_log', 'gene_count']
        #eQTL_gene_merged_df = eQTL_gene_merged_df.loc[:, select_cols]
        eQTL_unfilt_stats.append(eQTL_gene_merged_df)

    total_shared_eQTL_count = sum(shared_eQTL_count)
    print("\ntotal shared eQTL counts :  {}".format(total_shared_eQTL_count))
    eQTL_unfilt_stats_concat = pd.concat(eQTL_unfilt_stats, ignore_index=True)
    eQTL_unfilt_stats_concat.to_csv(join(output_path, tissue + "_shared_allcis_eQTL_unfilt_stats.txt"), sep="\t", index=False, header=True)
    return(eQTL_unfilt_stats_concat)


# filter shared eQTLs based on rsquare(i.e zscore test-stats):
def eQTL_rsquare_filt(cis_eQTL_raw_df, zscore_filt_thresh):
    unpassed_chrom = []

    # Apply hard-filtering:
    try:
        eQTL_filter_df = cis_eQTL_raw_df[cis_eQTL_raw_df["Zscore_rsquare"] >= zscore_filt_thresh]
        eQTL_filter_df.to_csv(join(output_path, tissue + "_shared_eQTL_filt_stats_{}.txt".format(zscore_filt_thresh)), sep="\t", index=False, header=True)
        #eQTL_filt_stats.append(eQTL_filter_df)

    except ValueError as e:
        unpassed_chrom.append(chrom)
        print(e,"\n Note : {} doesn't have enough genes with the required R2 similarity cutoff".format(chrom))

    unpassed_chrom_df = pd.DataFrame(unpassed_chrom, columns=["chrom"])
    print("Unpassed chroms for set R2 threshold : {}".format(unpassed_chrom))
    return(eQTL_filter_df)


# liftover chrom wise:
def liftover(shared_eqtl_filt_df):
    eqtl_statsfilt_df = shared_eqtl_filt_df["newsnp"].str.split(":", expand=True)
    eqtl_statsfilt_df.columns = ["chrom", "end", "ref", "alt"]
    eqtl_statsfilt_df["start"] = eqtl_statsfilt_df["end"].astype(int) - 1
    eqtl_statsfilt_df["newsnp"] = shared_eqtl_filt_df["newsnp"]
    select_cols = ["chrom", "start", "end", "newsnp"]
    select_eqtl_df = eqtl_statsfilt_df.loc[:,select_cols]
    pybed_df = pybedtools.BedTool.from_dataframe(select_eqtl_df)
    hg19_eqtl_pybed = pybed_df.liftover(chainfile)
    hg19_eqtl_df = pd.read_csv(hg19_eqtl_pybed.fn, sep="\t", header=None)
    hg19_eqtl_df.columns = ["chrom", "start", "end", "newsnp"]
    hg19_eqtl_df[["chr", "pos", "allele"]] = hg19_eqtl_df["newsnp"].str.split(":", 2, expand=True)
    hg19_eqtl_df["hg19_snp"] = hg19_eqtl_df["chrom"] + ":" +  hg19_eqtl_df["end"].astype(str) + ":" +  hg19_eqtl_df["allele"]
    select_cols = ['hg19_snp', 'newsnp']
    hg19_eqtl_final_df = hg19_eqtl_df.loc[:,select_cols]
    return(hg19_eqtl_final_df)


# Common GWAS with hg19 assembly:
def eQTL_gwas_loci(eqtl_filt_loci, gwas_summary_stats_file):
    gwas_df = pd.read_csv(gwas_summary_stats_file, sep="\t")
    select_cols = ["variant", "beta", "pval"]
    hg19_gwas = gwas_df.loc[:,select_cols]
    hg19_gwas["variant"] =  "chr" + hg19_gwas["variant"].astype(str)

    # merge eQTL GWAS:
    merged_final_df = pd.merge(eqtl_filt_loci, hg19_gwas, left_on=["hg19_snp"], right_on=["variant"], how="inner")
    
    # compute percent overlap:
    unique_eqtl_df = eqtl_filt_loci.drop_duplicates(subset=["hg19_snp"])
    unique_merged_df = merged_final_df.drop_duplicates(subset=["hg19_snp"])
    total_unique_eQTL_count = unique_eqtl_df.shape[0]
    total_eQTL_GWAS_count = merged_final_df.shape[0]
    percent = round((total_eQTL_GWAS_count/float(total_unique_eQTL_count))*100, 2)
    print("Total shared unique-eQTL counts :  {}".format(total_unique_eQTL_count))
    print("Total eQTL-GWAS counts :  {}".format(total_eQTL_GWAS_count))
    print("eQTL Percent Overlap with GWAS  :  {}%%".format(percent))
    return(unique_merged_df)


# PRS formatting:
def prs_format(eqtl_gwas_loci_df, zscore_thresh):
    eqtl_gwas_loci_df["hg19_snp"] = eqtl_gwas_loci_df["hg19_snp"].str.replace("chr", "")
    eqtl_gwas_loci_df[["chr", "pos", "ref", "alt"]] = eqtl_gwas_loci_df["hg19_snp"].str.split(":", expand=True)
    prs_format = eqtl_gwas_loci_df.loc[:, ["hg19_snp", "ref", "alt", "beta", "pval"]]
    prs_format.rename(columns={"hg19_snp" : "ID", "ref" : "REF", "alt" : "ALT"}, inplace=True)
    prs_format.to_csv(join(output_path, tissue + "_prs_format_{}.txt".format(zscore_thresh)), sep="\t", index=False, header=True)
    return(prs_format)


"""function call"""
if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    start_time = time.time()

    """Input parameters"""
    tissue = "Adipose_Subcutaneous" #tissue = "Heart_Atrial"
    chromrange = 22
    chrom_X = False
    pval = 10e-5
    zscorethresh = 0.1

   # Shared eQTLs and storage in compressed formats: read/write with pickle
    cis_eqtl_unfilt_df = pop_shared_variants_betabased(aa_filepath, ea_filepath, chromrange, pval, tissue, zscorethresh, chrom_X=True)    
    cis_eqtl_unfilt_df.to_pickle(join(output_path, tissue + "_shared_allcis_eQTL_unfilt_stats.pkl"))
    cis_eqtl_unfilt_df = pd.read_pickle(join(output_path, tissue + "_shared_allcis_eQTL_unfilt_stats.pkl"))

    # Filter eQTLs based on rsquare threshold based on zscBeta:
    eqtl_filt_df = eQTL_rsquare_filt(cis_eqtl_unfilt_df, 0.1)
    eqtl_unfilt_df = cis_eqtl_unfilt_df.copy()

    # Dframe based genomelift of hg38 GTEx tissues to hg19:
    hg19_eqtl_filt_df = liftover(eqtl_filt_df)
    hg19_eqtl_all_df = liftover(eqtl_unfilt_df)

    # eQTL loci associated to GWAS loci:
    filt_eQTL_gwas_loci  = eQTL_gwas_loci(hg19_eqtl_filt_df, gwas_stats)
    all_eQTL_gwas_loci  = eQTL_gwas_loci(hg19_eqtl_all_df, gwas_stats)

    # PRS format adaptation(threshold input for unique file-naming only):
    filename_thresh = 0.1 #for naming purpose only
    filt_eqtl_prs_format = prs_format(filt_eQTL_gwas_loci, filename_thresh)
    unfilt_eqtl_prs_format = prs_format(all_eQTL_gwas_loci, "unfilt")


    # Baseline GWAS hg19-variants:
    gwas_df = pd.read_csv(gwas_stats, sep="\t", usecols=["variant", "beta", "pval"])
    gwas_df[["CHR", "POS", "REF", "ALT"]] = gwas_df["variant"].str.split(":", expand=True)
    gwas_df.rename(columns={"variant" : "ID"}, inplace=True)
    default_format = gwas_df.loc[:, ["ID", "REF", "ALT", "beta", "pval"]]
    default_format.to_csv(join(output_path, tissue + "_prs_format_gwasbaseline.txt"), sep="\t", index=False, header=True)


    #################################
    #eQTL based weighted beta:
    
    # eQTL weighted GWAS loci - based on zscore_rsuare 
    all_eQTL_weight_df = pd.merge(eqtl_unfilt_df, hg19_eqtl_all_df, on="newsnp", how="inner")
    all_eQTL_select = all_eQTL_weight_df.loc[:, ["hg19_snp", "Zscore_rsquare", "wilcoxon_score_adj_log"]]
    all_eQTL_select["hg19_snp"] = all_eQTL_select["hg19_snp"].str.replace("chr", "") #match chr ID
    all_eQTL_weighted_gwas = pd.merge(all_eQTL_gwas_loci, all_eQTL_select, on="hg19_snp", how="inner")
    uniq_all_eQTL_merged_gwas = all_eQTL_weighted_gwas.drop_duplicates(subset=["hg19_snp"])
    
    #Check if all initial eQTL-GWAS loci have being weighted:
    all_eQTL_gwas_loci.shape[0] == uniq_all_eQTL_merged_gwas.shape[0]
    uniq_all_eQTL_merged_gwas["pears_weighted_beta"] =  uniq_all_eQTL_merged_gwas["beta"] * uniq_all_eQTL_merged_gwas["Zscore_rsquare"]
    uniq_all_eQTL_merged_gwas["wilcox_weighted_beta"] =  uniq_all_eQTL_merged_gwas["beta"]/uniq_all_eQTL_merged_gwas["wilcoxon_score_adj_log"]


    #preping prs format for pears weighted beta:
    debug_df = uniq_all_eQTL_merged_gwas.loc[:, ["hg19_snp", "beta", "Zscore_rsquare", "wilcox_weighted_beta"]]
    select_cols = ["hg19_snp", "newsnp", "variant", "pears_weighted_beta", "pval", "chr", "pos", "ref", "alt"]
    all_eQTL_gwas_wt_loci = uniq_all_eQTL_merged_gwas.loc[:, select_cols]
    all_eQTL_gwas_wt_loci.rename(columns={"pears_weighted_beta" : "beta"}, inplace=True)

    
    #preping prs format for wilcox weighted beta:
    debug_df = uniq_all_eQTL_merged_gwas.loc[:, ["hg19_snp", "beta", "Zscore_rsquare", "wilcox_weighted_beta"]]
    select_cols = ["hg19_snp", "newsnp", "variant", "wilcox_weighted_beta", "pval", "chr", "pos", "ref", "alt"]
    all_eQTL_gwas_wt_loci1 = uniq_all_eQTL_merged_gwas.loc[:, select_cols]
    all_eQTL_gwas_wt_loci1.rename(columns={"wilcox_weighted_beta" : "beta"}, inplace=True)


    # PRS format adaptation(threshold input for unique file-naming only):
    filename_thresh = "weighted" #for naming purpose only
    all_eqtl_prs_format = prs_format(all_eQTL_gwas_wt_loci, filename_thresh)

    filename_thresh = "wilcoxweighted" #for naming purpose only
    all_eqtl_prs_format = prs_format(all_eQTL_gwas_wt_loci1, filename_thresh)

    total_time = time.time() - start_time
    print("Time for analysis : {}".format(total_time))
    print("Task completed!")



# Weight effects on Beta:
# In [545]: all_eQTL_gwas_wt_loci1.head()
# Out[545]:
#           hg19_snp              newsnp             variant      beta          pval chr        pos ref alt
# 0  1:151010521:G:C  chr1:151038045:G:C  chr1:151010521:G:C -0.007960  8.502000e-11   1  151010521   G   C
# 1  1:150510633:C:A  chr1:150538157:C:A  chr1:150510633:C:A -0.000003  9.975320e-01   1  150510633   C   A
# 2  1:150514747:C:T  chr1:150542271:C:T  chr1:150514747:C:T  0.005301  2.650860e-03   1  150514747   C   T
# 3  1:150514741:T:G  chr1:150542265:T:G  chr1:150514741:T:G  0.001702  1.105090e-01   1  150514741   T   G
# 4  1:150514568:C:T  chr1:150542092:C:T  chr1:150514568:C:T -0.001904  2.314570e-01   1  150514568   C   T

# In [546]: all_eQTL_gwas_wt_loci.head()
# Out[546]:
#           hg19_snp              newsnp             variant      beta          pval chr        pos ref alt
# 0  1:151010521:G:C  chr1:151038045:G:C  chr1:151010521:G:C -0.012087  8.502000e-11   1  151010521   G   C
# 1  1:150510633:C:A  chr1:150538157:C:A  chr1:150510633:C:A -0.000005  9.975320e-01   1  150510633   C   A
# 2  1:150514747:C:T  chr1:150542271:C:T  chr1:150514747:C:T  0.008048  2.650860e-03   1  150514747   C   T
# 3  1:150514741:T:G  chr1:150542265:T:G  chr1:150514741:T:G  0.002584  1.105090e-01   1  150514741   T   G
# 4  1:150514568:C:T  chr1:150542092:C:T  chr1:150514568:C:T -0.002892  2.314570e-01   1  150514568   C   T

# In [547]: all_eQTL_gwas_loci.head()
# Out[547]:
#           hg19_snp              newsnp             variant      beta          pval chr        pos ref alt
# 0  1:151010521:G:C  chr1:151038045:G:C  chr1:151010521:G:C -0.017937  8.502000e-11   1  151010521   G   C
# 1  1:150510633:C:A  chr1:150538157:C:A  chr1:150510633:C:A -0.000007  9.975320e-01   1  150510633   C   A
# 2  1:150514747:C:T  chr1:150542271:C:T  chr1:150514747:C:T  0.011944  2.650860e-03   1  150514747   C   T
# 3  1:150514741:T:G  chr1:150542265:T:G  chr1:150514741:T:G  0.003834  1.105090e-01   1  150514741   T   G
# 4  1:150514568:C:T  chr1:150542092:C:T  chr1:150514568:C:T -0.004291  2.314570e-01   1  150514568   C   T


> summary(data.frame(mat_data)$PP.H4)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.00000 0.09392 0.18522 0.31605 0.42250 1.00000
> mat_data %>% data.frame %>% filter(PP.H4>0.25) %>% dim
[1] 4601    5
> mat_data %>% data.frame %>% filter(PP.H4>0.1) %>% dim
[1] 8719    5
