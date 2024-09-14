
#!/usr/bin/env python

import numpy as np, pandas as pd
import pybedtools
import os, re
from glob import glob
from os.path import basename, join, splitext
from collections import defaultdict
import time, warnings
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy import stats
import seaborn as sns

""" Input parameters """
aa_filepath = "/Users/suryachhetri/datasets/prs_project/gtex/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_AFR_eQTL_all_associations"
ea_filepath = "/Users/suryachhetri/datasets/prs_project/gtex/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations"
gwas_stats = "/Users/suryachhetri/datasets/prs_project/gwas_stats/bmi/21001_irnt.gwas.imputed_v3.both_sexes.tsv"
#lifted_gwas_stats = "/Users/suryachhetri/datasets/prs_project/gwas_stats/bmi/lifted/updated_output.bed"
tissue = "Adipose_Subcutaneous" #tissue = "Heart_Atrial"
chromrange = 22
chrom_X = False
pval = 10e-5

output_path = "/Users/suryachhetri/datasets/prs_project/final_hg19"
plot_path = join(output_path, "plots")
plot_path = join(output_path, "plots_chrX")
plot_path = join(output_path, "redo_chr1")

ldgenome_path = ""

#liftover="/Users/suryachhetri/tools/liftOver"
#Needs `liftOver` from UCSC to be on the path and a `chainfile`
#downloaded from UCSC.
input_bed="/Users/suryachhetri/datasets/encode_screen/GRCh38-ccREs.bed"
chainfile="/Users/suryachhetri/tools/chain_files/hg38ToHg19.over.chain"
output_bed=input_bed + ".hg19.test.bed"
unmapped_bed=input_bed + ".hg19.test.unmapped.bed"
cmds = ['liftOver', input_bed, chainfile, output_bed, unmapped_bed]
#os.system(' '.join(cmds))

""" create plot output dir"""
if not os.path.exists(plot_path):
    os.makedirs(plot_path)


# # Update GWAS with new hg38 assembly:
# gwas_df = pd.read_csv(gwas_stats, sep="\t")
# select_cols = ["variant", "beta", "pval"]
# gwas_select = gwas_df.loc[:,select_cols]

# lifted_gwas = pd.read_csv(lifted_gwas_stats, sep="\t", header=None)
# lifted_gwas.columns = ["chrom", "start", "end", "oldsnp", "newsnp"]
# lifted_gwas["newsnp"] =  "chr" + lifted_gwas["newsnp"].astype(str)
# lifted_gwas.reset_index(drop=True, inplace=True)

# merged_gwas = pd.merge(gwas_select, lifted_gwas, left_on=["variant"], right_on=["oldsnp"], how="inner")
# select_cols = ["newsnp", "beta", "pval"]
# updated_gwas = merged_gwas.loc[:,select_cols]

# Common GWAS with hg19 assembly:
gwas_df = pd.read_csv(gwas_stats, sep="\t")
select_cols = ["variant", "beta", "pval"]
hg19_gwas = gwas_df.loc[:,select_cols]
hg19_gwas["variant"] =  "chr" + hg19_gwas["variant"].astype(str)

# lifted_gwas = pd.read_csv(lifted_gwas_stats, sep="\t", header=None)
# lifted_gwas.columns = ["chrom", "start", "end", "oldsnp", "newsnp"]
# lifted_gwas["newsnp"] =  "chr" + lifted_gwas["newsnp"].astype(str)
# lifted_gwas.reset_index(drop=True, inplace=True)

# merged_gwas = pd.merge(hg19_gwas, lifted_gwas, left_on=["variant"], right_on=["oldsnp"], how="inner")
# select_cols = ["newsnp", "beta", "pval"]
# updated_gwas = merged_gwas.loc[:,select_cols]

def pop_shared_variants_betabased(aa_filepath, ea_filepath, chromrange, pval, tissue, chrom_X=False):
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

    """calculate variant-gene pairwise slope corr(r) 
        and coeff of determination (rsquare) for each gene"""

    def stats_linregress(df): 
        return stats.linregress(df["slope_AA"], df["slope_EA"])[2] #slope,intercept,r_value,p_value,std_err

    def stats_zlinregress(df): 
        return stats.linregress(df["ZScore_AA"], df["ZScore_EA"])[2] #slope,intercept,r_value,p_value,std_err

    def stats_pearsonr(df):
        return stats.pearsonr(df["ZScore_AA"], df["ZScore_EA"])[0]

    def stats_spearmanr(df):
        return stats.spearmanr(df["ZScore_AA"], df["ZScore_EA"])[0]

    """Distribution differentiation test-statistics"""
    def stats_KS_2samp(df):
        return stats.ks_2samp(df["ZScore_AA"], df["ZScore_EA"])[1]

    def stats_wilcoxon(df):
        return stats.wilcoxon(df["ZScore_AA"], df["ZScore_EA"])[1]

    def stats_KS_2samp_adj(df):
        return stats.ks_2samp(df["scaled_Z_AA"], df["scaled_Z_EA"])[1]

    def stats_wilcoxon_adj(df):
        return stats.wilcoxon(df["scaled_Z_AA"], df["scaled_Z_EA"])[1]

    """Mean-centered scaling - adjusting marked distribution diff with mean 0 per std dev 1"""
    def zscore_scaling(s):
        return (s - s.mean())/s.std()

    #shared_variant_set2 = []
    stats_filt_all = []
    prs_filt_stat = []
    aa_specific_df = []
    ea_specific_df = []
    shared_eQTL_count = []
    shared_eQTL_GWAS_count = []
    unpassed_chrom = []
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
        ea_genes = pd.Series(ea_thresh["phenotype_id"].unique(), name="phenotype_id")
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
        from functools import reduce 
        dfs = [stats_beta, stats_pearson, stats_spearman, stats_zscore, stats_kstest, stats_kstest_adj, stats_wilcoxontest, stats_wilcoxontest_adj]
        stats_filt_df = reduce(lambda left, right: pd.merge(left,right,on=["phenotype_id"], how="outer"), dfs)
        stats_filt_df["KS_score_log"] = -(np.log10(stats_filt_df["KS_score"]))
        stats_filt_df["KS_score_adj_log"] = -(np.log10(stats_filt_df["KS_score_adj"]))
        stats_filt_df["wilcoxon_score_log"] = -(np.log10(stats_filt_df["wilcoxon_score"]))
        stats_filt_df["wilcoxon_score_adj_log"] = -(np.log10(stats_filt_df["wilcoxon_score_adj"]))
        stats_filt_df["Beta_rsquare"] = stats_filt_df["Beta_corr"]**2
        stats_filt_df["Zscore_rsquare"] = stats_filt_df["Zscore_corr"]**2

        # select_cols = ["ZScore_AA",  "ZScore_EA",  "scaled_Z_AA",  "scaled_Z_EA"]
        # stats_filt_df[stats_filt_df["phenotype_id"] == "ENSG00000280670.2"].loc[:, 
        #    ["KS_score_adj", "KS_score_adj_log", "wilcoxon_score_adj", "wilcoxon_score_adj_log"]]

        stats_filt_df["chrom"] = chrom
        cols = stats_filt_df.columns.tolist()
        select_cols = cols[-1:] + cols[:-1] #Flip chrom to first pos
        stats_filt_df = stats_filt_df.loc[:,select_cols]
        stats_filt_all.append(stats_filt_df)

        """variant, genecount and median stats"""
        #variantcount.append(shared_variants_df.shape[0])
        #genecount.append(stats_filt_df.shape[0])
        corr_median = stats_filt_df["Beta_corr"].median()
        rsquare_median = stats_filt_df["Beta_rsquare"].median()
        Z_rsquare_median = stats_filt_df["Zscore_corr"].median()
        Z_corr_median = stats_filt_df["Zscore_rsquare"].median()
        print("Total filtered genes :  {}".format(stats_filt_df.shape[0]))

        # Merge dataframes:
        select_cols = ["phenotype_id", "Zscore_rsquare","wilcoxon_score_adj_log"]
        genewise_select = stats_filt_df.loc[:, select_cols]

        select_cols = ["phenotype_id", "variant_id", "pval_nominal_EA", "pval_nominal_EA", "scaled_Z_EA", "scaled_Z_AA"]
        shared_vars_select = shared_variants_df.loc[:, select_cols]
        shared_vars_select["newsnp"] = shared_vars_select["variant_id"].str.replace("_b38", "").str.replace("_", ":")
        merged_eqtl_stats = pd.merge(shared_vars_select, genewise_select, on="phenotype_id", how="inner")
        merged_eqtl_statsfilt = merged_eqtl_stats.loc[:,["newsnp", "phenotype_id", "Zscore_rsquare", "wilcoxon_score_adj_log"]]
        shared_eQTL_count.append(shared_vars_select.shape[0])

        # liftover chrom wise:
        merged_eqtl_df = merged_eqtl_statsfilt["newsnp"].str.split(":", expand=True)
        merged_eqtl_df.columns = ["chrom", "end", "ref", "alt"]
        merged_eqtl_df["start"] = merged_eqtl_df["end"].astype(int) - 1
        merged_eqtl_df["newsnp"] = merged_eqtl_statsfilt["newsnp"]
        select_cols = ["chrom", "start", "end", "newsnp"]
        merged_eqtl_df = merged_eqtl_df.loc[:,select_cols]
        pybed_df = pybedtools.BedTool.from_dataframe(merged_eqtl_df)
        hg19_eqtl_pybed = pybed_df.liftover(chainfile)
        hg19_eqtl_df = pd.read_csv(hg19_eqtl_pybed.fn, sep="\t", header=None)
        hg19_eqtl_df.columns = ["chrom", "start", "end", "newsnp"]
        hg19_eqtl_df[["chr", "pos", "allele"]] = hg19_eqtl_df["newsnp"].str.split(":", 2, expand=True)
        hg19_eqtl_df["snpidx"] = hg19_eqtl_df["chrom"] + ":" +  hg19_eqtl_df["end"].astype(str) + ":" +  hg19_eqtl_df["allele"]
        select_cols = ['snpidx', 'newsnp']
        hg19_eqtl_df = hg19_eqtl_df.loc[:,select_cols]

        # merge info for eQTl stats
        hg19_eqtl_statsfilt = pd.merge(merged_eqtl_statsfilt, hg19_eqtl_df, on="newsnp", how="inner")



        """merge eqtl and gwas"""
        hg19_gwas.rename(columns={"variant": "newsnp"}, inplace=True)
        hg19_gwas_eqtl = pd.merge(merged_eqtl_statsfilt, hg19_gwas, on="newsnp", how="inner")
        #gwas_eqtl.loc[:, ["phenotype_id", "newsnp", "Zscore_rsquare", "wilcoxon_score_adj", "wilcoxon_score_adj_log"]]
        shared_eQTL_GWAS_count.append(gwas_eqtl.shape[0])

        eqtl_percent = round((gwas_eqtl.shape[0]/float(shared_vars_select.shape[0]))*100, 2)
        print("eQTL counts :  {}".format(shared_eQTL_count))
        print("eQTL-GWAS counts :  {}".format(shared_eQTL_GWAS_count))
        print("eQTL perc overlap with GWAS  :  {}%%".format(eqtl_percent))
     

        df1 = gwas_eqtl.sort_values("Zscore_rsquare", ascending=False)
        grouped = df1.groupby("newsnp")
        output1 = grouped.head(1).reset_index(drop=True)

        df2 = gwas_eqtl.sort_values("wilcoxon_score_adj_log", ascending=True)
        grouped = df2.groupby("newsnp")
        output_1 = grouped.head(1).reset_index(drop=True)

        output_2 = grouped.size().reset_index(name="gene_count")
        merged_df = pd.merge(output1, output_1, on=["newsnp", "beta", "pval"], how="inner")
        final_df = pd.merge(merged_df, output_2, on="newsnp", how="inner")

        # Apply hard-filtering:
        try:
            filter_df = final_df[final_df["Zscore_rsquare_x"] >= 0.1]
            #filter_df = final_df.copy()
            prs_df = filter_df.loc[:, ["newsnp", "beta", "pval"]]
            prs_df["newsnp"] = prs_df["newsnp"].str.replace("chr", "")
            prs_df[["chr", "pos", "ref", "alt"]] = prs_df["newsnp"].str.split(":", expand=True)
            prs_format = prs_df.loc[:, ["newsnp", "ref", "alt", "beta", "pval"]]
            prs_format.rename(columns={"newsnp" : "ID", "ref" : "REF", "alt" : "ALT"}, inplace=True)
            prs_filt_stat.append(prs_format)

        except ValueError as e:
            unpassed_chrom.append(chrom)
            print(e,"\n Note : {} doesn't have enough genes with the required R2 similarity cutoff".format(chrom))



    """summing stats and df for all chroms"""
    #total_genecount = sum(genecount)
    #total_variantcount = sum(variantcount)

    total_eQTL_count = sum(shared_eQTL_count)
    total_eQTL_GWAS_count = sum(shared_eQTL_GWAS_count)
    percent = round((total_eQTL_GWAS_count/float(total_eQTL_count))*100, 2)
    print("Total eQTL counts :  {}".format(total_eQTL_count))
    print("Total eQTL-GWAS counts :  {}".format(total_eQTL_GWAS_count))
    print("eQTL Percent Overlap with GWAS  :  {}%%".format(percent))
    print("Unpassed chroms for set R2 threshold : {}".format(unpassed_chrom))
    prs_concat_stat = pd.concat(prs_filt_stat, ignore_index=True)
    prs_concat_stat.to_csv(join(output_path, tissue + "_prs_ashton_format_0.8.txt"), sep="\t", index=False, header=True)
    stats_filt_comp = pd.concat(stats_filt_all, ignore_index=True)
    stats_filt_comp.to_csv(join(output_path, tissue+"_shared_gene_corr_coeff_stats_0.8.txt"), sep="\t", index=False, header=True)
    unpassed_chrom_df = pd.DataFrame(unpassed_chrom, columns=["chrom"])
    stats_filt_comp.to_csv(join(output_path, tissue+"_R2_unpassed_chroms_0.8.txt"), sep="\t", index=False, header=True)

    return(prs_filt_stat, stats_filt_comp)

"""function call"""
if __name__ == "__main__":
    #import warnings
    warnings.filterwarnings("ignore")
    start_time = time.time()
#    chromrange = 1
    result_df = pop_shared_variants_betabased(aa_filepath, ea_filepath, chromrange, pval, tissue, chrom_X=True)
    total_time = time.time() - start_time
    print("Time for analysis : {}".format(total_time))
    print("Task completed!")

    """ boxplot for genomewide beta corr and coeff variation """
    df_result = result_df[2]
    df_new = df_result.loc[:, ["chrom", "beta_correlation", "rsquare"]]
    df_melt = pd.melt(df_new, "chrom", var_name="measurement")

    """ using seaborn backend and placed in subplot"""
    sns.set(style="ticks")
    fig, ax = plt.subplots()

    # Draw a nested boxplot:
    # fig.set_size_inches(8, 6)
    axs = sns.boxplot(x="chrom", y="value",
                hue="measurement", palette=["r", "b"],
                data=df_melt, \
                ax=ax, \
                #"""remove outlier"""
                flierprops = dict(markerfacecolor = '0.3', markersize = 1))

    sns.despine(offset=5, trim=True)

    # set axes and labels for sns subplot:
    axs.set_xlabel("Chromosome", fontsize=10)
    axs.set_ylabel("Correlation | Coeff of Variation", fontsize=10)
    axs.set_title("Genome Wide Distribution (r|rsquare)", weight='bold', fontsize=12)
    axs.tick_params(axis='x', labelrotation=45, labelcolor='black', labelsize=9)
    axs.tick_params(axis='y', labelrotation=0, labelcolor='black', labelsize=9)

    # Put the legend out of the figure
    ax.legend(loc=4, prop={'size': 10}, edgecolor="black")
    ax.axhline(0.8, color='tab:blue', linestyle='dashed', linewidth=2)

    # Save figure:
    fig.tight_layout()
    fig.savefig(join(plot_path, '{}_genomewide_betaStats_coeff_boxplot.pdf'.format(tissue)))
    plt.close()


aa_filt0 = pd.read_csv("Adipose_Subcutaneous_prs_format_AA_0.1.full_prs.scores_AA_BMI_age_r2counts.tsv", sep="\t")
aa_filt1 = pd.Series([546, 930, 2661, 5378, 12922, 37672, 86918, 115508], name="snp_count")
aa_filt = pd.concat([aa_filt0, pd.DataFrame(aa_filt1)], axis=1)

eur_filt0 = pd.read_csv("Adipose_Subcutaneous_prs_format_EUR_0.1.full_prs.scores_EUR_BMI_age_r2counts.tsv", sep="\t")
eur_filt1 = pd.Series([548, 938, 2684, 5427, 13041, 37844, 87018, 115508], name="snp_count")
eur_filt = pd.concat([eur_filt0, pd.DataFrame(eur_filt1)], axis=1)

aa_unfilt0 = pd.read_csv("Adipose_Subcutaneous_prs_format_unfilt_AA_unfilt.full_prs.scores_AA_BMI_age_r2counts.tsv", sep="\t")
aa_unfilt1 = pd.Series([823, 1475, 4503, 9455, 24003, 73016, 172639, 231011], name="snp_count")
aa_unfilt = pd.concat([aa_unfilt0, pd.DataFrame(aa_unfilt1)], axis=1)


eur_unfilt0 = pd.read_csv("Adipose_Subcutaneous_prs_format_unfilt_EUR_unfilt.full_prs.scores_EUR_BMI_age_r2counts.tsv", sep="\t")
eur_unfilt1 = pd.Series([826, 1484, 4528, 9504, 24136, 73223, 172697, 230867], name="snp_count")
eur_unfilt = pd.concat([eur_unfilt0, pd.DataFrame(eur_unfilt1)], axis=1)

aa_weighted0 = pd.read_csv("Adipose_Subcutaneous_prs_format_weighted_AA_weighted.full_prs.scores_AA_BMI_age_r2counts.tsv", sep="\t")
aa_weighted1 = pd.Series([823, 1475, 4503, 9455, 24003, 73016, 172639, 231011], name="snp_count")
aa_weighted = pd.concat([aa_weighted0, pd.DataFrame(aa_weighted1)], axis=1)


eur_weighted0 = pd.read_csv("Adipose_Subcutaneous_prs_format_weighted_EUR_weighted.full_prs.scores_EUR_BMI_age_r2counts.tsv", sep="\t")
eur_weighted1 = pd.Series([826, 1484, 4528, 9504, 24136, 73223, 172697, 230867], name="snp_count")
eur_weighted = pd.concat([eur_weighted0, pd.DataFrame(eur_weighted1)], axis=1)

aa_gwasbaseline0 = pd.read_csv("Adipose_Subcutaneous_prs_format_gwasbaseline_AA_baseline.full_prs.scores_AA_BMI_age_r2counts.tsv", sep="\t")
aa_gwasbaseline1 = pd.Series([1022, 1838, 5866, 12762, 33863, 108650, 266294, 363893], name="snp_count")
aa_gwasbaseline = pd.concat([aa_gwasbaseline0, pd.DataFrame(aa_gwasbaseline1)], axis=1)

eur_gwasbaseline0 = pd.read_csv("Adipose_Subcutaneous_prs_format_gwasbaseline_EUR_baseline.full_prs.scores_EUR_BMI_age_r2counts.tsv", sep="\t")
eur_gwasbaseline1 = pd.Series([1067, 1944, 6553, 15027, 43379, 153570, 405562, 568985], name="snp_count")
eur_gwasbaseline = pd.concat([eur_gwasbaseline0, pd.DataFrame(eur_gwasbaseline1)], axis=1)


frames = [aa_filt, eur_filt, aa_unfilt, eur_unfilt, aa_weighted, eur_weighted, aa_gwasbaseline, eur_gwasbaseline]
result = pd.concat(frames, keys=['aa_filt', 'eur_filt', 'aa_unfilt', 'eur_unfilt', 'aa_weighted', 'eur_weighted', 'aa_gwasbaseline', 'eur_gwasbaseline'])

final_df = result.reset_index().drop(columns="level_1", axis=0)
final_df = final_df.rename(columns={"level_0": "Category"})
final_select_df = final_df.loc[~final_df["Category"].str.contains("weighted")]


"""BMI"""
#bmi_df1 = bmi_df[~(bmi_df["Cohorts"] == "ALL")]

g = sns.barplot(x = 'pval_names', y = 'r2', hue = 'Category', data = final_df,
            palette = 'hls',
            #col="Cohorts",
            #order = ['aa_filt', 'eur_filt', 'aa_unfilt', 'eur_unfilt', 'aa_weighted', 'eur_weighted', 'aa_gwasbaseline', 'eur_gwasbaseline'],  
            capsize = 0.05,             
            saturation = 8,             
            errcolor = 'gray', errwidth = 2,  
            ci = 'sd'   
            )

box = g.get_position()
g.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # resize position

# Put a legend to the right side
g.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)

plot_path="/Users/suryachhetri/datasets/prs_project/final_hg19/outputs"
plt.savefig(join(plot_path, 'first_plot.pdf'))
plt.close()



"""BMI"""
g = sns.barplot(x = 'pval_names', y = 'R2perSNP', hue = 'Category', data = final_df,
            palette = 'hls',
            #col="Cohorts",
            #order = ['aa_filt', 'eur_filt', 'aa_unfilt', 'eur_unfilt', 'aa_weighted', 'eur_weighted', 'aa_gwasbaseline', 'eur_gwasbaseline'],  
            capsize = 0.05,             
            saturation = 8,             
            errcolor = 'gray', errwidth = 2,  
            ci = 'sd'   
            )

box = g.get_position()
g.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # resize position

# Put a legend to the right side
g.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)

plot_path="/Users/suryachhetri/datasets/prs_project/final_hg19/outputs"
plt.savefig(join(plot_path, 'second_plot.pdf'))
plt.close()


#########################

final_select_df = final_df.loc[~final_df["Category"].str.contains("weighted")]
final_select_df["Group"] = ""
final_select_df.loc[final_select_df["Category"].str.contains("_filt"), "Group"] = "eQTL-informed-genes"
final_select_df.loc[final_select_df["Category"].str.contains("_unfilt"), "Group"] = "cis-genes"
final_select_df.loc[final_select_df["Category"].str.contains("gwasbaseline"), "Group"] = "gwas-baseline"
final_select_df["R2perSNP"] = final_select_df["r2"]/final_df["snp_count"].astype(float)

#final_df1 = final_df[final_df["Category"] == "eQTL Informed"]
sns.catplot(x = 'pval_names', y = 'r2', hue = 'Group', data = final_select_df,
            palette = 'hls',
            col="cohort",
            kind="bar",
            height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
#plt.legend()
plt.savefig(join(plot_path, 'CombinedplotCohorts.pdf'))
plt.close()



sns.catplot(x = 'pval_names', y = 'R2perSNP', hue = 'Group', data = final_select_df,
            palette = 'hls',
            col="cohort",
            kind="bar",
            height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
#plt.legend()
plt.savefig(join(plot_path, 'CombinedplotCohorts_R2perSNP.pdf'))
plt.close()


final_select_df1 = final_select_df[~(final_select_df["cohort"] == "AA")]
sns.catplot(x = 'pval_names', y = 'snp_count', hue = 'Group', data = final_select_df1,
            palette = 'hls',
            col="cohort",
            kind="bar",
            height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
#plt.legend()
plt.savefig(join(plot_path, 'CombinedplotCohorts_EURsnpcount.pdf'))
plt.close()









#!scp prs_test_sample2.txt marcc:/work-zfs/abattle4/surya/datasets/prs_project/test_samp
#prs_test_sample2.txt
filt_path = "/Users/suryachhetri/datasets/prs_project/final/outputs/rsq_0.1/groupwise/plots"
bmi_df = pd.read_csv(join(filt_path, "combined_Adipose_Subcutaneous_prs_0.1_BMI.txt"), sep="\t", names=["P-Value", "R2", "Cohorts"]) 
bmi_df["Traits"] =  "BMI"
bmi_df["Category"] = "eQTL Informed"

bs_df = pd.read_csv(join(filt_path, "combined_Adipose_Subcutaneous_prs_0.1_BMI_SEX.txt"), sep="\t", names=["P-Value", "R2", "Cohorts"]) 
bs_df["Traits"] =  "BMI(Covars: SEX)"
bs_df["Category"] = "eQTL Informed"

bsa_df = pd.read_csv(join(filt_path, "combined_Adipose_Subcutaneous_prs_0.1_BMI_SEX_AGE.txt"), sep="\t", names=["P-Value", "R2", "Cohorts"]) 
bsa_df["Traits"] =  "BMI(Covars: SEX+AGE)"
bsa_df["Category"] = "eQTL Informed"

unfilt_path = "/Users/suryachhetri/datasets/prs_project/final/outputs/rsq_unfilt/plots"
bmi_unfilt_df = pd.read_csv(join(unfilt_path, "combined_Adipose_Subcutaneous_prs_unfilt_BMI.txt"), sep="\t", names=["P-Value", "R2", "Cohorts"])
bmi_unfilt_df["Traits"] =  "BMI"
bmi_unfilt_df["Category"] = "Baseline"

bs_unfilt_df = pd.read_csv(join(unfilt_path, "combined_Adipose_Subcutaneous_prs_unfilt_BMI_SEX.txt"), sep="\t", names=["P-Value", "R2", "Cohorts"])
bs_unfilt_df["Traits"] =  "BMI(Covars: SEX)"
bs_unfilt_df["Category"] = "Baseline"

bsa_unfilt_df = pd.read_csv(join(unfilt_path, "combined_Adipose_Subcutaneous_prs_unfilt_BMI_SEX_AGE.txt"), sep="\t", names=["P-Value", "R2", "Cohorts"])
bsa_unfilt_df["Traits"] =  "BMI(Covars: SEX+AGE)"
bsa_unfilt_df["Category"] = "Baseline"

dfs = [bmi_df, bs_df, bsa_df, bmi_unfilt_df, bs_unfilt_df, bsa_unfilt_df]
final_df = pd.concat(dfs, ignore_index=True)



"""BMI"""
bmi_df1 = bmi_df[~(bmi_df["Cohorts"] == "ALL")]
sns.barplot(x = 'P-Value', y = 'R2', hue = 'Cohorts', data = final_df,
            palette = 'hls',
            #col="Cohorts",
            #order = ['EA', 'AA', 'ALL'],  
            capsize = 0.05,             
            saturation = 8,             
            errcolor = 'gray', errwidth = 2,  
            ci = 'sd'   
            )

plot_path="/Users/suryachhetri/datasets/prs_project/final/outputs"
plt.savefig(join(plot_path, 'BMI_PRS_0.1.pdf'))
plt.close()



final_df_bmi = final_df[final_df["Traits"] == "BMI"]
baseline_df = final_df[final_df["Category"]=="Baseline"]
EA_AAonly_basedf = baseline_df[~(baseline_df["Cohorts"] == "ALL")]
plot_path = "/Users/suryachhetri/datasets/prs_project/final/outputs"

# sns.set_context('paper')
# sns.catplot(x = 'P-Value', y = 'R2', hue = 'Category', data = final_df_bmi,
#             palette = 'hls',
#             #col="Cohorts",
#             kind="bar",
#             height=4, aspect=1,
#             #order = ['EA', 'AA', 'ALL'],  
#             #capsize = 0.05,             
#             #saturation = 8,             
#             #errcolor = 'gray', errwidth = 2,  
#             #ci = 'sd'   
#             )
# plt.legend()


# sns.set_context('paper')
# sns.catplot(x = 'P-Value', y = 'R2', hue = 'Cohorts', data = EA_AAonly_basedf,
#             palette = 'hls',
#             col="Traits",
#             kind="bar",
#             height=4, aspect=1,
#             #order = ['EA', 'AA', 'ALL'],  
#             #capsize = 0.05,             
#             #saturation = 8,             
#             #errcolor = 'gray', errwidth = 2,  
#             #ci = 'sd'   
#             )
# plt.legend()

"""BMI"""
bmi_df1 = bmi_df[~(bmi_df["Cohorts"] == "ALL")]
sns.barplot(x = 'P-Value', y = 'R2', hue = 'Cohorts', data = bmi_df1,
            palette = 'hls',
            #col="Cohorts",
            #order = ['EA', 'AA', 'ALL'],  
            capsize = 0.05,             
            saturation = 8,             
            errcolor = 'gray', errwidth = 2,  
            ci = 'sd'   
            )

plot_path="/Users/suryachhetri/datasets/prs_project/final/outputs"
plt.savefig(join(plot_path, 'BMI_PRS_0.1.pdf'))
plt.close()

"""BMI"""
bs_df1 = bs_df[~(bs_df["Cohorts"] == "ALL")]
sns.barplot(x = 'P-Value', y = 'R2', hue = 'Cohorts', data = bs_df1,
            palette = 'hls',
            #col="Cohorts",
            #order = ['EA', 'AA', 'ALL'],  
            capsize = 0.05,             
            saturation = 8,             
            errcolor = 'gray', errwidth = 2,  
            ci = 'sd'   
            )

plot_path="/Users/suryachhetri/datasets/prs_project/final/outputs"
plt.savefig(join(plot_path, 'BMI_AGE_PRS_0.1.pdf'))
plt.close()

"""BMI"""
bsa_df1 = bsa_df[~(bsa_df["Cohorts"] == "ALL")]
sns.barplot(x = 'P-Value', y = 'R2', hue = 'Cohorts', data = bsa_df1,
            palette = 'hls',
            #col="Cohorts",
            #order = ['EA', 'AA', 'ALL'],  
            capsize = 0.05,             
            saturation = 8,             
            errcolor = 'gray', errwidth = 2,  
            ci = 'sd'   
            )

plot_path="/Users/suryachhetri/datasets/prs_project/final/outputs"
plt.savefig(join(plot_path, 'BMI_SEX_AGE_PRS_0.1.pdf'))
plt.close()

final_df1 = final_df[final_df["Category"] == "eQTL Informed"]
sns.catplot(x = 'P-Value', y = 'R2', hue = 'Traits', data = final_df1,
            palette = 'hls',
            col="Cohorts",
            kind="bar",
            height=4, aspect=1,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
plt.legend()
plt.savefig(join(plot_path, 'CombinedplotCohorts_BMI_AGE_PRS_0.1.pdf'))
plt.close()


final_df1 = final_df[final_df["Category"] == "eQTL Informed"]
sns.catplot(x = 'P-Value', y = 'R2', hue = 'Cohorts', data = final_df1,
            palette = 'hls',
            col="Traits",
            kind="bar",
            height=4, aspect=1,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
plt.legend()
plt.savefig(join(plot_path, 'CombinedplotTraits_BMI_AGE_PRS_0.1.pdf'))
plt.close()
    
# Put the legend out of the figure
# Put a legend to the right side
#plt.legend(loc='top', ncol=1)

# Plotting with baseline only:
baseline_df = final_df[final_df["Category"]=="Baseline"]
EA_AAonly_baselinedf = baseline_df[~(baseline_df["Cohorts"] == "ALL")]

sns.catplot(x = 'P-Value', y = 'R2', hue = 'Cohorts', data = EA_AAonly_baselinedf,
            palette = 'Set2',
            col="Traits",
            kind="bar",
            height=4, aspect=1,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
plt.legend()
plt.savefig(join(plot_path, 'Baseline_CombinedplotTraits_BMI_AGE_PRS_0.1.pdf'))
plt.close()

# Plotting the comparison between baseline and eQTL informed:
EA_AAonly = final_df[~(final_df["Cohorts"] == "ALL")]
final_EA_AAonly = EA_AAonly[~(EA_AAonly["P-Value"] == 3.000000e-01)] #filtering this weird date
bmi_final_EA_AA = final_EA_AAonly[final_EA_AAonly["Traits"] == "BMI"]
bs_final_EA_AA = final_EA_AAonly[final_EA_AAonly["Traits"] == "BMI(Covars: SEX)"]
bsa_final_EA_AA = final_EA_AAonly[final_EA_AAonly["Traits"] == "BMI(Covars: SEX+AGE)"]

sns.catplot(x = 'P-Value', y = 'R2', hue = 'Category', data = final_EA_AAonly,
            palette = 'Set2',
            col="Cohorts",
            kind="bar",
            height=4, aspect=1,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )

plt.legend()
plt.savefig(join(plot_path, 'Compare_eQTL_Baseline_CombinedplotTraits_BMI_AGE_PRS_0.1.pdf'))
plt.close()

# # Plotting the comparison between baseline and eQTL informed:
# EA_AAonly = final_df[~(final_df["Cohorts"] == "ALL")]
# final_EA_AAonly = EA_AAonly[~(EA_AAonly["P-Value"] == 3.000000e-01)] #filtering this weird date
# sns.catplot(x = 'P-Value', y = 'R2', hue = 'Category', data = final_EA_AAonly,
#             palette = 'Set3',
#             col="Traits",
#             kind="bar",
#             height=4, aspect=1,
#             #order = ['EA', 'AA', 'ALL'],  
#             #capsize = 0.05,             
#             #saturation = 8,             
#             #errcolor = 'gray', errwidth = 2,  
#             #ci = 'sd'   
#             )

# plt.legend()
# plt.savefig(join(plot_path, 'Compare_eQTL_Baseline_CombinedplotTraits_BMI_AGE_PRS_0.1.pdf'))
# plt.close()

#Individual level plot:
#BMI compare
sns.catplot(x = 'P-Value', y = 'R2', hue = 'Category', data = bmi_final_EA_AA,
            palette = 'Set2',
            col="Cohorts",
            kind="bar",
            height=4, aspect=1,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )

plt.legend()
plt.savefig(join(plot_path, 'BMICompare_eQTL_Baseline_CombinedplotTraits_BMI_AGE_PRS_0.1.pdf'))
plt.close()


#Individual level plot:
sns.catplot(x = 'P-Value', y = 'R2', hue = 'Category', data = bsa_final_EA_AA,
            palette = 'Set2',
            col="Cohorts",
            kind="bar",
            height=4, aspect=1,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )

plt.legend()
plt.savefig(join(plot_path, 'BMI_SEX_AGE_Compare_eQTL_Baseline_CombinedplotTraits_BMI_AGE_PRS_0.1.pdf'))
plt.close()

#############################
# Gene Datasets:
path="/Users/suryachhetri/datasets/prs_project/final"
genecorr_df = pd.read_csv("Adipose_Subcutaneous_shared_gene_corr_coeff_stats.05.txt", sep="\t")
tissue="Adipose"
""" boxplot for genomewide beta corr and coeff variation """
#df_result = result_df[2]
df_new = genecorr_df.loc[:, ["chrom", "Beta_rsquare", "Zscore_rsquare"]]
df_new = genecorr_df.loc[:, ["chrom", "Beta_corr", "Zscore_corr"]]
df_melt = pd.melt(df_new, "chrom", var_name="measurement")

""" using seaborn backend and placed in subplot"""
#sns.set(style="ticks")
sns.set_context('paper')
fig, ax = plt.subplots()

# Draw a nested boxplot:
# fig.set_size_inches(8, 6)
axs = sns.boxplot(x="chrom", y="value",
            hue="measurement", palette=["r", "b"],
            data=df_melt, \
            ax=ax, \
            #"""remove outlier"""
            flierprops = dict(markerfacecolor = '0.3', markersize = 1))

sns.despine(offset=5, trim=True)

# set axes and labels for sns subplot:
axs.set_xlabel("Chromosome", fontsize=10)
axs.set_ylabel("Correlation (n=15884 Unique Genes)", fontsize=10)
axs.set_title("Genome Wide Distribution (Beta|ZScore)", weight='bold', fontsize=12)
axs.tick_params(axis='x', labelrotation=45, labelcolor='black', labelsize=9)
axs.tick_params(axis='y', labelrotation=0, labelcolor='black', labelsize=9)

# Put the legend out of the figure
#ax.legend(loc=4, prop={'size': 10}, edgecolor="black")
#ax.axhline(0.8, color='tab:blue', linestyle='dashed', linewidth=2)

# Save figure:
fig.tight_layout()
fig.savefig(join(plot_path, 'Combinedplot_{}_genomewide_betaStats_boxplot.pdf'.format(tissue)))
plt.close()

