

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
tissue = "Adipose_Subcutaneous" #tissue = "Heart_Atrial"
chromrange = 22
chrom_X = False
pval = 10e-5 

output_path = "/Users/suryachhetri/datasets/prs_project/out_sharedvars"
plot_path = join(output_path, "plots")
plot_path = join(output_path, "plots_chrX")
ldgenome_path = ""

""" create plot output dir"""
if not os.path.exists(plot_path):
    os.makedirs(plot_path)

def pop_shared_variants_betabased(aa_filepath, ea_filepath, chromrange, pval, tissue, chrom_X=False):
    """ 
    Finds shared variants and specfic variants between 
    population cohorts without similar effects/eQTL info

    """
    ## Enter the chromosome no. range that you would like to analyse the data on. By default, it would take the autosomal 
    ## gene model from (chr1:chr22, excluding chrX & chrY), however, set chrom_XY=True for normal gene model analysis.
    chromrange = 1
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

    shared_variant_set1 = []
    shared_variant_set2 = []
    stats_filt_all = []
    aa_specific_df = []
    ea_specific_df = []
    genecount = []
    variantcount = []

    for chrom in chrom_list:
        print("\nProcessing : {} ...\n".format(chrom))

        """ parse afr cohorts """
        read_aa = pd.read_parquet(aa_dict.get(chrom), engine="fastparquet")
        read_aa.reset_index(inplace=True, drop=True)
        sorted_aa = read_aa.sort_values(["variant_id", "pval_nominal"], ascending=True).reset_index(drop=True)
        sorted_aa.drop_duplicates(["variant_id"], inplace=True, keep="first") #drop based on most sig pval

        """ parse eur cohorts"""
        read_ea = pd.read_parquet(ea_dict.get(chrom), engine="fastparquet")
        read_ea.reset_index(inplace=True, drop=True)
        sorted_ea = read_ea.sort_values(["variant_id", "pval_nominal"], ascending=True).reset_index(drop=True)
        sorted_ea.drop_duplicates(["variant_id"], inplace=True, keep="first") #drop based on most sig pval


                """ merge afr eur cohorts """
        merged_df = pd.merge(sorted_aa, sorted_ea, on=["variant_id", "phenotype_id"], how = "outer", suffixes=["_AA", "_EA"], indicator=True)
        shared_variants = merged_df[merged_df["_merge"] == "both"] #filter based on sig pval 0.0001 of at least in one cohort
        print("total shared variants across population without thresholds : {}".format(shared_variants.shape[0]))
        shared_vars_thres = shared_variants.loc[(shared_variants['pval_nominal_AA'] <= pval) | (shared_variants['pval_nominal_EA'] <= pval)]
        aa_specific = merged_df[merged_df["_merge"] == "left_only"]
        aa_specific_thresh = aa_specific.loc[aa_specific['pval_nominal_AA'] <= pval]
        ea_specific = merged_df[merged_df["_merge"] == "right_only"]
        ea_specific_thresh = ea_specific.loc[ea_specific['pval_nominal_EA'] <= pval]
        
        """Basic summary stats"""
        print('''\tTotal unique variants detected across cohorts ::
        European cohorts : {0}
        African cohorts : {1}'''.format(sorted_aa.shape[0], sorted_ea.shape[0]))
        print('''\t\nVariants dist without thresholds ::
        Shared : {0}
        AFR specific : {1}
        EUR specific : {2}'''.format(shared_variants.shape[0],aa_specific.shape[0],ea_specific.shape[0]))
        print('''\t\nVariants with pval {0} thresh ::
        Shared (at least one pop thresh) : {1}
        AFR specific : {2}
        EUR specific : {3}'''.format(pval, shared_vars_thres.shape[0],aa_specific_thresh.shape[0],ea_specific_thresh.shape[0]))

        """calculate correlation"""
        shared_vars_thres["logpval_nominal_AA"] = -(np.log10(shared_vars_thres["pval_nominal_AA"]))
        shared_vars_thres["logpval_nominal_EA"] = -(np.log10(shared_vars_thres["pval_nominal_EA"]))

        """calculate variant-gene pairwise slope corr(r) and coeff of var(rsquare) for each gene"""
        def scipystats_linregress(df): return stats.linregress(df["slope_AA"], df["slope_EA"])[2] #slope,intercept,r_value,p_value,std_err
        def scipystats_zscorelinregress(df): return stats.linregress(df["ZScore_AA"], df["ZScore_EA"])[2] #slope,intercept,r_value,p_value,std_err
        stats_df1 = merged_df1.groupby(["phenotype_id"]).apply(scipystats_zscorelinregress)

        stats_df = shared_vars_thres.groupby(["phenotype_id"]).apply(scipystats_linregress)
        stats_df = shared_vars_thres.groupby(["phenotype_id"]).apply(scipystats_linregress)
        stats_df = stats_df.reset_index(name="beta_correlation")
        stats_filt_df = stats_df[~(stats_df["beta_correlation"] == 0)].reset_index(drop=True)
        stats_filt_df["rsquare"] = stats_filt_df["beta_correlation"]**2
        stats_filt_df["chrom"] = chrom
        cols = stats_filt_df.columns.tolist()
        select_cols = cols[-1:] + cols[:-1] #Flip chrom to first pos
        stats_filt_df = stats_filt_df.loc[:,select_cols]
        stats_filt_all.append(stats_filt_df)

        """variant, genecount and median stats"""
        variantcount.append(shared_vars_thres.shape[0])
        genecount.append(stats_filt_df.shape[0])
        corr_median = stats_filt_df["beta_correlation"].median()
        rsquare_median = stats_filt_df["rsquare"].median()
        print("total genes, filtered genes :  {}, {}".format(stats_df.shape[0], stats_filt_df.shape[0]))

        def pandas_linregress_corr(df):
            corr = df.groupby("phenotype_id")[['slope_AA','slope_EA']].corr()
            df_corr = corr.iloc[0::2,-1]
            corr_df = df_corr.reset_index()
            nans = lambda df1: df1.loc[df1.isnull().any(axis=1)]
            nancorr_rows = nans(corr_df)
            corr_df.dropna(inplace=True)
            corr_df.columns = ["phenotype_id", "slope_level", "beta_correlation"]
            corr_df = corr_df.loc[:, ["phenotype_id", "beta_correlation"]]
            return corr_df

        """ calculate correlation using pandas"""
        # corr_df = pandas_linregress_corr(shared_vars_thres)

        """ plot histogram for correlation distribution """
        fig, ax = plt.subplots(1,1) # set up figure and axes; concise than : fig = plt.figure(); ax = fig.add_subplot(111)
        ax.hist(stats_filt_df["beta_correlation"], color="red", bins=20)
        plt.grid(True)
        plt.axvline(corr_median, color='b', linestyle='dashed', linewidth=2)
        plt.xlabel('corr(beta_AA vs beta_EA)  |  SharedVars(n = {})'.format(str(shared_vars_thres.shape[0])))
        plt.title('Beta Correlation for ' + str(chrom) + "(n = {} genes)".format(str(stats_filt_df.shape[0])))
        median = "Median : {}".format(str(round(corr_median, 2)))
        median_text = AnchoredText(median, loc=2)
        ax.add_artist(median_text)
        plt.savefig(join(plot_path, '{}_betaStats_corr.pdf'.format(chrom))) #saves current figure
        plt.close()

        fig, ax = plt.subplots(1,1)
        ax.hist(stats_filt_df["rsquare"], color="red", bins=20)
        plt.grid(True)
        plt.axvline(rsquare_median, color='b', linestyle='dashed', linewidth=2)
        plt.xlabel('coeff(beta_AA vs beta_EA)  |  SharedVars(n = {})'.format(str(shared_vars_thres.shape[0])))
        plt.title('Beta Rsquare for ' + str(chrom) + "(n = {} genes)".format(str(stats_filt_df.shape[0])))
        median = "Median : {}".format(str(round(rsquare_median, 2)))
        median_text = AnchoredText(median, loc=2)
        ax.add_artist(median_text)
        plt.savefig(join(plot_path, '{}_betaStats_coeff.pdf'.format(chrom))) #saves current figure
        plt.close()
        
        """ merge stats filt df for retention of variants and tss info"""
        shared_df = shared_vars_thres.reset_index(drop=True)

        #####Debug#####
        merged_df1 = pd.merge(neg_stats_df, shared_df, on=["phenotype_id"])
        select_cols = ["phenotype_id", "variant_id", "tss_distance_AA", "tss_distance_EA", "slope_AA", "slope_EA",
                        "logpval_nominal_AA", "logpval_nominal_EA", "beta_correlation", 'rsquare']
        shared_var_set1 = merged_df1.loc[:,select_cols]
        shared_var_set1[["tss_distance_AA", "tss_distance_EA"]] = shared_var_set1[["tss_distance_AA", "tss_distance_EA"]].astype(int)

        def scipystats_zscorelinregress(df): return stats.linregress(df["ZScore_AA"], df["ZScore_EA"])[2] #slope,intercept,r_value,p_value,std_err
        stats_df1 = merged_df1.groupby(["phenotype_id"]).apply(scipystats_zscorelinregress)

        ######Debug####

        merged_df1 = pd.merge(stats_filt_df, shared_df, on=["phenotype_id"])
        merged_df1["ZScore_AA"] = merged_df1["slope_AA"]/merged_df1["slope_se_AA"]
        merged_df1["ZScore_EA"] = merged_df1["slope_EA"]/merged_df1["slope_se_EA"]
        select_cols = ["phenotype_id", "variant_id", "tss_distance_AA", "tss_distance_EA", "slope_AA", "slope_EA",
                        "logpval_nominal_AA", "logpval_nominal_EA", "beta_correlation", 'rsquare']
        shared_var_set1 = merged_df1.loc[:,select_cols]
        shared_var_set1[["tss_distance_AA", "tss_distance_EA"]] = shared_var_set1[["tss_distance_AA", "tss_distance_EA"]].astype(int)
        shared_variant_set1.append(shared_var_set1)

        #filtered merged datasets: # coeff of var of at least 0.5
        stats_05 = stats_filt_df[stats_filt_df["rsquare"] >= 0.5]
        merged_df2 = pd.merge(stats_05, shared_df, on=["phenotype_id"])
        shared_var_set2 = merged_df2.loc[:,select_cols]
        shared_var_set2[["tss_distance_AA", "tss_distance_EA"]] = shared_var_set2[["tss_distance_AA", "tss_distance_EA"]].astype(int)
        shared_variant_set2.append(shared_var_set2)

        """ find variant to highlight """
        pheno_idx = shared_var_set1.groupby("phenotype_id").size().sort_values().idxmax()
        
        ###############
        variant_df = shared_var_set1[shared_var_set1["phenotype_id"] == pheno_idx] #pheno_idx = "ENSG00000001460.17"
        variant_df = merged_df1[merged_df1["phenotype_id"] == "ENSG00000097046.12"] #pheno_idx = "ENSG00000001460.17"
        ###############

        # generate scatter plot:
        cm = plt.cm.get_cmap('RdYlBu_r') #["autumn", "viridis"]
        plt.subplot(211)
        plt.scatter(x=variant_df["tss_distance_AA"], y=variant_df["logpval_nominal_AA"], 
                    c=variant_df["slope_AA"], cmap=cm, s=4, alpha=0.8)
        plt.colorbar().set_label('Beta value')
        text = AnchoredText("AFR Cohort", loc=1)
        plt.gca().add_artist(text)
        plt.ylabel("-Log10(pval)")
        #plt.legend(loc=1, prop={'size': 8}) #requires label="" in plt.scatter()

        plt.subplot(212)
        plt.scatter(x=variant_df["tss_distance_EA"], y=variant_df["logpval_nominal_EA"], 
                    c=variant_df["slope_EA"], cmap=cm, s=4, alpha=0.8)
        plt.colorbar().set_label('Beta value')
        text = AnchoredText("EUR Cohort", loc=1)
        plt.gca().add_artist(text)
        plt.xlabel("TSS Dist(bp)")
        plt.ylabel("-Log10(pval)")
        #plt.legend(loc=1, prop={'size': 8}) #requires label="" in plt.scatter()
        
        # super title
        plt.suptitle('''Shared Variants Dist (n = {0})
        {1}'''.format(variant_df.shape[0], str(pheno_idx)), fontsize=8, y=0.98, ha='center')
        plt.savefig(join(plot_path, '{}_debug__variant_scatter.pdf'.format("ENSG00000097046.12"))) #saves current figure
        plt.close()

    """summing stats and df for all chroms"""
    total_genecount = sum(genecount)
    total_variantcount = sum(variantcount)
    shared_set1 = pd.concat(shared_variant_set1, ignore_index=True)
    shared_set1.to_csv(join(output_path, tissue+"_shared_gene_variants_set1.txt"), sep="\t", index=False, header=True)
    shared_set2 = pd.concat(shared_variant_set2, ignore_index=True)
    shared_set2.to_csv(join(output_path, tissue+"_shared_gene_variants_set2_rsq.05.txt"), sep="\t", index=False, header=True)
    stats_filt_comp = pd.concat(stats_filt_all, ignore_index=True)
    stats_filt_comp.to_csv(join(output_path, tissue+"_shared_gene_corr_coeff_stats.05.txt"), sep="\t", index=False, header=True)


    """ plot histogram for correlation distribution """
    corr_median_comp = round(stats_filt_comp["beta_correlation"].median(), 2)
    rsquare_median_comp = round(stats_filt_comp["rsquare"].median(), 2)
    print("\n\n===========================================================\n\n")
    print("Total Genes, Total SharedVariants :  {}, {}\n".format(total_genecount, total_variantcount))
    print("Average Genome wide beta Corr, beta coeff :  {}, {}\n".format(corr_median_comp, rsquare_median_comp))
    print("\n\n===========================================================\n\n")

    fig, ax = plt.subplots(1,1) # set up figure and axes; concise than : fig = plt.figure(); ax = fig.add_subplot(111)
    ax.hist(stats_filt_comp["beta_correlation"], color="red", bins=20)
    plt.grid(True)
    plt.axvline(corr_median_comp, color='b', linestyle='dashed', linewidth=2)
    plt.xlabel('corr(beta_AA vs beta_EA)  |  SharedVars(n = {})'.format(str(total_variantcount)))
    plt.title('Avg Beta Correlation Genome-wide' + "(n = {} genes)".format(str(total_genecount)))
    median = "Median : {}".format(str(round(corr_median_comp, 2)))
    median_text = AnchoredText(median, loc=2)
    ax.add_artist(median_text)
    plt.savefig(join(plot_path, '{}_allchrom_betaStats_corr.pdf'.format(tissue))) #saves current figure
    plt.close()

    fig, ax = plt.subplots(1,1)
    ax.hist(stats_filt_comp["rsquare"], color="red", bins=20)
    plt.grid(True)
    plt.axvline(rsquare_median_comp, color='b', linestyle='dashed', linewidth=2)
    plt.xlabel('coeff(beta_AA vs beta_EA)  |  SharedVars(n = {})'.format(str(total_variantcount)))
    plt.title('Avg Beta Rsquare Genome-wide' + "(n = {} genes)".format(str(total_genecount)))
    median = "Median : {}".format(str(round(rsquare_median_comp, 2)))
    median_text = AnchoredText(median, loc=2)
    ax.add_artist(median_text)
    plt.savefig(join(plot_path, '{}_allchrom_betaStats_coeff.pdf'.format(tissue))) #saves current figure
    plt.close()

    return(shared_set1, shared_set2, stats_filt_comp)

"""function call"""
if __name__ == "__main__":
    #import warnings
    warnings.filterwarnings("ignore")
    start_time = time.time()
    chromrange = 1
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



###########################################################################
###########################################################################

def prs_calc():

shared_vars = "/Users/suryachhetri/datasets/prs_project/out_sharedvars/Adipose_Subcutaneous_shared_gene_variants_set1.txt"
gwas_stats_forBeta[i] = "/Users/suryachhetri/datasets/prs_project/gwas_stats/bmi/21001_irnt.gwas.imputed_v3.both_sexes.tsv"
lifted_gwas_stats = ""

genotype_forGeno[ij] = "vcf|plink .pvar" #(and find 0/1/2 cols)


ldgenome_path = ""
lifted_gwas_stats = ""
snp_list_for_ldclump = ""






