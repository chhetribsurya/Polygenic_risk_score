

import numpy as np, pandas as pd
import pybedtools
import os, re
from glob import glob
from os.path import basename, join, splitext
from collections import defaultdict
import time, warnings
import matplotlib; #matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy import stats
import seaborn as sns
    

""" Input parameters """
aa_filepath = "/Users/suryachhetri/datasets/prs_project/gtex/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_AFR_eQTL_all_associations"
ea_filepath = "/Users/suryachhetri/datasets/prs_project/gtex/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations"
gwas_stats = "/Users/suryachhetri/datasets/prs_project/gwas_stats/bmi/21001_irnt.gwas.imputed_v3.both_sexes.tsv"
lifted_gwas_stats = "/Users/suryachhetri/datasets/prs_project/gwas_stats/bmi/lifted/updated_output.bed"
tissue = "Adipose_Subcutaneous" #tissue = "Heart_Atrial"
chromrange = 22
chrom_X = False
pval = 10e-5

output_path = "/Users/suryachhetri/datasets/prs_project/out_sharedvars"
plot_path = join(output_path, "plots")
plot_path = join(output_path, "plots_chrX")
plot_path = join(output_path, "redo_chr1")

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

        """Mean-centered scaling - adjusting marked distribution diff with mean 0 per std dev 1"""
        def zscore_scaling(s):
            return (s - s.mean())/s.std()

        shared_variants_df.reset_index(drop=True, inplace=True)
        shared_variants_df["scaled_Z_EA"] = shared_variants_df.groupby(["phenotype_id"])["ZScore_EA"].transform(zscore_scaling)
        shared_variants_df["scaled_Z_AA"] = shared_variants_df.groupby(["phenotype_id"])["ZScore_AA"].transform(zscore_scaling)

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
        stats_wilcoxon = shared_variants_df.groupby(["phenotype_id"]).apply(stats_wilcoxon)
        stats_wilcoxontest = stats_wilcoxon.reset_index(name="wilcoxon_score")
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
        stats_filt_df[stats_filt_df["phenotype_id"] == "ENSG00000280670.2"].loc[:, 
            ["KS_score_adj", "KS_score_adj_log", "wilcoxon_score_adj", "wilcoxon_score_adj_log"]]

        stats_filt_df["chrom"] = chrom
        cols = stats_filt_df.columns.tolist()
        select_cols = cols[-1:] + cols[:-1] #Flip chrom to first pos
        stats_filt_df = stats_filt_df.loc[:,select_cols]
        stats_filt_all.append(stats_filt_df)

        """variant, genecount and median stats"""
        variantcount.append(shared_variants_df.shape[0])
        genecount.append(stats_filt_df.shape[0])
        corr_median = stats_filt_df["Beta_corr"].median()
        rsquare_median = stats_filt_df["Beta_rsquare"].median()
        Z_rsquare_median = stats_filt_df["Zscore_corr"].median()
        Z_corr_median = stats_filt_df["Zscore_rsquare"].median()

        print("Total filtered genes :  {}".format(stats_filt_df.shape[0]))

        ##############################
        # Merge dataframes:
        ##############################
        select_cols = ["phenotype_id", "Zscore_rsquare","wilcoxon_score_adj_log"]
        genewise_select = stats_filt_df.loc[:, select_cols]

        select_cols = ["phenotype_id", "variant_id", "pval_nominal_EA", "pval_nominal_EA", "scaled_Z_EA", "scaled_Z_AA"]
        shared_vars_select = shared_variants_df.loc[:, select_cols]
        shared_vars_select["newsnp"] = shared_vars_select["variant_id"].str.replace("_b38", "").str.replace("_", ":")

        merged_eqtl_stats = pd.merge(shared_vars_select, genewise_select, on="phenotype_id", how="inner")

        gwas_df = pd.read_csv(gwas_stats, sep="\t")
        select_cols = ["variant", "beta", "pval"]
        gwas_select = gwas_df.loc[:,select_cols]

        lifted_gwas = pd.read_csv(lifted_gwas_stats, sep="\t", header=None)
        lifted_gwas.columns = ["chrom", "start", "end", "oldsnp", "newsnp"]
        lifted_gwas["newsnp"] =  "chr" + lifted_gwas["newsnp"].astype(str)
        lifted_gwas.reset_index(drop=True, inplace=True)

        merged_gwas = pd.merge(gwas_select, lifted_gwas, left_on=["variant"], right_on=["oldsnp"], how="inner")
        #updated_gwas = merged_gwas.reset_index(drop=True)

        select_cols = ["newsnp", "beta", "pval"]
        updated_gwas = merged_gwas.loc[:,select_cols]

        """merge"""
        merged_eqtl_statsfilt = merged_eqtl_stats.loc[:,["newsnp", "phenotype_id", "Zscore_rsquare", "wilcoxon_score_adj_log"]]
        #merged_df = pd.merge(shared_vars_select, updated_gwas_df, on="newsnp", how="outer")
        gwas_eqtl = pd.merge(merged_eqtl_statsfilt, updated_gwas, on="newsnp", how="inner")
        gwas_eqtl.loc[:, ["phenotype_id", "newsnp", "Zscore_rsquare", "wilcoxon_score_adj", "wilcoxon_score_adj_log"]]

        df1 = gwas_eqtl.sort_values("Zscore_rsquare", ascending=False)
        grouped = df1.groupby("newsnp")
        output1 = grouped.head(1).reset_index(drop=True)

        df2 = gwas_eqtl.sort_values("wilcoxon_score_adj_log", ascending=True)
        grouped = df2.groupby("newsnp")
        output_1 = grouped.head(1).reset_index(drop=True)

        output_2 = grouped.size().reset_index(name="gene_count")
        merged_df = pd.merge(output1, output_1, on=["newsnp", "beta", "pval"], how="inner")
        final_df = pd.merge(merged_df, output_2, on="newsnp", how="inner")

        #Apply hard-filtering:
        filter_df = final_df[final_df["Zscore_rsquare_x"] >= 0.1]
        prs_df = filter_df.loc[:, ["newsnp", "beta", "pval"]]
        prs_df1 = prs_df.copy(); 
        #prs_df1["newsnp"] = prs_df1["newsnp"].str.replace(":", "_")
        #prs_df2 = prs_df.copy()

        prs_df1[["chr1", "pos", "ref", "alt"]] = prs_df1["newsnp"].str.split(":", expand=True)
        df_test = prs_df1.loc[:, ["newsnp", "ref", "alt", "beta", "pval"]]
        df_test.rename(columns={"newsnp" : "ID", "ref" : "REF", "alt" : "ALT"}, inplace=True)

        #output_2 = grouped.size().reset_index(name="genecount")
        #chr1:758351:A:G   ENSG00000235098.8        0.305608                1.648058 -0.001071  0.781921

        """ plot histogram for correlation distribution """
        fig, ax = plt.subplots(1,1) # set up figure and axes; concise than : fig = plt.figure(); ax = fig.add_subplot(111)
        ax.hist(stats_filt_df["Beta_corr"], color="red", bins=20)
        plt.grid(True)
        plt.axvline(corr_median, color='b', linestyle='dashed', linewidth=2)
        plt.xlabel('corr(beta_AA vs beta_EA)  |  SharedVars(n = {})'.format(str(shared_variants_df.shape[0])))
        plt.title('Beta Correlation for ' + str(chrom) + "(n = {} genes)".format(str(stats_filt_df.shape[0])))
        median = "Median : {}".format(str(round(corr_median, 2)))
        median_text = AnchoredText(median, loc=2)
        ax.add_artist(median_text)
        plt.savefig(join(plot_path, '{}_betaStats_corr.pdf'.format(chrom))) #saves current figure
        plt.close()

        fig, ax = plt.subplots(1,1)
        ax.hist(stats_filt_df["Beta_rsquare"], color="red", bins=20)
        plt.grid(True)
        plt.axvline(rsquare_median, color='b', linestyle='dashed', linewidth=2)
        plt.xlabel('coeff(beta_AA vs beta_EA)  |  SharedVars(n = {})'.format(str(shared_variants_df.shape[0])))
        plt.title('Beta Rsquare for ' + str(chrom) + "(n = {} genes)".format(str(stats_filt_df.shape[0])))
        median = "Median : {}".format(str(round(rsquare_median, 2)))
        median_text = AnchoredText(median, loc=2)
        ax.add_artist(median_text)
        plt.savefig(join(plot_path, '{}_betaStats_coeff.pdf'.format(chrom))) #saves current figure
        plt.close()

        """ plot histogram for correlation distribution """
        fig, ax = plt.subplots(1,1) # set up figure and axes; concise than : fig = plt.figure(); ax = fig.add_subplot(111)
        ax.hist(stats_filt_df["Zscore_corr"], color="red", bins=20)
        plt.grid(True)
        plt.axvline(corr_median, color='blue', linestyle='dashed', linewidth=2)
        plt.xlabel('corr(Zscore_AA vs Zscore_EA)  |  SharedVars(n = {})'.format(str(shared_variants_df.shape[0])))
        plt.title('Zscore Correlation for ' + str(chrom) + "(n = {} genes)".format(str(stats_filt_df.shape[0])))
        median = "Median : {}".format(str(round(Z_rsquare_median, 2)))
        median_text = AnchoredText(median, loc=2)
        ax.add_artist(median_text)
        plt.savefig(join(plot_path, '{}_ZscoreStats_corr.pdf'.format(chrom))) #saves current figure
        plt.close()

        fig, ax = plt.subplots(1,1)
        ax.hist(stats_filt_df["Zscore_rsquare"], color="red", bins=20)
        plt.grid(True)
        plt.axvline(rsquare_median, color='blue', linestyle='dashed', linewidth=2)
        plt.xlabel('coeff(Zscore_AA vs Zscore_EA)  |  SharedVars(n = {})'.format(str(shared_variants_df.shape[0])))
        plt.title('Zscore Rsquare for ' + str(chrom) + "(n = {} genes)".format(str(stats_filt_df.shape[0])))
        median = "Median : {}".format(str(round(Z_rsquare_median, 2)))
        median_text = AnchoredText(Z_rsquare_median, loc=2)
        ax.add_artist(median_text)
        plt.savefig(join(plot_path, '{}_ZscoreStats_coeff.pdf'.format(chrom))) #saves current figure
        plt.close()

        # ##############
        # #Scatter plots (Solo)
        samp_1 = stats_filt_df[stats_filt_df["Beta_corr"] > 0.2].sample(1, random_state=1)["phenotype_id"].tolist()[0]
        samp_2 = stats_filt_df[stats_filt_df["Beta_corr"] < -0.2].sample(1, random_state=1)["phenotype_id"].tolist()[0]
        samp_3 = stats_filt_df.loc[stats_filt_df["Beta_corr"].idxmax()]["phenotype_id"]

        #Extract the variants to analyze on:
        samp_1_df = shared_variants_df[shared_variants_df["phenotype_id"] == samp_1]
        samp_2_df = shared_variants_df[shared_variants_df["phenotype_id"] == samp_2]
        samp_3_df = shared_variants_df[shared_variants_df["phenotype_id"] == samp_3]

        #####################
        # ### KS test:
        # ks_2_sample_test  = stats.ks_2samp(samp_1_df["ZScore_EA"], samp_1_df["ZScore_AA"])
        # ks_pval_prom = ks_2_sample_test[1]

        # if ks_pval_prom == 0:
        #     new_ks_pval_prom = 1e-323
        #     print "KS pvalue for promoter assoc :", new_ks_pval_prom
        # else:
        #     new_ks_pval_prom = ks_pval_prom
        #     print "KS pvalue for promoter assoc :", new_ks_pval_prom

        import numpy as np, pandas as pd
        import seaborn as sns
        sns.set(style="ticks", color_codes=True)
        from scipy import stats

        def r(x, y):
            return stats.pearsonr(x, y)[0]


        def r2(x, y):
            return stats.pearsonr(x, y)[0] ** 2


        sns.jointplot(x="ZScore_AA", y="ZScore_EA", data=samp_1_df, kind="reg", stat_func=r, color="red", 
                    scatter_kws={"s": 15, "alpha" : 0.6}, height=6)
        
        plt.tight_layout()
        plt.savefig(join(plot_path, '{}_Gene1_Scatter_histogram.pdf'.format(chrom)))
        plt.close()


        sns.jointplot(x="ZScore_AA", y="ZScore_EA", data=samp_2_df, kind="reg", stat_func=r, color="red", 
                    scatter_kws={"s": 15, "alpha" : 0.6}, height=6)
        
        plt.tight_layout()
        plt.savefig(join(plot_path, '{}_Gene2_Scatter_histogram.pdf'.format(chrom)))
        plt.close()


        sns.jointplot(x="ZScore_AA", y="ZScore_EA", data=samp_3_df, kind="reg", stat_func=r, color="red", 
                    scatter_kws={"s": 15, "alpha" : 0.6}, height=6)
        
        plt.tight_layout()
        plt.savefig(join(plot_path, '{}_Gene3_Scatter_histogram.pdf'.format(chrom)))
        plt.close()


        sns.jointplot(x="Zscore_corr", y="KS_score_adj_log", data=stats_filt_df, kind="reg", stat_func=r, color="red", 
                    scatter_kws={"s": 15, "alpha" : 0.6}, height=6)
        plt.axhline(5, color='black', linestyle='dashed', linewidth=1)
        
        plt.tight_layout()
        plt.savefig(join(plot_path, 'KS_adj_ZScoreCorr_scatter.pdf'.format(chrom)))
        plt.close()


        sns.jointplot(x="Zscore_rsquare", y="wilcoxon_score_adj_log", data=stats_filt_df, kind="reg", stat_func=r, color="orangered", 
                    scatter_kws={"s": 15, "alpha" : 0.6}, height=6)

        plt.axhline(-np.log10(3.125e-05), color='black', linestyle='dashed', linewidth=1)
        
        plt.tight_layout()
        plt.savefig(join(plot_path, 'Wilcoxon_adj_log_ZScoreRsquare_scatter.pdf'.format(chrom)))
        plt.close()

        sns.jointplot(x="Spearman_corr", y="wilcoxon_score_log", data=stats_filt_df, kind="reg", stat_func=r, color="orangered", 
                    scatter_kws={"s": 15, "alpha" : 0.6}, height=6)
        
        plt.tight_layout()
        plt.savefig(join(plot_path, '{}_Wilcoxon_SpearmanCorr_scatter.pdf'.format(chrom)))
        plt.close()

        # sns.jointplot(x="Beta_corr", y="wilcoxon_score_log", data=stats_filt_df, kind="reg", stat_func=r, color="orangered", 
        #             scatter_kws={"s": 15, "alpha" : 0.6}, height=6)
        
        # plt.tight_layout()
        # plt.savefig(join(plot_path, '{}_Wilcoxon_SpearmanCorr_scatter.pdf'.format(chrom)))
        # plt.close()

        sns.jointplot(x="Pearson_corr", y="Beta_corr", data=stats_filt_df, kind="reg", stat_func=r, color="orangered", 
                    scatter_kws={"s": 15, "alpha" : 0.6}, height=6)
        
        plt.tight_layout()
        plt.savefig(join(plot_path, '{}_PearsonoCorr_ZScoreCorr_scatter.pdf'.format(chrom)))
        plt.close()

        #sns.despine(offset=5, trim=True)
        # plt.hist(x=stats_filt_df["KS_score_log"], color="r", bins=150)

        # In [190]: plt.gca().spines['right'].set_visible(False)

        # In [191]: plt.gca().spines['top'].set_visible(False)

        # In [192]: plt.gca().spines['left'].set_smart_bounds(True)

        # In [193]: plt.gca().spines['bottom'].set_smart_bounds(True)

        # In [194]: plt.gca().spines['right'].set_smart_bounds(True)

        # In [195]: plt.gca().xaxis.set_ticks_position('bottom')

        # In [196]: plt.gca().yaxis.set_ticks_position('left')


        #####################

        plt.scatter(y=samp_1_df["slope_EA"], x=samp_1_df["slope_AA"], c=samp_1_df["tss_dist_mb"], s=4, alpha=0.7)
        plt.colorbar().set_label('TSS Dist (Mb)')
        plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
        plt.axhline(0, color='black', linestyle='dashed', linewidth=1)
        plt.ylabel('Slope EA')
        plt.xlabel('Slope AA (n = {} variants)'.format(samp_1_df.shape[0]))
        plt.title('Gene : {} ({})'.format(samp_1, chrom))
        plt.savefig(join(plot_path, '{}_Gene1_BetaScatter_tss.pdf'.format(chrom)))
        plt.close()

        plt.scatter(y=samp_2_df["slope_EA"], x=samp_2_df["slope_AA"], c=samp_2_df["tss_dist_mb"], s=4, alpha=0.7)
        plt.colorbar().set_label('TSS Dist (Mb)')
        plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
        plt.axhline(0, color='black', linestyle='dashed', linewidth=1)
        plt.ylabel('Slope EA')
        plt.xlabel('Slope AA (n = {} variants)'.format(samp_2_df.shape[0]))
        plt.title('Gene : {} ({})'.format(samp_2, chrom))
        plt.savefig(join(plot_path, '{}_Gene2_BetaScatter_tss.pdf'.format(chrom)))
        plt.close()

        plt.scatter(y=samp_3_df["slope_EA"], x=samp_3_df["slope_AA"], c=samp_3_df["tss_dist_mb"], s=4, alpha=0.7)
        plt.colorbar().set_label('TSS Dist (Mb)')
        plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
        plt.axhline(0, color='black', linestyle='dashed', linewidth=1)
        plt.ylabel('Slope EA')
        plt.xlabel('Slope AA (n = {} variants)'.format(samp_3_df.shape[0]))
        plt.title('Gene : {} ({})'.format(samp_3, chrom))
        plt.savefig(join(plot_path, '{}_Gene3_BetaScatter_tss.pdf'.format(chrom)))
        plt.close()

        ## Legend with pval ratio:
        plt.scatter(y=samp_1_df["slope_EA"], x=samp_1_df["slope_AA"], c=samp_1_df["pval_nominal_AA"], s=4, alpha=0.7)
        plt.colorbar().set_label('P Value AA')
        plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
        plt.axhline(0, color='black', linestyle='dashed', linewidth=1)
        plt.ylabel('Slope EA')
        plt.xlabel('Slope AA (n = {} variants)'.format(samp_1_df.shape[0]))
        plt.title('Gene : {} ({})'.format(samp_1, chrom))
        plt.savefig(join(plot_path, '{}_Gene1_BetaScatter_PVal_AA.pdf'.format(chrom)))
        plt.close()

        plt.scatter(y=samp_2_df["slope_EA"], x=samp_2_df["slope_AA"], c=samp_2_df["pval_nominal_AA"], s=4, alpha=0.7)
        plt.colorbar().set_label('P Value AA')
        plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
        plt.axhline(0, color='black', linestyle='dashed', linewidth=1)
        plt.ylabel('Slope EA')
        plt.xlabel('Slope AA (n = {} variants)'.format(samp_2_df.shape[0]))
        plt.title('Gene : {} ({})'.format(samp_2, chrom))
        plt.savefig(join(plot_path, '{}_Gene2_BetaScatter_PVal_AA.pdf'.format(chrom)))
        plt.close()

        plt.scatter(y=samp_3_df["slope_EA"], x=samp_3_df["slope_AA"], c=samp_3_df["pval_nominal_AA"], s=4, alpha=0.7)
        plt.colorbar().set_label('P Value AA')
        plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
        plt.axhline(0, color='black', linestyle='dashed', linewidth=1)
        plt.ylabel('Slope EA')
        plt.xlabel('Slope AA (n = {} variants)'.format(samp_3_df.shape[0]))
        plt.title('Gene : {} ({})'.format(samp_3, chrom))
        plt.savefig(join(plot_path, '{}_Gene3_BetaScatter_PVal_AA.pdf'.format(chrom)))
        plt.close()

        ## ZScore based scatter plots:
        plt.scatter(y=samp_1_df["ZScore_EA"], x=samp_1_df["ZScore_AA"], c=samp_1_df["pval_nominal_AA"], s=4, alpha=0.7)
        plt.colorbar().set_label('P Value AA')
        plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
        plt.axhline(0, color='black', linestyle='dashed', linewidth=1)
        plt.ylabel('ZScore EA')
        plt.xlabel('ZScore AA (n = {} variants)'.format(samp_1_df.shape[0]))
        plt.title('Gene : {} ({})'.format(samp_1, chrom))
        plt.savefig(join(plot_path, '{}_Gene1_ZscoreScatter_PVal_AA.pdf'.format(chrom)))
        plt.close()

        plt.scatter(y=samp_2_df["ZScore_EA"], x=samp_2_df["ZScore_AA"], c=samp_2_df["pval_nominal_AA"], s=4, alpha=0.7)
        plt.colorbar().set_label('P Value AA')
        plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
        plt.axhline(0, color='black', linestyle='dashed', linewidth=1)
        plt.ylabel('ZScore EA')
        plt.xlabel('ZScore AA (n = {} variants)'.format(samp_2_df.shape[0]))
        plt.title('Gene : {} ({})'.format(samp_2, chrom))
        plt.savefig(join(plot_path, '{}_Gene2_ZscoreScatter_PVal_AA.pdf'.format(chrom)))
        plt.close()

        plt.scatter(y=samp_3_df["ZScore_EA"], x=samp_3_df["ZScore_AA"], c=samp_3_df["pval_nominal_AA"], s=4, alpha=0.7)
        plt.colorbar().set_label('P Value AA')
        plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
        plt.axhline(0, color='black', linestyle='dashed', linewidth=1)
        plt.ylabel('ZScore EA')
        plt.xlabel('ZScore AA (n = {} variants)'.format(samp_3_df.shape[0]))
        plt.title('Gene : {} ({})'.format(samp_3, chrom))
        plt.savefig(join(plot_path, '{}_Gene3_ZscoreScatter_PVal_AA.pdf'.format(chrom)))
        plt.close()

        #############
        # Samp1 : Scatter plots combined:
        cm = plt.cm.get_cmap('RdYlBu_r') #["autumn", "viridis"]

        plt.subplot(211)
        plt.scatter(x=samp_1_df["tss_dist_mb"], y=samp_1_df["logpval_nominal_AA"], 
                    c=samp_1_df["slope_AA"], cmap=cm, s=4, alpha=0.8)
        plt.colorbar().set_label('Beta value')
        text = AnchoredText("AFR Cohort", loc=2)
        plt.gca().add_artist(text)
        plt.ylabel("-Log10(pval)")
        plt.ylim((0,5))
        plt.axhline(2, color='black', linestyle='dashed', linewidth=1)

        plt.subplot(212)
        plt.scatter(x=samp_1_df["tss_dist_mb"], y=samp_1_df["logpval_nominal_EA"], 
                    c=samp_1_df["slope_EA"], cmap=cm, s=4, alpha=0.8)
        plt.colorbar().set_label('Beta value')
        text = AnchoredText("EUR Cohort", loc=2)
        plt.gca().add_artist(text)
        plt.ylabel("-Log10(pval)")
        plt.ylim((0,5))
        plt.axhline(2, color='black', linestyle='dashed', linewidth=1)

        #Common X axis:
        plt.xlabel("TSS Dist (Mbp)")
        #plt.legend(loc=1, prop={'size': 8}) #requires label="" in plt.scatter()
        
        # super title
        plt.suptitle('''Shared Variants Dist (n = {0})
        {1} ({2})'''.format(samp_1_df.shape[0], samp_1, chrom), fontsize=8, y=0.98, ha='center')
        plt.savefig(join(plot_path, 'Gene1_YScaled_variant_TSS_scatter_AA_EA.pdf')) #saves current figure
        plt.close()

        # Samp 2: Scatter plots combined:
        cm = plt.cm.get_cmap('RdYlBu_r') #["autumn", "viridis"]
        
        plt.subplot(211)
        plt.scatter(x=samp_2_df["tss_dist_mb"], y=samp_2_df["logpval_nominal_AA"], 
                    c=samp_2_df["slope_AA"], cmap=cm, s=4, alpha=0.8)
        plt.colorbar().set_label('Beta value')
        text = AnchoredText("AFR Cohort", loc=2)
        plt.gca().add_artist(text)
        plt.ylabel("-Log10(pval)")
        plt.ylim((0,5))
        plt.axhline(2, color='black', linestyle='dashed', linewidth=1)

        plt.subplot(212)
        plt.scatter(x=samp_2_df["tss_dist_mb"], y=samp_2_df["logpval_nominal_EA"], 
                    c=samp_2_df["slope_EA"], cmap=cm, s=4, alpha=0.8)
        plt.colorbar().set_label('Beta value')
        text = AnchoredText("EUR Cohort", loc=2)
        plt.gca().add_artist(text)
        plt.ylabel("-Log10(pval)")
        plt.ylim((0,5))
        plt.axhline(2, color='black', linestyle='dashed', linewidth=1)

        #Common X axis:
        plt.xlabel("TSS Dist (Mbp)")
        #plt.legend(loc=1, prop={'size': 8}) #requires label="" in plt.scatter()
        
        # super title
        plt.suptitle('''Shared Variants Dist (n = {0})
        {1} ({2})'''.format(samp_1_df.shape[0], samp_1, chrom), fontsize=8, y=0.98, ha='center')
        plt.savefig(join(plot_path, 'Gene2_YScaled_variant_TSS_scatter_AA_EA.pdf')) #saves current figure
        plt.close()

        # Samp 3: Scatter plots combined:
        cm = plt.cm.get_cmap('RdYlBu_r') #["autumn", "viridis"]
        
        plt.subplot(211)
        plt.scatter(x=samp_3_df["tss_dist_mb"], y=samp_3_df["logpval_nominal_AA"], 
                    c=samp_3_df["slope_AA"], cmap=cm, s=4, alpha=0.8)
        plt.colorbar().set_label('Beta value')
        text = AnchoredText("AFR Cohort", loc=2)
        plt.gca().add_artist(text)
        plt.ylabel("-Log10(pval)")
        plt.ylim((0,15))
        plt.axhline(2, color='black', linestyle='dashed', linewidth=1)

        plt.subplot(212)
        plt.scatter(x=samp_3_df["tss_dist_mb"], y=samp_3_df["logpval_nominal_EA"], 
                    c=samp_3_df["slope_EA"], cmap=cm, s=4, alpha=0.8)
        plt.colorbar().set_label('Beta value')
        text = AnchoredText("EUR Cohort", loc=2)
        plt.gca().add_artist(text)
        plt.ylabel("-Log10(pval)")
        plt.ylim((0,15))
        plt.axhline(2, color='black', linestyle='dashed', linewidth=1)

        #Common X axis:
        plt.xlabel("TSS Dist (Mbp)")
        #plt.legend(loc=1, prop={'size': 8}) #requires label="" in plt.scatter()
        
        # super title
        plt.suptitle('''Shared Variants Dist (n = {0})
        {1} ({2})'''.format(samp_3_df.shape[0], samp_3, chrom), fontsize=8, y=0.98, ha='center')
        plt.savefig(join(plot_path, 'Gene3_re_Yscaled_variant_TSS_scatter_AA_EA.pdf')) #saves current figure
        plt.close()

        """ Variant to highlight with highest number of variants"""
        pheno_idx = shared_variants_df.groupby("phenotype_id").size().sort_values().idxmax()
        variant_df = shared_variants_df[shared_variants_df["phenotype_id"] == pheno_idx] #pheno_idx = "ENSG00000001460.17"


        # Gene with highest number of variants: Scatter plots combined:
        cm = plt.cm.get_cmap('RdYlBu_r') #["autumn", "viridis"]

        plt.subplot(211)
        plt.scatter(x=variant_df["tss_dist_mb"], y=variant_df["logpval_nominal_AA"], 
                    c=variant_df["slope_AA"], cmap=cm, s=4, alpha=0.8)
        plt.colorbar().set_label('Beta value')
        text = AnchoredText("AFR Cohort", loc=2)
        plt.gca().add_artist(text)
        plt.ylabel("-Log10(pval)")
        #plt.ylim((0,110))
        plt.axhline(2, color='black', linestyle='dashed', linewidth=1)

        plt.subplot(212)
        plt.scatter(x=variant_df["tss_dist_mb"], y=variant_df["logpval_nominal_EA"], 
                    c=variant_df["slope_EA"], cmap=cm, s=4, alpha=0.8)
        plt.colorbar().set_label('Beta value')
        text = AnchoredText("EUR Cohort", loc=2)
        plt.gca().add_artist(text)
        plt.ylabel("-Log10(pval)")
        #plt.ylim((0,110))
        plt.axhline(2, color='black', linestyle='dashed', linewidth=1)

        #Common X axis:
        plt.xlabel("TSS Dist (Mbp)")
        #plt.legend(loc=1, prop={'size': 8}) #requires label="" in plt.scatter()
        
        # super title
        plt.suptitle('''Shared Variants Dist (n = {0})
        {1} ({2})'''.format(variant_df.shape[0], pheno_idx, chrom), fontsize=8, y=0.98, ha='center')
        plt.savefig(join(plot_path, 'Gene_withhighest_numVariants_TSS_scatter_AA_EA.pdf')) #saves current figure
        plt.close()

        ######################
        # Seaborn plots:
        
        # ## Legend with pval ratio:
        # plt.scatter(y=samp_1_df["slope_EA"], x=samp_1_df["slope_AA"], c=samp_1_df["pval_nominal_AA"], s=4, alpha=0.7)
        # plt.colorbar().set_label('P Value AA')
        # plt.axvline(0, color='black', linestyle='dashed', linewidth=1)
        # plt.axhline(0, color='black', linestyle='dashed', linewidth=1)
        # plt.ylabel('Slope EA')
        # plt.xlabel('Slope AA (n = {} variants)'.format(samp_1_df.shape[0]))
        # plt.title('Gene : {} ({})'.format(samp_1, chrom))
        # plt.savefig(join(plot_path, '{}_Gene1_BetaScatter_PVal_AA.pdf'.format(chrom)))
        # plt.close()

        import numpy as np, pandas as pd
        import seaborn as sns
        sns.set(style="ticks", color_codes=True)
        from scipy import stats

        def r(x, y):
            return stats.pearsonr(x, y)[0]

        def r2(x, y):
            return stats.pearsonr(x, y)[0] ** 2

        sns.jointplot(x="slope_AA", y="slope_EA", data=samp_1_df, kind="reg", stat_func=r, color="red", 
                    scatter_kws={"s": 15, "alpha" : 0.6}, height=6)
        
        plt.tight_layout()
        plt.savefig(join(plot_path, '{}_Gene1_Scatter_histogram.pdf'.format(chrom)))
        plt.close()

        #sns.despine(offset=5, trim=True)


#################################################
#################################################
#################################################


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



#################
#################
## Test
#################
#################

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

chrom="chr1"
for chrom in chrom_list:
    print("\nProcessing : {} ...\n".format(chrom))

    """ parse afr cohorts """
    read_aa = pd.read_parquet(aa_dict.get(chrom), engine="fastparquet")
    read_aa.reset_index(inplace=True, drop=True)
    sorted_aa = read_aa.sort_values(["variant_id", "pval_nominal"], ascending=True).reset_index(drop=True)
    sorted_aa.drop_duplicates(["phenotype_id", "variant_id"], inplace=True, keep="first") #drop based on most sig pval

    aa_thresh = sorted_aa.groupby(["phenotype_id"]).apply(lambda X:  X[X["pval_nominal"] <= pval])
    aa_thresh.reset_index(drop=True, inplace=True)
    len(aa_thresh["phenotype_id"].unique())
    aa_thresh.drop_duplicates(["phenotype_id"], inplace=True, keep="first")
    aa_df = read_aa[read_aa["phenotype_id"] =="ENSG00000001461.16"]

    """ parse eur cohorts"""
    read_ea = pd.read_parquet(ea_dict.get(chrom), engine="fastparquet")
    read_ea.reset_index(inplace=True, drop=True)
    sorted_ea = read_ea.sort_values(["variant_id", "pval_nominal"], ascending=True).reset_index(drop=True)
    sorted_ea.drop_duplicates(["phenotype_id", "variant_id"], inplace=True, keep="first") #drop based on most sig pval

    ea_thresh = sorted_ea.groupby(["phenotype_id"]).apply(lambda X:  X[X["pval_nominal"] <= pval])
    ea_thresh.reset_index(drop=True, inplace=True)
    len(ea_thresh["phenotype_id"].unique())
    ea_thresh.drop_duplicates(["phenotype_id"], inplace=True, keep="first")
    aa_df = read_ea[read_ea["phenotype_id"] =="ENSG00000001461.16"]
    
    merged_df = pd.merge(aa_thresh, ea_thresh, on=["phenotype_id"], how = "outer", suffixes=["_AA", "_EA"], indicator=True)


############################
############################
############################
        # grouped = df.groupby(['user_nm', 'halves'])
        # my_lambda = lambda x, y: reduce(pd.merge(x, y, on=["phenotype_id"], how="outer"))
        # output = grouped.aggregate({'unique_ips': my_lambda,
        #                             'shifted_ips': my_lambda})
        # stats_scores = pd.merge(stats_beta, stats_zscore, on=["phenotype_id"], how="outer")
        # stats_filt_df = pd.merge(stats_scores, stats_kstest, on=["phenotype_id"], how="outer")

