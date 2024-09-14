#!/usr/bin/env python
import pandas as pd
import numpy as np
import pybedtools
from glob import glob
import os, re, pickle
from os.path import basename, join, splitext


# #split files into chrom files
# for chrom in {1..22};do 
#     echo "processing chr:${chrom}..."; 
#     grep "chr${chrom}_" ../sorted.all.AFR.Meta.eQTL.genPC.nominal.maf05.mvmeta.fe.txt > AA/chr${chrom}_montgoeQTL_AA.txt; 
# done


# for chrom in {1..22};do 
#     echo "processing chr:${chrom}..."; 
#     grep "^\b${chrom}\b" ../sorted.dist.hwe.af.EUR_META.eQTL.nominal.hg38a.copy.txt > EUR/chr${chrom}_montgoeQTL_EUR.txt; 
# done


"""supress warnings"""
import time, warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None

"""tissue-type"""
tissue = "EBV_LCL"

"""input files"""
aa_filepath = "/scratch16/abattle4/surya/datasets/prs_project/montgomery_lab/AFGR/prs_project/AFR"
ea_filepath = "/scratch16/abattle4/surya/datasets/prs_project/montgomery_lab/AFGR/prs_project/EUR"
gwas_stats = "/scratch16/abattle4/surya/datasets/prs_project/montgomery_lab/AFGR/prs_project/4080_irnt.sbp.gwas.imputed_v3.both_sexes.tsv"
#gwas_stats = "/Users/suryachhetri/datasets/prs_project/final_hg19/EBV_LCL/4080_irnt.sbp.gwas.imputed_v3.both_sexes.tsv"

output_path = "/scratch16/abattle4/surya/datasets/prs_project/montgomery_lab/AFGR/prs_project/output"
plot_path = join(output_path, "plots")

#requires `liftOver`from UCSC to be on the path and a `chainfile`
chainfile = "/scratch16/abattle4/surya/tools/hg38ToHg19.over.chain"

"""create plot output dir"""
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

    # tissue = "Cells_EBV-transformed_lymphocytes" #tissue = "Heart_Atrial"
    # chromrange = 22
    # chrom_X = False
    pval = 1e-5
    # pval = 1e-4

    chrom_list = []
    for chrom_no in range(chromrange):
        chrom_list.append("chr" + str(chrom_no + 1))
    if chrom_X:
        chrom_list.append("chrX")
        #chrom_list.append("chrY")

    # read chromwise for AFR and EUR:
    #aa_filelist = glob(join(aa_filepath, str(tissue) + "*.parquet"))
    #aa_dict = {basename(file).split(".")[4] :file for file in aa_filelist}
    aa_filelist = glob(join(aa_filepath, str("chr") + "*.txt"))
    aa_dict = {basename(file).split("_")[0] :file for file in aa_filelist}
    ea_filelist = glob(join(ea_filepath, str("chr") + "*.txt"))
    ea_dict = {basename(file).split("_")[0] :file for file in ea_filelist}

    #eQTL_filt_stats = []
    pop_shared_betas_ses = [] 
    aa_betas_ses = []
    ea_betas_ses = []
    shared_cis_genes = []
    sig_shared_cis_genes = []
    total_cis_genes = []

    #genecount = []
    #variantcount = []

    for chrom in chrom_list:
        #chrom="chr10"
        print("\nProcessing : {} ...\n".format(chrom))

        """ parse afr cohorts """
        col_names = ["Variant", "Transcript", "BETA", "SE", "ZSTAT", "P", "95PCT_CI_LOWER", \
                    "95PCT_CI_UPPER", "Q", "DF", "P_HET"]
        
        read_afa = pd.read_csv(aa_dict.get(chrom), sep="\t", names=col_names)
        select_cols = ["Transcript", "Variant", "P", "BETA", "SE"]
        read_aa = read_afa.loc[:, select_cols]
        read_aa.rename(columns={"Transcript" : "phenotype_id", "Variant" : "variant_id", "P": "pval_nominal", \
                                "BETA": "slope", "SE" :"slope_se"}, inplace=True)
        #read_aa["maf"] = 0.051 # dummy maf
        #aa_thresh = read_aa.groupby(["phenotype_id"]).apply(lambda X:  X[X["pval_nominal"] <= pval])
        #aa_df = read_aa[read_aa["phenotype_id"] =="ENSG00000001461.16"]


        """ parse eur cohorts"""
        col_names = ["chr", "snp_pos", "snp_pos2", "ref", "alt", "variant", "effect_af_eqtl", "p_hwe", "feature", 
                    "geneSymbol", "geneType", "dist_start", "dist_end", "pvalue", "beta", "se", "qstat", "p_het", "df"]

        read_eua = pd.read_csv(ea_dict.get(chrom), sep="\t", names=col_names)
        select_cols = ["feature", "variant", "pvalue", "beta", "se"]
        read_ea = read_eua.loc[:, select_cols]
        read_ea.rename(columns={"feature" : "phenotype_id", "variant" : "variant_id", "pvalue": "pval_nominal", \
                                "beta": "slope", "se" :"slope_se"}, inplace=True)
        read_ea["variant_id"] = "chr" + read_ea["variant_id"]
        
        #read_ea.reset_index(inplace=True, drop=True)
        #read_ea["variant_id"] = read_ea["variant_id"].str.replace("_b38", "")
        #ea_thresh = read_ea.groupby(["phenotype_id"]).apply(lambda X:  X[X["pval_nominal"] <= pval])
        #ea_df = read_ea[read_ea["phenotype_id"] =="ENSG00000001461.16"]
        
        #list of gene indexes with at least one variant having significant assoc:
        aa_genes = pd.Series(read_aa["phenotype_id"].unique(), name="phenotype_id")
        ea_genes = pd.Series(read_ea["phenotype_id"].unique(), name="phenotype_id")
        #aa_genes = aa_genes.str.replace("\\..*", "") #remove version numbers
        #ea_genes = ea_genes.str.replace("\\..*", "") #remove version numbers
        
        #union gives list of genes where at least one cohort has one variant per gene with sig assoc:
        ea_aa_union = pd.merge(aa_genes, ea_genes, how="outer")
        ea_aa_intersect = pd.merge(aa_genes, ea_genes, how="inner")
        total_cis_genes.append(ea_aa_union.shape[0])
        shared_cis_genes.append(ea_aa_intersect.shape[0])

        # select overlap/shared genes b/w cohorts
        aa_filt = pd.merge(read_aa, ea_aa_intersect, how="inner") 
        ea_filt = pd.merge(read_ea, ea_aa_intersect, how="inner")

        """find shared gene-variants between EA and AA pops on same gene"""
        """merge afr eur cohorts based on selected gene idx"""
        shared_vars = pd.merge(ea_filt, aa_filt, on=["phenotype_id", "variant_id"], how = "outer", suffixes=["_EA", "_AA"], indicator=True)
        shared_variants = shared_vars[shared_vars["_merge"]=="both"]        
        aa_variants = shared_vars[shared_vars["_merge"]=="right_only"]
        ea_variants = shared_vars[shared_vars["_merge"]=="left_only"]


        #filter based on sig pval_nominal 0.0001 of at least in one cohort
        shared_vars_thres = shared_variants.loc[(shared_variants['pval_nominal_AA'] <= pval) | (shared_variants['pval_nominal_EA'] <= pval)]
        aa_variants_thresh = aa_variants.loc[aa_variants['pval_nominal_AA'] <= pval]
        ea_variants_thresh = ea_variants.loc[ea_variants['pval_nominal_EA'] <= pval]
        
        #shared cis genes
        shared_genes = pd.Series(shared_variants["phenotype_id"].unique(), name="phenotype_id")
        
        #significant shared and unique cis genes 
        shared_threshpval_genes = pd.Series(shared_vars_thres["phenotype_id"].unique(), name="phenotype_id")
        aa_threshpval_genes = pd.Series(aa_variants_thresh["phenotype_id"].unique(), name="phenotype_id")
        ea_threshpval_genes = pd.Series(ea_variants_thresh["phenotype_id"].unique(), name="phenotype_id")

        #toal shared cis genes vs significant shared cis genes
        shared_genecount = len(shared_variants["phenotype_id"].unique())
        shared_threshpval_genecount = len(shared_vars_thres["phenotype_id"].unique())
        
        #find uniq genes to ea and aa not in shared thresholded genes list
        aa_filt_genes = set(aa_genes) - set(shared_threshpval_genes)
        ea_filt_genes = set(ea_genes) - set(shared_threshpval_genes)

        #Select significant shared genes in at least one cohort
        shared_pop_betas_ses = pd.merge(shared_variants, shared_threshpval_genes, how="inner")
        aa_pop_betas_ses = pd.merge(read_aa,  pd.DataFrame(aa_filt_genes).rename(columns={0: "phenotype_id"}), how="inner")
        ea_pop_betas_ses = pd.merge(read_ea,  pd.DataFrame(ea_filt_genes).rename(columns={0: "phenotype_id"}), how="inner")

        #append gene-variants pairs classified as eQTLgene in either one cohort 
        pop_shared_betas_ses.append(shared_pop_betas_ses)
        aa_betas_ses.append(aa_pop_betas_ses)
        ea_betas_ses.append(ea_pop_betas_ses)
        
        #append significant shared cis genes
        sig_shared_cis_genes.append(shared_threshpval_genecount)

        """basic summary stats chromwise"""
        print('''\n
        Total variants across cohorts ::
        European cohorts : {0}
        African cohorts : {1}'''.format(read_ea.shape[0], read_aa.shape[0]))

        print('''\n
        Total Shared cis gene-variants across cohorts ::
        Shared cis-genes : {0}
        Shared cis-variants : {1}'''.format(shared_genecount, shared_variants.shape[0]))
        
        print('''\n
        Significant gene count with at least 1 significant hit in either pop ::
        Sig Shared Genecount : {0}
        Sig Shared Cis-variants: {1}'''.format(shared_threshpval_genecount, pop_betas_ses.shape[0]))

    betas_ses_concat = pd.concat(pop_shared_betas_ses, ignore_index=True)
    #betas_ses_concat.to_csv(join(output_path, tissue + "_coloc_betas_ses_eqtl.hg38.shared.txt"), sep="\t", index=False, header=True )
    aa_ses_concat = pd.concat(aa_betas_ses, ignore_index=True)
    #aa_ses_concat.to_csv(join(output_path, tissue + "_coloc_betas_ses_eqtl.hg38.aa.txt"), sep="\t", index=False, header=True )
    ea_ses_concat = pd.concat(ea_betas_ses, ignore_index=True)
    #ea_ses_concat.to_csv(join(output_path, tissue + "_coloc_betas_ses_eqtl.hg38.ea.txt"), sep="\t", index=False, header=True )
    
    # summary stats whole genome
    print("\nTotal Sig shared eQTL counts :  {}".format(betas_ses_concat.shape[0]))
    print("\nTotal Sig shared cis_genes counts :  {}".format(sum(sig_shared_cis_genes)))
    print("\nTotal Shared cis_genes counts :  {}".format(sum(shared_cis_genes)))
    print("\nTotal cis_genes(unique + shared) :  {}".format(sum(total_cis_genes)))
    return(betas_ses_concat)


# liftover chrom wise:
def liftover(eqtl_dataframe):

    #column list to select cols for retention 
    #eqtl_dataframe = cis_eqtl_coloc_df.copy()
    select_cols =  eqtl_dataframe.columns.tolist()[1:-1] + [eqtl_dataframe.columns.tolist()[0]]
    new_df = eqtl_dataframe.loc[:, select_cols]
    new_df["phenotype_id"] = new_df["phenotype_id"].str.replace("\\..*", "")
    new_df.rename(columns={new_df.columns[0]: "newsnp"}, inplace=True)
    new_df["newsnp"] = new_df["newsnp"].str.replace("_", ":")
    new_df["idx"] =  np.arange(len(new_df))
    dfsubset = new_df["newsnp"].str.split(":", expand=True)

    #new_df["newsnp"] = new_df["newsnp"].str.replace("_b38", "").str.replace("_", ":")
    #df_test[["chrom", "end", "ref", "alt"]] = pd.new_df([ x.split(':') for x in list(df_test['newsnp']) ])
    #df_test[["chrom", "end", "ref", "alt"]] = df_test["newsnp"].str.split(":", return_type='frame')
    
    dfsubset.columns = ["chrom", "end", "ref", "alt"]
    new_df["ref"] = dfsubset["ref"]
    new_df["alt"] = dfsubset["alt"]
    dfsubset["start"] = dfsubset["end"].astype(int) - 1
    dfsubset["idx"] = new_df["idx"] #.str.cat(new_df[new_df.columns[1:].values], sep="_")
    select_cols = ["chrom", "start", "end", "idx"]
    subset_df = dfsubset.loc[:,select_cols]

    # perform lift chromosome wise to save memory kill
    chrom_list = []
    for chrom_no in range(22):
       chrom_list.append("chr" + str(chrom_no + 1))

    liftover_list = []
    for chrom in chrom_list:
        print("\nProcessing Liftover for {} ...".format(chrom))
        nsubset_df = subset_df.loc[subset_df["chrom"] == chrom]
        pybed_df = pybedtools.BedTool.from_dataframe(nsubset_df)
        hg19_pybed = pybed_df.liftover(chainfile)
        print("\nLiftover intermediate")
        hg19_df = pd.read_csv(hg19_pybed.fn, sep="\t", header=None)
        hg19_df.columns = ["chrom", "start", "end", "idx"]
        df_merged = pd.merge(hg19_df, new_df, on="idx", how="inner")
        print("Liftover dataframe merged")

        #df_merged[["chr", "pos", "allele"]] = df_merged["newsnp"].str.split(":", 2, expand=True)
        df_merged["hg19_snp"] = df_merged["chrom"] + ":" +  df_merged["end"].astype(str) + ":" +  \
                                df_merged["ref"] + ":" + df_merged["alt"]
        select_cols = ['hg19_snp', 'newsnp']+ new_df.columns[1:-3].values.tolist()
        hg19_final_df = df_merged.loc[:,select_cols]
        liftover_list.append(hg19_final_df)

    #combine the chromosomes
    hg19_final_df = pd.concat(liftover_list, ignore_index=True)
    #hg19_final_df.to_csv(join(output_path, tissue + "_coloc_betas_ses_eqtl.hg19.chromwise.txt"), 
    #                            sep="\t", index=False, header=True )
    print("\nLiftover completed")
    return(hg19_final_df)

#hg19_eqtl_coloc_df = liftover(cis_eqtl_coloc_df)


# liftover whole genome instead of chrom wise:
def liftover_wholegenome(eqtl_dataframe):

    #column list to select cols for retention 
    #eqtl_dataframe = cis_eqtl_coloc_df.copy()
    select_cols =  eqtl_dataframe.columns.tolist()[1:-1] + [eqtl_dataframe.columns.tolist()[0]]
    new_df = eqtl_dataframe.loc[:, select_cols]
    new_df["phenotype_id"] = new_df["phenotype_id"].str.replace("\\..*", "")
    new_df.rename(columns={new_df.columns[0]: "newsnp"}, inplace=True)
    new_df["newsnp"] = new_df["newsnp"].str.replace("_", ":")
    new_df["idx"] =  np.arange(len(new_df))
    dfsubset = new_df["newsnp"].str.split(":", expand=True)

    #new_df["newsnp"] = new_df["newsnp"].str.replace("_b38", "").str.replace("_", ":")
    #df_test[["chrom", "end", "ref", "alt"]] = pd.new_df([ x.split(':') for x in list(df_test['newsnp']) ])
    #df_test[["chrom", "end", "ref", "alt"]] = df_test["newsnp"].str.split(":", return_type='frame')
    
    dfsubset.columns = ["chrom", "end", "ref", "alt"]
    new_df["ref"] = dfsubset["ref"]
    new_df["alt"] = dfsubset["alt"]
    dfsubset["start"] = dfsubset["end"].astype(int) - 1
    dfsubset["idx"] = new_df["idx"] #.str.cat(new_df[new_df.columns[1:].values], sep="_")
    select_cols = ["chrom", "start", "end", "idx"]
    subset_df = dfsubset.loc[:,select_cols]

    pybed_df = pybedtools.BedTool.from_dataframe(subset_df)
    hg19_pybed = pybed_df.liftover(chainfile)
    print("\nLiftover intermediate")
    hg19_df = pd.read_csv(hg19_pybed.fn, sep="\t", header=None)
    hg19_df.columns = ["chrom", "start", "end", "idx"]
    df_merged = pd.merge(hg19_df, new_df, on="idx", how="inner")
    print("Liftover dataframe merged")

    #df_merged[["chr", "pos", "allele"]] = df_merged["newsnp"].str.split(":", 2, expand=True)
    df_merged["hg19_snp"] = df_merged["chrom"] + ":" +  df_merged["end"].astype(str) + ":" +  \
                            df_merged["ref"] + ":" + df_merged["alt"]
    select_cols = ['hg19_snp', 'newsnp']+ new_df.columns[1:-3].values.tolist()
    hg19_final_df = df_merged.loc[:,select_cols]
    liftover_list.append(hg19_final_df)

    print("\nLiftover completed")
    return(hg19_final_df)


# Common GWAS with hg19 assembly:
def eQTL_gwas_overlap(eqtl_filt_loci, gwas_summary_stats_file):
    gwas_df = pd.read_csv(gwas_summary_stats_file, sep="\t")
    select_cols = ["variant", "beta", "se", "pval", "minor_AF"]
    hg19_gwas = gwas_df.loc[:,select_cols]
    hg19_gwas["variant"] =  "chr" + hg19_gwas["variant"].astype(str)

    # merge eQTL GWAS:
    merged_final_df = pd.merge(eqtl_filt_loci, hg19_gwas, left_on=["hg19_snp"], right_on=["variant"], how="inner")
    #merged_final_df = pd.merge( hg19_eqtl_coloc_df, hg19_gwas, left_on=["hg19_snp"], right_on=["variant"], how="inner")
    
    # compute percent overlap:
    unique_merged_df = merged_final_df.drop_duplicates(subset=["hg19_snp"])
    unique_eqtl_df = eqtl_filt_loci.drop_duplicates(subset=["hg19_snp"])
    total_eQTL_GWAS_count = unique_merged_df.shape[0]
    total_unique_eQTL_count = unique_eqtl_df.shape[0]
    percent = round((total_eQTL_GWAS_count/float(total_unique_eQTL_count))*100, 2)
    
    # basic summary stats
    print("Total eQTL-GWAS counts :  {}".format(total_eQTL_GWAS_count))
    print("Total unique-eQTL counts :  {}".format(total_unique_eQTL_count))
    print("eQTL Percent Overlap with GWAS  :  {}%%".format(percent)) 
    return(merged_final_df)

#all_eQTL_gwas_loci = eQTL_gwas_colocformat(hg19_eqtl_coloc_df, gwas_stats)

# For baseline study
def unique_eQTL_gwas_overlap(eqtl_filt_loci, gwas_summary_stats_file):
    gwas_df = pd.read_csv(gwas_summary_stats_file, sep="\t")
    gwas_df = pd.read_csv(gwas_stats, sep="\t")
    select_cols = ["variant", "beta", "se", "pval", "minor_AF"]
    hg19_gwas = gwas_df.loc[:,select_cols]
    hg19_gwas["variant"] =  "chr" + hg19_gwas["variant"].astype(str)

    # merge eQTL GWAS:
    merged_final_df = pd.merge(eqtl_filt_loci, hg19_gwas, left_on=["hg19_snp"], right_on=["variant"], how="inner")
    
    # compute percent overlap:
    unique_merged_df = merged_final_df.drop_duplicates(subset=["hg19_snp"])
    unique_eqtl_df = eqtl_filt_loci.drop_duplicates(subset=["hg19_snp"])
    total_eQTL_GWAS_count = merged_final_df.shape[0]
    total_unique_eQTL_count = unique_eqtl_df.shape[0]
    percent = round((total_eQTL_GWAS_count/float(total_unique_eQTL_count))*100, 2)
    
    #basic summary stats
    print("Total eQTL-GWAS counts :  {}".format(total_eQTL_GWAS_count))
    print("Total unique-eQTL counts :  {}".format(total_unique_eQTL_count))
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
    #tissue = "Adipose_Subcutaneous" #tissue = "Heart_Atrial"
    tissue = "EBV_LCL" #tissue = "Heart_Atrial"
    chromrange = 22
    chrom_X = False
    pval = 1e-5

    # Shared eQTLs and storage in compressed formats: read/write with pickle
    cis_eqtl_coloc_df = pop_shared_variants_colocformat(aa_filepath, ea_filepath, chromrange, pval, tissue, chrom_X=False)    
    cis_eqtl_coloc_df.to_csv(join(output_path, tissue + "_coloc_betas_ses_eqtl.hg38.txt"), sep="\t", index=False, header=True )

    # Dframe based genomelift of hg38 genome to hg19:
    hg19_eqtl_coloc_df = liftover(cis_eqtl_coloc_df)
    hg19_eqtl_coloc_df.to_csv(join(output_path, tissue + "_coloc_betas_ses_eqtl.hg19.txt"), sep="\t", index=False, header=True )

    # # eQTL loci associated to GWAS loci:
    eQTL_gwas_loci  = eQTL_gwas_overlap(hg19_eqtl_coloc_df, gwas_stats)
    eQTL_gwas_loci.drop(["newsnp", "variant"], inplace=True, axis=1)
    eQTL_gwas_loci.to_csv(join(output_path, tissue + "_coloc_betas_ses.eqtl_gwas.hg19.txt"), sep="\t", index=False, header=True )

    end_time = time.time()
    total_time = end_time - start_time
    print("Time for analysis : {}".format(total_time))
    print("Task completed!")



###############################
###############################
###############################
#!/usr/bin/env Rscript

suppressMessages({

    #load libraries:
    library(data.table)
    library(tidyverse)
    library(coloc)
    library(optparse, quietly=TRUE)

})

#timeit
"https://www.cell.com/ajhg/pdf/S0002-9297(21)00271-8.pdf"
start_time <- Sys.time()

#input filenames:
snpcount <- 10000
tissue <- "EBV_LCL"
#tissue <- "Cells_EBV-transformed_lymphocytes"
output_path <- "/Users/suryachhetri/datasets/prs_project/montgomery_lab/AFGR/prs_project/output"
#output_path <- "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL"

#load files:
colocfile <- file.path(output_path, paste0(tissue, "_coloc_betas_ses_eqtl.hg19.txt"))
#colocfile <- file.path(output_path, paste0(tissue, "_coloc_betas_ses_eqtl.txt"))
coloc_df <- fread(colocfile, sep="\t")
#ncoloc_df <- fread(colocfile, sep="\t")

#create unique genelist:
#coloc_df <- coloc_df %>% rename(phenotype_id = gene)
#coloc_df <- coloc_df %>% rename(gene = phenotype_id)

genevect <- coloc_df %>% distinct(phenotype_id) %>% unlist(use.name=F)
#genevect <- coloc_df %>% distinct(gene) %>% unlist(use.name=F)

#distinct gene count:
coloc_df %>% distinct(phenotype_id) %>% nrow
#coloc_df %>% distinct(gene) %>% nrow

#create empty list for rbind:
gene_coloc <- list()
snp_coloc <- list()

idx_vect <- c(1:length(genevect))
#masterlist <- list(genevect, idx_vect)

for(i in idx_vect) {
    
    cat("\n")
    print(paste("Processing gene ..:", i))
    cat("\n")

    #gene <- masterlist[[1]][i]
    gene <- genevect[i] 
    print(paste0("genename: ", gene));
    gene_df <- coloc_df %>% filter(phenotype_id == gene)
    gene_df <- gene_df %>% drop_na
    #snpid <- list(snpid=gene_df$hg19_snp, snp=paste0("SNP.", seq(1:nrow(gene_df))))
    snpid <- data.frame(snpid=gene_df$hg19_snp, snp=paste0("SNP.", seq(1:nrow(gene_df))))
    
    betas <- gene_df %>% 
        select(slope_EA, slope_AA) %>%
        rename(EA = slope_EA, AA = slope_AA)

    ses <- gene_df %>% 
        select(slope_se_EA, slope_se_AA) %>%
        rename(EA = slope_se_EA, AA = slope_se_AA)

    # maf <- gene_df %>% 
    #     select(maf_EA, maf_AA) %>%
    #     rename(EA = maf_EA, AA = maf_AA)

    #colocalization:
    # res1 <- coloc.abf(dataset1=list(beta=betas$EA, varbeta=ses$EA^2, N=nrow(betas),sdY=1,type="quant"),
    #                     dataset2=list(beta=betas$AA, varbeta=ses$AA^2, N=nrow(betas),sdY=1,type="quant"),
    #                     MAF=maf$EA, p1=1e-04, p2=1e-04, p12=5e-05) # p12=1e-05

    res <- coloc.abf(dataset1=list(beta=betas$EA, varbeta=ses$EA^2, sdY=1, type="quant"),
                        dataset2=list(beta=betas$AA, varbeta=ses$AA^2, sdY=1, type="quant")) # p12=1e-05

    #colocalization probablility:
    local_list <- list()
    local_list[[gene]] <- res$summary[c(2:6)]

    #create summary dataframe
    res_sum_df <- data.frame(local_list) %>% t() %>% data.frame %>% rownames_to_column("gene")
    
    #merge snpid dataframe with summary dataframe
    snp_sum_df <- data.frame(snpid, res_sum_df)

    #extract variant wise pp
    snp_pp_df <- res$result %>% select(snp, SNP.PP.H4)

    #combine gene pp and variant wise pp
    concat_df <- left_join(snp_pp_df, snp_sum_df, by = c("snp" = "snp"))

    #select and clean up columns
    new_df <- concat_df %>% select(c(-snp, -gene))
    final_df <- new_df %>% select(PP.H0.abf:PP.H4.abf,snpid,SNP.PP.H4)
    
    #append coloc result and snp pp info list
    gene_coloc[[gene]] <- res$summary[c(2:6)]
    snp_coloc[[gene]] <- final_df

    #compute posterior probability that the variant is causal
    causal_ppH4_df <- final_df %>% mutate(SNP.PP.H4.Causal = PP.H4.abf*SNP.PP.H4, phenotype_id=gene)

    #write snp dataframe at gene level to concat from terminal later
    outdir <- file.path(output_path, "Genewise_colocsnpDF")
    dir.create(outdir)
    write.table(causal_ppH4_df, file.path(outdir, paste0(i, "_", gene, ".txt")), quote=F, sep="\t", col.names=F, row.names=F)


}


#concat dataframe
gene_coloc_df <- data.frame(gene_coloc) %>% t() %>% data.frame %>% rownames_to_column("gene") %>% data.table
#snp_pp_df <- bind_rows(snp_coloc, .id="phenotype_id") %>% data.table

#label filenames
output_genelevel <- file.path(output_path, paste0(tissue, "_coloc_result_Genelevel.txt"))
output_snplevel <- file.path(output_path, paste0(tissue, "_coloc_result_Snplevel.txt"))

#instead of bind_rows do concat of files
outdir <- file.path(output_path, "Genewise_colocsnpDF")
setwd(outdir)
system(paste0("cat ", "*.txt", " >", " ../merged_test.txt"))

#write file
write.table(gene_coloc_df, output_genelevel, quote=F, sep="\t", col.names=T, row.names=F)
#write.table(snp_pp_df, output_snplevel, quote=F, sep="\t", col.names=T, row.names=F)
#genecoloc_df <- bind_rows(!!!gene_coloc, .id="phenotype_id") %>% as_tibble()

#load file
gene_coloc_df <- fread(output_genelevel, sep="\t")

#save R object to a file
# output_genelevel <- file.path(output_path, paste0(tissue, "_coloc_result_Genelevel.txt"))
# output_snplevel <- file.path(output_path, paste0(tissue, "_coloc_result_Snplevel.txt"))

saveRDS(gene_coloc_df, file=file.path(output_path, "gene_coloc_df.RData"))
coloc_df <- readRDS(file.path(output_path, "gene_coloc_df.RData"))


saveRDS(snp_coloc, file=file.path(output_path, "snp_coloc_list.RData"))
snp_coloc <- readRDS(file.path(output_path, "snp_coloc_list.RData"))

#write dataframes to file
write.table(coloc_df, output_genelevel, quote=F, sep="\t", col.names=T, row.names=F)
write.table(snp_coloc, output_snplevel, quote=F, sep="\t", col.names=T, row.names=F)


###########################
# Complex heatmap Library
###########################
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(data.table)
library(dplyr)


##########################################################################
##Coloc profile with each coloc category:
##########################################################################

#input filenames:
#tissue <- "Cells_EBV-transformed_lymphocytes"
tissue <- "EBV_LCL"
#output_path <- "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL"
output_path <- "/Users/suryachhetri/datasets/prs_project/montgomery_lab/AFGR/prs_project/output"
input_file <- file.path(output_path, paste0(tissue, "_coloc_result_Genelevel.txt"))

#load input file:
read_df <- fread(input_file, header=T) %>% data.frame
data <- read_df[,2:ncol(read_df)]
mat_data <- data %>% as.matrix

## Rowwise scaling - interested in meth variance across cis-reg regions for any samples
#mat_data <- t(apply(mat_data, 1, scale))
mat_data[is.na(mat_data)] <- 0

## Naming of cols and rows
colnames(mat_data) <- c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")
rownames(mat_data) <- read_df[,1]
rownames(mat_data)
colnames(mat_data)

output_file <- file.path(output_path, paste0(tissue, "_coloc_result_heatmap.pdf"))
pdf(output_file)

#custom_col=colorRamp2(c(0,0.5,1), c("blue", "red", "white"))
#custom_col = colorRamp2(c(min(mat_data), 0, max(mat_data)), c("blue", "white", "red"))
ht7  <- Heatmap(mat_data, name="Posterior Prob", 
    column_title="Colocalization Posterior",
    row_title="Ensemble Genes (EA vs AA)",
    na_col = "orange",
    #col=custom_col,
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 0),
    column_names_gp = gpar(fontsize = 10),
    cluster_columns = FALSE
    #km=7
    # show_column_dend = FALSE #(Hiding column clusters)
    #show_row_dend = FALSE #(Hiding column clusters)
    ## For heatmap and annotation legends (add new-legends)
    # heatmap_legend_param = list(
    # title_gp = gpar(fontsize = 10), 
    # labels_gp = gpar(fontsize = 6),
    # legend_direction = "horizontal",
    # legend_width = unit(3, "cm"), title_position = "topcenter"
    # )
) 

ht7 <- draw(ht7)
dev.off()



##########################
#threshold cutoff PP 0.75
##########################

tissue <- "EBV_LCL"
output_path <- "/Users/suryachhetri/datasets/prs_project/montgomery_lab/AFGR/prs_project/output"

outdir = file.path(output_path, "PPH_coloc_files")
dir.create(outdir)

#write geneset files for annotation prep and ldsc runs:
PP.H0 = file.path(outdir, "PP.H0.coloc.abf.txt")
PP.H1 = file.path(outdir, "PP.H1.coloc.abf.txt")
PP.H2 = file.path(outdir, "PP.H2.coloc.abf.txt")
PP.H3 = file.path(outdir, "PP.H3.coloc.abf.txt")
PP.H4 = file.path(outdir, "PP.H4.coloc.abf.txt")
fullENSG = file.path(outdir, "PP.fullENSG.coloc.abf.txt")
#Null = file.path(outdir, "Null.coloc.abf.txt")

#filter coloc evidence and write files
read_df %>% filter(PP.H0.abf>0.75) %>%
    write.table(., PP.H0, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H1.abf>0.75) %>%
    write.table(., PP.H1, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H2.abf>0.75) %>%
    write.table(., PP.H2, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H3.abf>0.75) %>%   
    write.table(., PP.H3, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H4.abf>0.75) %>%
    write.table(., PP.H4, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>%
    write.table(., fullENSG, quote=F, sep="\t", row.names=F, col.names=T)

# #generate random null gene file:
# set.seed(1)
# read_df1 %>% sample_n(.,size) %>%
#     write.table(., Null, quote=F, row.names=F, col.names=T)


###################################################
# continuum/range of PPH4 colocalization evidence:
###################################################


outdir = file.path(output_path, "PPH4_coloc_files")
dir.create(outdir)

PP.H4.05 = file.path(outdir, "PP.H4.05.coloc.abf.txt")
PP.H4.1 = file.path(outdir, "PP.H4.1.coloc.abf.txt")
PP.H4.2 = file.path(outdir, "PP.H4.2.coloc.abf.txt")
PP.H4.3 = file.path(outdir, "PP.H4.3.coloc.abf.txt")
PP.H4.4 = file.path(outdir, "PP.H4.4.coloc.abf.txt")
PP.H4.5 = file.path(outdir, "PP.H4.5.coloc.abf.txt")
PP.H4.6 = file.path(outdir, "PP.H4.6.coloc.abf.txt")
PP.H4.7 = file.path(outdir, "PP.H4.7.coloc.abf.txt")
PP.H4.75 = file.path(outdir, "PP.H4.75.coloc.abf.txt")
PP.H4.8 = file.path(outdir, "PP.H4.8.coloc.abf.txt")
PP.H4.9 = file.path(outdir, "PP.H4.9.coloc.abf.txt")
PP.H4.95 = file.path(outdir, "PP.H4.95.coloc.abf.txt")
PP.H4.99 = file.path(outdir, "PP.H4.99.coloc.abf.txt")


#filter coloc evidence and write files
read_df %>% filter(PP.H4.abf>0.05) %>%
    write.table(., PP.H4.05, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H4.abf>0.1) %>%
    write.table(., PP.H4.1, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H4.abf>0.2) %>%
    write.table(., PP.H4.2, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H4.abf>0.3) %>%
    write.table(., PP.H4.3, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H4.abf>0.4) %>%
    write.table(., PP.H4.4, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H4.abf>0.5) %>%
    write.table(., PP.H4.5, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H4.abf>0.6) %>%
    write.table(., PP.H4.6, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H4.abf>0.7) %>%
    write.table(., PP.H4.7, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H4.abf>0.75) %>%
    write.table(., PP.H4.75, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H4.abf>0.8) %>%
    write.table(., PP.H4.8, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H4.abf>0.9) %>%
    write.table(., PP.H4.9, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H4.abf>0.95) %>%
    write.table(., PP.H4.95, quote=F, sep="\t", row.names=F, col.names=T)

read_df %>% filter(PP.H4.abf>0.99) %>%
    write.table(., PP.H4.99, quote=F, sep="\t", row.names=F, col.names=T)



##############################################

#Fetch TOP-eQTL snps from each coloc category

##############################################

library(data.table)
library(tidyverse)
# library(dplyr)
# library(purrr)
# library(tidyr)

#print full path:
#printf '%s\n' $PWD/PP*.75*txt
#system('printf "%s\n" $PWD/PP*.75*txt > PP_coloc_filelist.txt')

topeQTL_snps <- function(coloc_info_df, ppH_genepath, output_dir){

  #coloc_info <- fread(coloc_file, sep="\t")
  #coloc_info <- fread("/Users/suryachhetri/datasets/prs_project/montgomery_lab/AFGR/prs_project/output/EBV_LCL_coloc_betas_ses_eqtl.hg19.txt", sep="\t")
  # coloc_info <- coloc_info %>% rename(gene = phenotype_id, variant_id = hg19_snp)
  ppH4_genes <- fread(ppH_genepath, sep="\t")
  #ppH4_genes <- fread("/Users/suryachhetri/datasets/prs_project/montgomery_lab/AFGR/prs_project/output/PPH_coloc_files/PP.fullENSG.coloc.abf.txt", sep="\t")
  ppH4_genes <- ppH4_genes %>% select("gene")
  #names(ppH4_genes) <- "gene"

  df_info <- inner_join(ppH4_genes, coloc_info)
  cat("\n")

  df_new <- df_info %>% 
    mutate(abs_ztest_EA=abs(slope_EA/slope_se_EA),
          abs_ztest_AA=abs(slope_AA/slope_se_AA),
          pvalEA=2*pnorm(abs_ztest_EA, lower.tail = FALSE), #factor 2 given 2 tail distbn
          pvalAA=2*pnorm(abs_ztest_AA, lower.tail = FALSE)  #factor 2 given 2 tail distbn
          )

  df_EA_topeqtl <- df_new %>% 
    group_by(gene) %>%
    arrange(desc(abs_ztest_EA)) %>%
    slice(1) %>%
    mutate(metapval=pchisq(sum(-2*log(c(pvalEA,pvalAA))),df = 2,lower.tail = FALSE),
      pop = "EA",
      snp_id=paste(gene,variant_id,pop, sep="_")
    )

  df_AA_topeqtl <- df_new %>% 
    group_by(gene) %>%
    arrange(desc(abs_ztest_AA)) %>%
    slice(1) %>%
    mutate(metapval=pchisq(sum(-2*log(c(pvalEA,pvalAA))),df = 2,lower.tail = FALSE),
      pop = "AA",
      snp_id=paste(gene,variant_id,pop, sep="_")
    )


  top_eQTL <- function(gene_id){
    df_EA <- df_EA_topeqtl %>% filter(gene==gene_id) %>% data.frame
    df_AA <- df_AA_topeqtl %>% filter(gene==gene_id) %>% data.frame
    snp_id <- c(df_EA$snp_id, df_AA$snp_id)
    metapval <- c(df_EA$metapval,df_AA$metapval)
    min_idx <- which(metapval==min(metapval))
    min_metapval <- metapval[min_idx]
    topeqtl <- snp_id[min_idx]
    return(paste(topeqtl, min_metapval, sep="_"))
  }

  #sanity check with purrr:map
  #print(map(c("ENSG00000077935","ENSG00000188976"), top_eQTL))
  #print(map(c("ENSG00000155542","ENSG00000165171"), top_eQTL))

  #apply on all gene lists (eqiv to lapply)
  map_list <- map(ppH4_genes$gene, top_eQTL)

  #generate snp dataframe from purrr:map list
  snp_alldf <- data.table(unlist(map_list))
  names(snp_alldf) <- "snp_id"

  #split snp_id to dataframe columns
  snp_split_df <- snp_alldf %>% 
    separate(snp_id, sep="_",
    into = c("gene", "variant_id", "pop", "pval")
    )

  #deduplicate the snp_id based on rs_id and gene
  snp_finaldf <- snp_split_df %>% distinct(variant_id, gene, .keep_all = TRUE)

  #furthermore deduplicate genes for separate tie breaker pvalue - unique eQTL 
  #print(map(c("ENSG00000155542","ENSG00000165171"), top_eQTL))
  snp_finaldf <- snp_finaldf %>% distinct(gene, .keep_all = TRUE)
  snp_list <- data.table(snp_finaldf$variant_id)
  snp_finallist <- snp_list %>% distinct(V1) #select unique rsid

  #write dataframe and snp list to file
  outlist <- unlist(strsplit(ppH_genepath, "/|.txt"))
  outfile <- outlist[length(outlist)]
  outfile1 <- paste0(outfile, ".topsnpDF_forLDSC.txt")
  outfile2 <- paste0(outfile,".topsnpList_forLDSC.txt")

  snpDF <- file.path(output_dir, "snpDF")
  if (!dir.exists(snpDF)){
    dir.create(snpDF)
    }

  fwrite(snp_finaldf, file.path(snpDF, outfile1), sep="\t", row.names=F, col.names=T)
  fwrite(snp_finallist, file.path(output_dir, outfile2), sep="\t", row.names=F, col.names=T)

}


main_func <- function(){
  #print full path:
  # system('printf "%s\n" $PWD/PP.H*.75.txt > PP_coloc_filelist.txt')
  # system('printf "%s\n" $PWD/PP.H*.txt > PP_coloc_fullfilelist.txt')

  # ls /Users/suryachhetri/datasets/prs_project/montgomery_lab/AFGR/prs_project/output/PPH_coloc_files > ../PPH_coloc_filelist.txt
  # ls /Users/suryachhetri/datasets/prs_project/montgomery_lab/AFGR/prs_project/output/PPH4_evidence_range > ../PPH4_coloc_filelist.txt
  
  #work_dir <- "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL"
  work_dir <- "/Users/suryachhetri/datasets/prs_project/montgomery_lab/AFGR/prs_project/output"
  coloc_filelist <- file.path(work_dir, "PPH_coloc_filelist.txt")
  coloc_filelist <- file.path(work_dir, "PPH4_coloc_filelist.txt")
  coloc_file <- file.path(work_dir, "EBV_LCL_coloc_betas_ses_eqtl.hg19.txt")

  output_dir <- file.path(work_dir, "PPH_coloc_topsnps")
  output_dir <- file.path(work_dir, "PPH4_coloc_topsnps")
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }

  #output_dir1 <- file.path(work_dir, "coloc_topsnps", "snpDF")
  # snpDF <- file.path(output_dir, "snpDF")
  # if (!dir.exists(snpDF)){
  #   dir.create(snpDF)
  #   }

  #read original coloc snp file
  coloc_info <- fread(coloc_file, sep="\t")
  coloc_info <- coloc_info %>% rename(gene = phenotype_id, variant_id = hg19_snp)

  filelist <- readLines(coloc_filelist)
  #filelist <- "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/PP.H4.75.coloc.abf.txt"
  for (ppH_genepath in filelist){
    print(paste("Processing Coloc-file:", ppH_genepath))
    topeQTL_snps(coloc_info, ppH_genepath, output_dir)
    }
}

main_func()

#merge coloc results to PPH3|PPH4 SNP_DF annotation 
# 1. PPH4
# 2. PPH3
# 3. PPH4/(PPH3 + PPH4)
# 4. PPH3/(PPH3 + PPH4)

#tissue <- "Cells_EBV-transformed_lymphocytes"
tissue <- "EBV_LCL"
#output_path <- "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL"
output_path <- "/Users/suryachhetri/datasets/prs_project/montgomery_lab/AFGR/prs_project/output"
work_dir <- output_path

#final coloc result file with genes by PPH posterior probability
#coloc_result_file <- file.path(output_path, paste0(tissue, "_coloc_result.txt"))
coloc_result_file <- file.path(output_path, paste0(tissue, "_coloc_result_Genelevel.txt"))

################
#Extract prior prob for PPH3
PPH3_file <- file.path(work_dir,"PPH_coloc_topsnps/snpDF/PP.H3.coloc.abf.topsnpDF_forLDSC.txt")
PPH4_file <- file.path(work_dir,"PPH_coloc_topsnps/snpDF/PP.H4.coloc.abf.topsnpDF_forLDSC.txt")

#load coloc results due to odd col format with "number" -- setting names of cols
coloc_df <- fread(coloc_result_file, sep="\t")
coloc_df <- coloc_df %>% rename(phenotype_id = gene)
#names(coloc_df) <- c("phenotype_id", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
#coloc_df <- coloc_df %>% select(-num)

#load snp dataframe containing genes info
snp_df_pph3 <- fread(PPH3_file, sep="\t")
snp_df_pph3 <- snp_df_pph3 %>% rename(phenotype_id = gene)
#names(snp_df_pph3) <- c("rs_id", "phenotype_id", "snp_id", "cohort", "pval")

#merge with final coloc df
merged_df_pph3 <- left_join(snp_df_pph3, coloc_df, by=c("phenotype_id"))

merged_final_pph3 <- merged_df_pph3 %>% 
    separate(variant_id, into=c("chr", "bp", "ref", "alt"), "[:]") %>% 
    mutate(snp_id=paste0(chr,":",bp)) %>% 
    select(c(snp_id, contains("PP.H")), phenotype_id) %>% 
    #mutate(PP.H3.ratio.abf=PP.H3.abf/(PP.H3.abf+PP.H4.abf)) %>%
    #mutate(PP.H4.ratio.abf=PP.H4.abf/(PP.H3.abf+PP.H4.abf))
    select(snp_id, PP.H3.abf, phenotype_id) %>%
    mutate(snp_id = str_replace(snp_id, "chr", "")) #format sync with bim file


merged_final_pph3_div_h3h4 <- merged_df_pph3 %>% 
    separate(variant_id, into=c("chr", "bp", "ref", "alt"), "[:]") %>% 
    mutate(snp_id=paste0(chr,":",bp)) %>% 
    select(c(snp_id, contains("PP.H")), phenotype_id) %>% 
    mutate(PP.H3.ratio.abf=PP.H3.abf/(PP.H3.abf+PP.H4.abf)) %>%
    select(snp_id, PP.H3.ratio.abf, phenotype_id) %>%
    mutate(snp_id = str_replace(snp_id, "chr", "")) #format sync with bim file


fwrite(merged_final_pph3, file.path(output_path, "PPH3_EBV_coloc_prior.txt"), sep="\t", col.names=T, row.names=F)
fwrite(merged_final_pph3_div_h3h4, file.path(output_path, "PPH3ratio_EBV_coloc_prior.txt"), sep="\t", col.names=T, row.names=F)

# #create control file
# merged_final_pph3$PP.H3.control <- 1
# pph3_control <- merged_final_pph3 %>% select(snp_id, PP.H3.control, phenotype_id)
# fwrite(pph3_control, file.path(output_path, "PPH3control_EBV_coloc_prior.txt"), sep="\t", col.names=T, row.names=F)

# #identify duplicated df only
# duplicate_df <- merged_final_pph3 %>% group_by(phenotype_id) %>% filter(n()>1) %>% data.frame


################
#Repeat extraction of prior prob for PPH4
PPH4_file <- file.path(work_dir,"PPH_coloc_topsnps/snpDF/PP.H4.coloc.abf.topsnpDF_forLDSC.txt")
#PPH4_file <- file.path(work_dir,"PPH_coloc_topsnps/snpDF/PP.H4.coloc.abf.topsnpDF_forLDSC.txt")

#load coloc results due to odd col format with "number" -- setting names of cols
coloc_df <- fread(coloc_result_file, sep="\t")
coloc_df <- coloc_df %>% rename(phenotype_id = gene)
# names(coloc_df) <- c("num", "phenotype_id", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
# coloc_df <- coloc_df %>% select(-num)

#load snp dataframe containing genes info
snp_df_pph4 <- fread(PPH4_file, sep="\t")
snp_df_pph4 <- snp_df_pph4 %>% rename(phenotype_id = gene)
#names(snp_df_pph4) <- c("rs_id", "phenotype_id", "snp_id", "cohort", "pval")

#merge with final coloc df
merged_df_pph4 <- left_join(snp_df_pph4, coloc_df, by=c("phenotype_id"))

merged_final_pph4 <- merged_df_pph4 %>% 
    separate(variant_id, into=c("chr", "bp", "ref", "alt"), "[:]") %>% 
    mutate(snp_id=paste0(chr,":",bp)) %>% 
    select(c(snp_id, contains("PP.H")), phenotype_id) %>% 
    #mutate(PP.H3.ratio.abf=PP.H3.abf/(PP.H3.abf+PP.H4.abf)) %>% 
    #mutate(PP.H4.ratio.abf=PP.H4.abf/(PP.H3.abf+PP.H4.abf)) %>%
    select(snp_id, PP.H4.abf, phenotype_id) %>%
    mutate(snp_id = str_replace(snp_id, "chr", "")) #format sync with bim file


merged_final_pph4_div_h3h4 <- merged_df_pph4 %>% 
    separate(variant_id, into=c("chr", "bp", "ref", "alt"), "[:]") %>% 
    mutate(snp_id=paste0(chr,":",bp)) %>% 
    select(c(snp_id, contains("PP.H")), phenotype_id) %>% 
    mutate(PP.H4.ratio.abf=PP.H4.abf/(PP.H3.abf+PP.H4.abf)) %>%
    select(snp_id, PP.H4.ratio.abf, phenotype_id) %>%
    mutate(snp_id = str_replace(snp_id, "chr", "")) #format sync with bim file

fwrite(merged_final_pph4, file.path(output_path, "PPH4_EBV_coloc_prior.txt"), sep="\t", col.names=T, row.names=F)
fwrite(merged_final_pph4_div_h3h4, file.path(output_path, "PPH4ratio_EBV_coloc_prior.txt"), sep="\t", col.names=T, row.names=F)

# #create control file
# merged_final_pph4$PP.H4.control <- 1
# pph4_control <- merged_final_pph4 %>% select(snp_id, PP.H4.control)
# fwrite(pph4_control, file.path(output_path, "PPH4control_EBV_coloc_prior.txt"), sep="\t", col.names=T, row.names=F)

# #duplicated df only
# duplicate_df <- merged_final_pph4 %>% group_by(phenotype_id) %>% filter(n()>1) %>% data.frame


#################################
##prepare gene based annot file for ldsc run:
#################################
#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool
import gzip
import os
from os.path import join

#input filenames:
tissue = "Cells_EBV-transformed_lymphocytes"
output_path = "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL"
#input_file = join(output_path, tissue + "_coloc_result.txt")

annotpath = join(output_path, "ldsc_annot")
#annotpath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/my_annot"
#annotpath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/my_annot"

# prefix="EBV_LCL_Binary_100kb"
# prefix="EBV_LCL_Binary_100kb_minusnull"
# prefix="EBV_LCL_Binary_100kb_fullENSGnull"

# prefix="EBV_LCL_Binary_500kb_minusnull"
# prefix="EBV_LCL_Binary_500kb"
# prefix="EBV_LCL_Binary_500kb_fullENSGnull"

# prefix="EBV_LCL_Binary_100kb_sequential"
# prefix="EBV_LCL_Binary_100kb_sequentialfull"
# prefix="EBV_LCL_Binary_100kb_colocthresh08"
# prefix="EBV_LCL_Binary_100kb_colocthresh07"

prefix="EBV_LCL_Binary_100kb_colocthresh075"
# prefix="EBV_LCL_Binary_50kb_colocthresh075"
# prefix="EBV_LCL_Binary_10kb_colocthresh075"
# prefix="EBV_LCL_Binary_500kb_colocthresh075"
# prefix="EBV_LCL_Binary_1MB_colocthresh075"
# prefix="EBV_LCL_Binary_2kb_colocthresh075"
# prefix="EBV_LCL_Binary_5kb_colocthresh075"
# prefix="EBV_LCL_Binary_25kb_colocthresh075"
# prefix="EBV_LCL_Binary_75kb_colocthresh075"
# prefix="EBV_LCL_Binary_200kb_colocthresh075"
# prefix="EBV_LCL_Binary_300kb_colocthresh075"
# prefix="EBV_LCL_Binary_400kb_colocthresh075"
# prefix="EBV_LCL_Binary_600kb_colocthresh075"
# prefix="EBV_LCL_Binary_700kb_colocthresh075"
# prefix="EBV_LCL_Binary_800kb_colocthresh075"
# prefix="EBV_LCL_Binary_900kb_colocthresh075"

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
args_gene_coord_file = "/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/make_annot_sample_files/ENSG_coord.txt"
#args_gene_coord_file = "/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/make_annot_sample_files/ENSG_coord.txt"

args_windowsize = 50000
args_windowsize = 500000
args_windowsize = 1000000
args_windowsize = 900000

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
    args_bimfile="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_EUR_Phase3_plink/1000G.EUR.QC.{}.bim".format(chrom)
    #args_bimfile="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_EUR_Phase3_plink/1000G.EUR.QC.{}.bim".format(chrom)

    outfilename= prefix + ".{}.annot.gz".format(chrom)
    args_annot_file=join(annotation_dir, outfilename)

    df_list = []
    # annotation_list = ["Null", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]
    # annotation_list = ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]
    # annotation_list = ["fullENSG", "Null", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]
    
    annotation_list = ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]
    # annotation_list = ["Null", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]
    # annotation_list = ["fullENSG", "Null", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]
    
    # annotation_list = ["Null", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "PP.H4.6", "PP.H4.7", "PP.H4.8", "PP.H4.9"]
    # annotation_list = ["Null", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "PP.H4.6", "PP.H4.7", "PP.H4.8", 
    #                     "PP.H4.9", "PP.H4.95", "PP.H4.99", "PP.H4.999", "PP.H4.9999"]
    # annotation_list = ["PP.H0.8", "PP.H1.8", "PP.H2.8", "PP.H3.8", "PP.H4.8", "PP.H4.9", "PP.H4.95"]
    # annotation_list = ["PP.H0.7", "PP.H1.7", "PP.H2.7", "PP.H3.7", "PP.H4.7", "PP.H4.8", "PP.H4.9", "PP.H4.95"]
    # annotation_list = ["PP.H0.75", "PP.H1.75", "PP.H2.75", "PP.H3.75", "PP.H4.75", "PP.H4.8", "PP.H4.9", "PP.H4.95"]
    
    for annotation in annotation_list:
        #annotation = "Null"
        args_gene_set_file = join(output_path, "{}.coloc.abf.txt".format(annotation))
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

    # df_merged = df_merge_idx.fillna(0)
    # df_merged.to_csv(args_annot_file, sep = "\t", index = False, compression="gzip")

    df_merged = df_merge_idx.fillna(0)
    print("\n\nmerging custom annot with baselineLD annotation ...\n\n")
    #df_merged["base"] = 1
    #select_cols = ["CHR", "BP", "SNP", "CM"] + ["base"] + annotation_list
    #df_merged_f = df_merged.loc[:,select_cols]
    baselineLD_dir = "/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/baselineLD_annot"
    df_baselineLD = pd.read_csv(join(baselineLD_dir, "baselineLD.{}.annot.gz".format(chrom)), sep="\t")
    df_merged_baselineLD = pd.merge(df_baselineLD, df_merged, how='left', on=["CHR", "BP", "SNP", "CM"])
    df_merged_baselineLD.to_csv(args_annot_file, sep = "\t", index = False, compression="gzip")

print("Task completed...")

#Sanity check for annotations picked:
# dfnull = df_list[0]
# dfnull[dfnull["Null"] >0]

# dfh0 = df_list[1]
# dfh0[dfh0["PP.H0"] >0]

# dfh1 = df_list[2]
# dfh1[dfh1["PP.H1"] >0]

# dfh2 = df_list[3]
# dfh2[dfh2["PP.H2"] >0]

# dfh3 = df_list[4]
# dfh3[dfh3["PP.H3"] >0]

# dfh4 = df_list[5]
# dfh4[dfh4["PP.H4"] >0]


#############################################
##prepare snp based annot file for ldsc run:
#############################################

#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool
import gzip
import os
from os.path import join

#input filenames:
tissue = "Cells_EBV-transformed_lymphocytes"
work_dir = "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps"
work_dir = "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/high_ldsnps"
#input_file = join(work_dir, tissue + "_coloc_result.txt")

annotpath = join(work_dir, "ldsc_annot")
# prefix="EBV_LCL_Binary_100kb_colocthresh075"
# prefix="EBV_LCL_Binary_100kb_colocthresh07"
prefix="EBV_LCL_Binary_100kb_colocthresh075_nonull"
# prefix="EBV_LCL_Binary_100kb_colocthresh07_nonull"

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
    args_bimfile="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_EUR_Phase3_plink/1000G.EUR.QC.{}.bim".format(chrom)
    #args_bimfile="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_EUR_Phase3_plink/1000G.EUR.QC.{}.bim".format(chrom)

    outfilename= prefix + ".{}.annot.gz".format(chrom)
    args_annot_file=join(annotation_dir, outfilename)

    df_list = []
    # annotation_list = ["PP.H0.75", "PP.H1.75", "PP.H2.75", "PP.H3.75", "PP.H4.75"]
    # annotation_list = ["PP.H0.7", "PP.H1.7", "PP.H2.7", "PP.H3.7", "PP.H4.7"]
    annotation_list = ["PP.H1.75", "PP.H2.75", "PP.H3.75", "PP.H4.75"]
    # annotation_list = ["PP.H1.7", "PP.H2.7", "PP.H3.7", "PP.H4.7"]
    
    for annotation in annotation_list:
        #annotation = "PP.H4.75"
        coloc_geneqtl_file = join(work_dir, "{}.coloc.abf.topsnpList_forLDSC.txt".format(annotation))
        coloc_geneqtl_file = join(work_dir, "{}.coloc.abf.topsnpList_forLDSC.txt.ldSnps.txt".format(annotation))
        #args_gene_set_file = "/home-4/schhetr1@jhu.edu/surya/datasets/prs_project/coloc_output/Coloc_Geneset_EBV_LCC_{}.abf.txt".format(fileidx)

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
    df_merged.to_csv(args_annot_file, sep = "\t", index = False, compression="gzip")

print("Task completed...")


####################
##run s-LDSC program:
####################

freqpath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_Phase3_frq"
#freqpath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_Phase3_frq"

regweights_path="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_Phase3_weights_hm3_no_MHC"
#regweights_path="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_Phase3_weights_hm3_no_MHC"

#############################################
#parameter changes
work_dir="/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps"
work_dir="/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/high_ldsnps"
#############################################

annotpath="${work_dir}/ldsc_annot"

#############################################
#parameter changes
#prefix="EBV_LCL_Binary_100kb_colocthresh075"
#prefix="EBV_LCL_Binary_100kb_colocthresh07"
prefix="EBV_LCL_Binary_100kb_colocthresh075_nonull"
#prefix="EBV_LCL_Binary_100kb_colocthresh07_nonull"
#############################################

annotation_dir="${annotpath}/${prefix}"


outputdir="${annotation_dir}/results"
if [[ ! -d ${outputdir} ]];then mkdir -p ${outputdir}; fi


scriptpath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc"
#scriptpath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc"

bfile_path="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_EUR_Phase3_plink"
#bfile_path="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_EUR_Phase3_plink"

hapmap3_snpspath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/hapmap3_snps"
#hapmap3_snpspath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/hapmap3_snps"

#compute ldscore estimates
for chr in {1..22}; do
    python ${scriptpath}/ldsc.py \
      --l2 \
      --bfile ${bfile_path}/1000G.EUR.QC.${chr} \
      --ld-wind-cm 1 \
      --print-snps ${hapmap3_snpspath}/hm.${chr}.snp \
      --annot ${annotation_dir}/${prefix}.${chr}.annot.gz \
      --out ${annotation_dir}/${prefix}.${chr}
done


sumstats_path="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/my_annot/EBV_LCL_Binary"
#sumstats_path="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/my_annot/EBV_LCL_Binary/"

sumstatsfile="${sumstats_path}/4080_irnt.sbp.ldsc.imputed_v3.both_sexes.tsv.bgz"

genome_1kg_prefix="1000G.EUR.QC"
hm3snp_prefix="weights.hm3_noMHC"

#compute proportion of heritability explained by snps
python ${scriptpath}/ldsc.py \
    --h2 ${sumstatsfile} \
    --ref-ld-chr ${annotation_dir}/${prefix}. \
    --frqfile-chr ${freqpath}/${genome_1kg_prefix}. \
    --w-ld-chr ${regweights_path}/${hm3snp_prefix}. \
    --overlap-annot \
    --print-coefficients \
    --print-delete-vals \
    --out ${outputdir}/${prefix}.ldsc


#####################
#python plotting:
#coloc snps within hign LD of topeQTLs(0.8)

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from os.path import join

sns.set_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8})
sns.set_style("white", {"xtick.major.size": 8, "ytick.major.size": 8})

work_dir = "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/high_ldsnps/ldsc_annot/\
EBV_LCL_Binary_100kb_colocthresh07_nonull/results"
work_dir = "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/high_ldsnps/ldsc_annot/\
EBV_LCL_Binary_100kb_colocthresh075_nonull/results"

infile = join(work_dir, "EBV_LCL_Binary_100kb_colocthresh07_nonull.ldsc.results")
infile = join(work_dir, "EBV_LCL_Binary_100kb_colocthresh075_nonull.ldsc.results")

df = pd.read_csv(infile, sep="\t")
df["Enrichment_p"] = -np.log10(df["Enrichment_p"])
select_cols = ['Category', 'Prop._h2', 'Enrichment', 'Enrichment_p', 'Coefficient_z-score']
df_plot = df.loc[:,select_cols]

category_order = df_plot.sort_values(["Enrichment"]).Category.values.tolist()
#category_order = df_plot.Category.values.tolist()
formatted_order = [each.split(".7")[0] for each in category_order]

fig, axes = plt.subplots(1,4, figsize=(20,5), sharey=True)
fig, axes = plt.subplots(2,4, figsize=(20,5))
axes = axes.flatten()

tips = sns.load_dataset("tips")
iris = sns.load_dataset("iris")
tips.head()
iris.head()

#current_palette = sns.color_palette('colorblind')
sns.barplot(x='Prop._h2', y='Category', data=df_plot, ax=axes[0], order=category_order)
sns.barplot(x='Enrichment', y='Category', data=df_plot, ax=axes[1], order=category_order)

sns.scatterplot(data=tips, x="total_bill", y="tip", ax=axes[2], hue="time")
sns.scatterplot(data=tips, x="total_bill", y="tip", ax=axes[3], hue="size")

df_test = pd.pivot_table(data=sns.load_dataset("flights"),
                    index='month',
                    values='passengers',
                    columns='year')
df_test.head()
sns.heatmap(df_test, ax=axes[4], cmap='coolwarm')
sns.heatmap(df_test, ax=axes[6], cmap='coolwarm', annot=True, cbar=False, annot_kws={'size':1})
sns.heatmap(df_test, ax=axes[5])
sns.scatterplot(data=iris, x="sepal_length", y="sepal_width", ax=axes[3])
sns.scatterplot(data=iris, x="sepal_length", y="petal_length", ax=axes[3], hue="species")
sns.scatterplot(data=iris, x="sepal_length", y="petal_length", ax=axes[5], hue="species")
sns.scatterplot(data=iris, x="sepal_length", y="petal_length", ax=axes[6], hue="species")
sns.scatterplot(data=iris, x="sepal_length", y="petal_length", ax=axes[7], hue="species")

#axes[5].set_axis_off() #remove or delete specific subplot
axes[5].set_visible(False) #remove or delete specific subplot
axes[2].set_visible(False) #remove or delete specific subplot
axes[4].remove()
axes[2].remove()
#fig.delaxes(axes[4])

# Remove the legend and add a colorbar
#https://stackoverflow.com/questions/62884183/trying-to-add-a-colorbar-to-a-seaborn-scatterplot
ax.get_legend().remove()
ax.figure.colorbar(sm)

#adding new subplots directly:
fig.add_subplot(2,4,3)

#Adjusting the spacing of margins and subplots using subplots_adjust.
fig.subplots_adjust(hspace=0, wspace=0.05)

sns.barplot(x='Enrichment_p', y='Category', data=df_plot, ax=axes[7], order=category_order)
sns.barplot(x='Coefficient_z-score', y='Category', data=df_plot, ax=axes[6], order=category_order)

#formatting each axes:
axes[0].set_ylabel("")
axes[1].set_yticklabels("")
axes[2].set_ylabel("")
axes[3].set_ylabel("")
axes[0].set_xlabel("Prop.(heritability)")
axes[2].set_xlabel("-log10(Enrichment Pvalue)")
axes[3].set_xlabel("Coefficient(Z-score)")

# axes[2].legend(fontsize = 5, \
#                bbox_to_anchor= (1.03, 1), \
#                title="Delivery Type", \
#                title_fontsize = 5, \
#                shadow = True, \
#                facecolor = 'white');

axes[2].legend(frameon = False)
axes[3].legend(title="Size", frameon = False)
sns.despine(trim=True, ax=axes[2])
sns.despine(trim=True, offset=True, ax=axes[3])

axes[0].set_yticklabels(formatted_order)

#Create legend handles for single plot
Ax = axes[0]
import matplotlib
Boxes = [item for item in Ax.get_children()
         if isinstance(item, matplotlib.patches.Rectangle)][:-1]

# There is no labels, need to define the labels
legend_labels  = ['EA|AA', 'EA|AA Signal', 'AA Signal', 'EA Signal']
legend_labels  = ['EA|AA', 'AA Signal', 'EA|AA Signal', 'EA Signal']

# Create the legend patches
legend_patches = [matplotlib.patches.Patch(color=C, label=L) for
                  C, L in zip([item.get_facecolor() for item in Boxes],
                              legend_labels)]

# Plot the legend inside plot
plt.legend(handles=legend_patches, frameon = False)

# Plot the legend outside plot
plt.legend(handles=legend_patches, frameon = False, bbox_to_anchor=(1.05, 1), loc='upper left')


output_dir = "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/"
fig.savefig(join(output_dir, 'Figure_highLD_EBV_LCL_Binary_100kb_colocthresh07.pdf'))
fig.savefig(join(output_dir, 'Figure_highLD_EBV_LCL_Binary_100kb_colocthresh075.pdf'))
#plt.tight_layout()

#####################
#python plotting:
#coloc snps of topeQTLs

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from os.path import join

sns.set_style("ticks", {"xtick.major.size": 8, "ytick.major.size": 8})
sns.set_style("white", {"xtick.major.size": 8, "ytick.major.size": 8})

work_dir = "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/\
EBV_LCL_Binary_100kb_colocthresh07_nonull/results"
work_dir = "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/\
EBV_LCL_Binary_100kb_colocthresh075_nonull/results"

infile = join(work_dir, "EBV_LCL_Binary_100kb_colocthresh07_nonull.ldsc.results")
infile = join(work_dir, "EBV_LCL_Binary_100kb_colocthresh075_nonull.ldsc.results")

df = pd.read_csv(infile, sep="\t")
df["Enrichment_p"] = -np.log10(df["Enrichment_p"])
select_cols = ['Category', 'Prop._h2', 'Enrichment', 'Enrichment_p', 'Coefficient_z-score']
df_plot = df.loc[:,select_cols]

category_order = df_plot.sort_values(["Enrichment"]).Category.values.tolist()
#category_order = df_plot.Category.values.tolist()
formatted_order = [each.split(".7")[0] for each in category_order]

fig, axes = plt.subplots(1,4, figsize=(20,5), sharey=True)
axes = axes.flatten()

sns.barplot(x='Prop._h2', y='Category', data=df_plot, ax=axes[0], order=category_order)
sns.barplot(x='Enrichment', y='Category', data=df_plot, ax=axes[1], order=category_order)
sns.barplot(x='Enrichment_p', y='Category', data=df_plot, ax=axes[2], order=category_order)
sns.barplot(x='Coefficient_z-score', y='Category', data=df_plot, ax=axes[3], order=category_order)

#formatting each axes:
axes[1].set_ylabel("")
axes[2].set_ylabel("")
axes[3].set_ylabel("")
axes[0].set_xlabel("Prop.(heritability)")
axes[2].set_xlabel("-log10(Enrichment Pvalue)")
axes[3].set_xlabel("Coefficient(Z-score)")
#sns.despine(offset=True)
axes[0].set_yticklabels(formatted_order)

#Create legend handles for single plot
Ax = axes[0]
import matplotlib
Boxes = [item for item in Ax.get_children()
         if isinstance(item, matplotlib.patches.Rectangle)][:-1]

# There is no labels, need to define the labels
legend_labels  = ['EA|AA', 'EA|AA Signal', 'AA Signal', 'EA Signal']
legend_labels  = ['EA|AA', 'AA Signal', 'EA|AA Signal', 'EA Signal']

# Create the legend patches
legend_patches = [matplotlib.patches.Patch(color=C, label=L) for
                  C, L in zip([item.get_facecolor() for item in Boxes],
                              legend_labels)]

# Plot the legend inside plot
plt.legend(handles=legend_patches, frameon = False)

# Plot the legend outside plot
plt.legend(handles=legend_patches, frameon = False, bbox_to_anchor=(1.05, 1), loc='upper left')

output_dir = "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/"
fig.savefig(join(output_dir, 'Figure_withoutLD_EBV_LCL_Binary_100kb_colocthresh07.pdf'))
fig.savefig(join(output_dir, 'Figure_withoutLD_EBV_LCL_Binary_100kb_colocthresh075.pdf'))
plt.tight_layout()



########################################################################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################


####################
##run s-LDSC program:
####################

freqpath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_Phase3_frq"
#freqpath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_Phase3_frq"

regweights_path="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_Phase3_weights_hm3_no_MHC"
#regweights_path="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_Phase3_weights_hm3_no_MHC"

output_path="/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL"
annotpath="${output_path}/ldsc_annot"

prefix="EBV_LCL_Binary_100kb"
prefix="EBV_LCL_Binary_100kb_minusnull"
prefix="EBV_LCL_Binary_100kb_fullENSGnull"

prefix="EBV_LCL_Binary_500kb_minusnull"
prefix="EBV_LCL_Binary_500kb"
prefix="EBV_LCL_Binary_500kb_fullENSGnull"

prefix="EBV_LCL_Binary_100kb_sequential"
prefix="EBV_LCL_Binary_100kb_sequentialfull"
prefix="EBV_LCL_Binary_100kb_colocthresh08"
prefix="EBV_LCL_Binary_100kb_colocthresh07"

prefix="EBV_LCL_Binary_100kb_colocthresh075"
prefix="EBV_LCL_Binary_50kb_colocthresh075"
prefix="EBV_LCL_Binary_10kb_colocthresh075"
prefix="EBV_LCL_Binary_500kb_colocthresh075"
prefix="EBV_LCL_Binary_1MB_colocthresh075"
prefix="EBV_LCL_Binary_2kb_colocthresh075"
prefix="EBV_LCL_Binary_5kb_colocthresh075"
prefix="EBV_LCL_Binary_25kb_colocthresh075"
prefix="EBV_LCL_Binary_75kb_colocthresh075"
prefix="EBV_LCL_Binary_200kb_colocthresh075"
prefix="EBV_LCL_Binary_300kb_colocthresh075"
prefix="EBV_LCL_Binary_400kb_colocthresh075"
prefix="EBV_LCL_Binary_600kb_colocthresh075"
prefix="EBV_LCL_Binary_700kb_colocthresh075"
prefix="EBV_LCL_Binary_800kb_colocthresh075"
prefix="EBV_LCL_Binary_900kb_colocthresh075"

annotation_dir="${annotpath}/${prefix}"

outputdir="${annotation_dir}/results"
if [[ ! -d ${outputdir} ]];then mkdir -p ${outputdir}; fi


scriptpath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc"
#scriptpath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc"

bfile_path="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_EUR_Phase3_plink"
#bfile_path="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_EUR_Phase3_plink"

hapmap3_snpspath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/hapmap3_snps"
#hapmap3_snpspath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/hapmap3_snps"

#compute ldscore estimates
for chr in {1..22}; do
    python ${scriptpath}/ldsc.py \
      --l2 \
      --bfile ${bfile_path}/1000G.EUR.QC.${chr} \
      --ld-wind-cm 1 \
      --print-snps ${hapmap3_snpspath}/hm.${chr}.snp \
      --annot ${annotation_dir}/${prefix}.${chr}.annot.gz \
      --out ${annotation_dir}/${prefix}.${chr}
done


sumstats_path="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/my_annot/EBV_LCL_Binary"
#sumstats_path="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/my_annot/EBV_LCL_Binary/"

sumstatsfile="${sumstats_path}/4080_irnt.sbp.ldsc.imputed_v3.both_sexes.tsv.bgz"

genome_1kg_prefix="1000G.EUR.QC"
hm3snp_prefix="weights.hm3_noMHC"

#compute proportion of heritability explained by snps
python ${scriptpath}/ldsc.py \
    --h2 ${sumstatsfile} \
    --ref-ld-chr ${annotation_dir}/${prefix}. \
    --frqfile-chr ${freqpath}/${genome_1kg_prefix}. \
    --w-ld-chr ${regweights_path}/${hm3snp_prefix}. \
    --overlap-annot \
    --print-coefficients \
    --print-delete-vals \
    --out ${outputdir}/${prefix}.ldsc



# #############################
# #inputs for ldscore estimates
# #############################

# scriptpath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc"
# #scriptpath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc"

# bfile_path="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_EUR_Phase3_plink"
# #bfile_path="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_EUR_Phase3_plink"

# hapmap3_snpspath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/hapmap3_snps"
# #hapmap3_snpspath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/hapmap3_snps"


# #for chr in {1..22}; do python ${scriptpath}/ldsc.py --l2 --bfile ${bfile_path}/1000G.EUR.QC.${chr} --ld-wind-cm 1 --print-snps ${hapmap3_snpspath}/hm.${chr}.snp --annot ${annotation_dir}/${prefix}.${chr}.annot.gz --out ${annotation_dir}/${prefix}.${chr}; done
# annotpath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/my_annot"
# annotpath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/my_annot"
# prefix="EBV_LCL_Binary"
# prefix="EBV_LCL_Binary_100kb"
# prefix="EBV_LCL_Binary_10kb"
# annotation_dir="${annotpath}/${prefix}"

# # #create annotation dir
# #if [[ ! -d ${annotation_dir} ]];then mkdir -p ${annotation_dir}; fi

# for chr in {1..22}; do
#     python ${scriptpath}/ldsc.py \
#       --l2 \
#       --bfile ${bfile_path}/1000G.EUR.QC.${chr} \
#       --ld-wind-cm 1 \
#       --print-snps ${hapmap3_snpspath}/hm.${chr}.snp \
#       --annot ${annotation_dir}/${prefix}.${chr}.annot.gz \
#       --out ${annotation_dir}/${prefix}.${chr}
# done

cat PP.H.7/PP.H2.7.coloc.abf.topsnpList_forLDSC.txt <(tail -n+2 /home-4/schhetr1@jhu.edu/surya/datasets/1KG/coloc_topsnps/PP.H.7/coloc_topsnps_ldval0.8_afr.ld| awk '{print $3}') <(tail -n+2 /home-4/schhetr1@jhu.edu/surya/datasets/1KG/coloc_topsnps/PP.H.7/coloc_topsnps_ldval0.8_afr.ld| awk '{print $6}') | sort| uniq> coloc_topsnps_ldval0.8_uniqtotals
nps_afr.txt



################################
#render latex equation with python and matplotlib

import matplotlib.pyplot as plt
#\\ indicates new line in latex
a = r'f(x) = Testing: \frac{\exp(-x^2/2)}{\sqrt{2*\pi}}\sum_{n=1}^\infty\frac{a}{a+b+c} \\ \frac{\exp(-x^2/2)}{\sqrt{2*\pi}} \frac{a}{b}'
b = r'f(x) = Newone: \frac{\exp(-x^2/2)}{\sqrt{2*\pi}}\sum_{n=1}^\infty \frac{\exp(-x^2/2)}{\sqrt{2*\pi}} \frac{a}{b}'
ax = plt.axes([0,0,0.3,0.3]) #left,bottom,width,height
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')
plt.text(0.4,1.2,'$%s\\\\%s$'%(a,b), size=15,color="black") #\\\\ escaping each newline char(\\) with 1 slash each
plt.text(0.4,0.5,'$%s\\\\%s$'%(a,b), size=15,color="red")


################################

#################################

#LDPRED FUNCT FILE h2snp COMPUTE

#################################
#FOR single trait at a time:


library(data.table)
library(tidyverse)

#trait of interest
traitname <- "systolic"

#read annotation coefficient file
annocoeff_dir <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot"
df_anno <- fread(file.path(annocoeff_dir, "combined_ALL_results_EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_WITHbaseONLYwithPPH0to4.txt"), sep="\t")

#read ldsc summarystats(1.1M SNPs) to find sample size N and use it for automated trait selection later
ldsc_sumstats_path <- "/work-zfs/abattle4/surya/datasets/prs_project/ldsc/summary_stats_ldsc_format"
ldsctrait_filelist <- list.files(ldsc_sumstats_path, "*.tsv.bgz", full.names=TRUE, recursive = TRUE)
#ldsctrait_file <- "systolic_bloodPressure.4080_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz"
ldsctrait_file <- ldsctrait_filelist[grepl(traitname, ldsctrait_filelist)]
traitname <- paste(unlist(strsplit(basename(ldsctrait_file), "\\."))[1:2], collapse=".")

infile <- file.path(ldsctrait_file)
ldsctrait_df <- read.table(gzfile(infile), header=T, sep="\t") #read gzip file with read.table
N_samplesize <- ldsctrait_df %>% distinct(N) %>% as.numeric

#read full summarystats file (13M SNPs) derived by ldpredfunct_format.py script
#this is to derive the pvalues for all the 1KG annotation SNPs(MC)for later filters(i.e P+T thresholding)
#sumstatsfile <- "/work-zfs/abattle4/surya/datasets/WHIMS/UKBB_sumstats/SBP_eQTLprsformat_summarystats_baseline_ldpredFormat.txt"
sumstatsfile <- "/work-zfs/abattle4/surya/datasets/WHIMS/UKBB_sumstats/SBP_Unique_summarystats_baseline_ldpredFormat.txt"
ssf_df <- fread(sumstatsfile, sep=" ", fill=TRUE) #fill NA for empty elements
sumstats_df <- ssf_df %>% drop_na

#read argsvariable
annotpath <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"
prefix <- "EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_WITHbaseONLYwithPPH0to4"

#ldsc binary annotation 1KG files for regression Coefficient (multiply with MC annotation)
refanno_dir <- file.path(annotpath, prefix)
annotation_list <- list.files(refanno_dir, ".annot.gz", full.names=TRUE, recursive = TRUE)
print(paste("Annotation list:", annotation_list))

#read annotation coefficient file
# df_snp <- fread(file.path(annocoeff_dir, "combined_ALL_results_EBV_LCL_Binary_snpbased_1bp_colocthresh075_WITHbaseONLYwithPPH0to4.txt"), sep="\t")
# df_snp_pph34 <- fread(file.path(annocoeff_dir, "combined_ALL_results_EBV_LCL_Binary_snpbased_1bp_colocthresh075_WITHbaseONLYwithPPH3plus4.txt"), sep="\t")
# df_snp_pph4 <- fread(file.path(annocoeff_dir, "combined_ALL_results_EBV_LCL_Binary_snpbased_1bp_colocthresh075_WITHbaseONLYwithPPH4.txt"), sep="\t")

#df_anno <- fread(file.path(annocoeff_dir, "combined_ALL_results_EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_WITHbaseONLYwithPPH0to4.txt"), sep="\t")
#df_anno_pph34 <- fread(file.path(annocoeff_dir, "combined_ALL_results_EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_WITHbaseONLYwithPPH3plus4.txt"), sep="\t")
#df_anno_pph4 <- fread(file.path(annocoeff_dir, "combined_ALL_results_EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_WITHbaseONLYwithPPH4.txt"), sep="\t")

#extract genomewide heritability h2g
#read string file for genomewide heritabiltiy
# refanno_dir <- file.path(annotpath,prefix)
# stringfile <- "/results_systolic_bloodPressure.4080_irnt/EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_WITHbaseONLYwithPPH0to4.ldsc.log" 
# stringfile <- paste0(refanno_dir, stringfile)
# readfile <- read_file(stringfile)

# #R string capture (Regex group capture in R with multiple capture-groups)
# cisherit <- readfile %>% str_match("scale h2: (0\\.[0-9]+)")
# cisherit <- cisherit[2] %>% as.numeric


#process sldsc log files for cis-genome heritability
filelist <- list.files(file.path(annotpath, prefix), "*ldsc.log", full.names=TRUE, recursive = TRUE)
print(paste("Results_dirlist:", filelist))

trait_list <- list()

for(strfile in filelist) {

    filename <- unlist(strsplit(dirname(strfile), "results_"))[2]
    print(paste("Processing...", filename))

    #load heritabilitly log file for cisherit capture
    #result_file <- file.path(res_dir, paste0(prefix, ".ldsc.log"))
    readfile <- read_file(strfile)

    #String capture 
    #cisherit <- readfile %>% str_match("scale h2: (0\\.[0-9]+)")
    cisherit <- readfile %>% str_match("scale h2: (-?0\\.[0-9]+)") #-? pick if contains -ve val
    cisherit <- cisherit[2] %>% as.numeric
    trait_list[[filename]] <- cisherit

}


#concat or cbind to dataframe
trait_herit_df <- cbind(trait_list) %>% data.frame
names(trait_herit_df) <- "cish2g" #wholegenome
print(head(trait_herit_df))

trait_herit_df$trait <- rownames(trait_herit_df)
cisH2_filter <- trait_herit_df %>% filter(grepl(traitname, trait_herit_df$trait))
cisH2 <- cisH2_filter$cish2g[[1]]

df_anno_merged <- left_join(df_anno, trait_herit_df, by = "trait")
df_anno_h2snp <- df_anno_merged %>% mutate(h2snp=as.numeric(Coefficient)/as.numeric(cish2g))

#compute C1 matrix -- used for muliplication with MC matrix(i.e. SnpsbyAnnotation from the *.annot.gz files)
#filter annotation and construct the matrix (MC and C1)
annocoef_matrix <- df_anno_h2snp %>% filter(grepl(traitname, trait))  %>% select(h2snp) %>% as.matrix
traitname_df <- df_anno_h2snp %>% filter(grepl(traitname, trait))
rownames(annocoef_matrix) <- traitname_df$Category
trait_filename <- traitname_df %>% distinct(trait) %>% as.character


#refanno_dir <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/...
#...cisregion_based/EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_WITHbaseONLYwithPPH0to4"

#process ldsc binary annotation files for regression Coefficient (multiply with MC annotation)
# refanno_dir <- file.path(annotpath,prefix)

# annotation_list <- list.files(refanno_dir, ".annot.gz", full.names=TRUE, recursive = TRUE)
# print(paste("Annotation list:", annotation_list))


#process ldsc binary annotation files for regression Coefficient (multiply with MC annotation)
annot_list <- list()

#computes ldsc funct file persnp-heritability and merges with original ssfs to derive the pvalues for
#1KG annotated SNPs (MCmatrix) for later filtering purpos(similar to P+T thresholding)
for (anno_binaryfile in annotation_list) {

    annotname <- paste0((unlist(strsplit(anno_binaryfile, "\\."))[-1]), collapse=".")
    print(paste("Processing...", annotname))
    #anno_binaryfile <- "EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_WITHbaseONLYwithPPH0to4.20.annot.gz"

    #install.packages('R.utils') #for fread compressed files
    df_snpanno <- fread(anno_binaryfile, sep="\t")
    df_snpanno <- df_snpanno %>% rename(RSID=SNP) %>% mutate(SNP=paste0(CHR, ":", BP))
    snpanno_matrix <- df_snpanno %>% select(-CHR,-BP,-CM,-RSID, -SNP) %>% as.matrix
    rownames(snpanno_matrix) <- df_snpanno$SNP

    #snpanno_matrix %*% annocoef_matrix %>% tail

    #LDmatrix multiplication(MC X C1)
    ldpredfunct_matrix <- snpanno_matrix %*% annocoef_matrix
    ldpredfunct_df <- ldpredfunct_matrix %>% data.frame
    ldpredfunct_df$SNP <- df_snpanno$SNP

    #merge with sumstats to filter based on pvals
    #merges with original ssfs to derive the pvalues for all 1KG annotated SNPS
    merged_df <- inner_join(ldpredfunct_df, sumstats_df, by="SNP")
    annot_list[[annotname]] <- merged_df

}

#concat to dataframe
concat_df <- bind_rows(annot_list) #.id="" for annotname

#filter the final cols for ldpredfunct h2snp file
final_df  <- concat_df %>% select(SNP, h2snp, P, BETA, Z)

#create new dir funct_file
output_dir <- file.path(refanno_dir, "funct_file")

if (!dir.exists(output_dir)){
dir.create(output_dir)
}

#out filename
funct_filename <- file.path(output_dir, paste0(trait_filename, ".ldpredfunct_file.txt"))

#write ldpred funct file
fwrite(final_df, outfile, row.names=F, col.names=T, sep="\t")


#################################################################
#prepare vars for passing args to bash script
scriptpath <- "/work-zfs/abattle4/surya/datasets/prs_project/scripts/ldpredfunct"
ldpredfunct_AA <- "new_ldpredfunct_AA.sh"
ldpredfunct_EUR <- "new_ldpredfunct_EUR.sh"

script_AA <- file.path(scriptpath,ldpredfunct_AA)
script_EUR <- file.path(scriptpath,ldpredfunct_EUR)
funct_file <- funct_filename
h2 <- cisH2

outdir <- annotpath
outfilename <- paste0(prefix, "-", traitname) 

CMD_AA <- paste(script_AA, funct_file, N_samplesize, h2, outdir, outfilename)
CMD_EUR <- paste(script_EUR, funct_file, N_samplesize, h2, outdir, outfilename)

#sanity check
#cat(CMD_EUR)
#cat(CMD_AA)

#run the script
system(CMD_AA)
system(CMD_EUR)





#################################

#LDPRED FUNCT FILE h2snp COMPUTE

#################################
#for looping through summary stats and file annotations

library(data.table)
library(tidyverse)

#trait of interest
traitname <- "systolic"

#read argsvariable
annotpath <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"
#prefix <- "EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_WITHbaseONLYwithPPH0to4"

#read annotation coefficient file
annocoeff_dir <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot"
annocoeff_list <- list.files(annocoeff_dir, "combined.*cisgenebased_500kb.*", full.name=TRUE)

#read ldsc summarystats(1.1M SNPs) to find sample size N and use it for automated trait selection later
ldsc_sumstats_path <- "/work-zfs/abattle4/surya/datasets/prs_project/ldsc/summary_stats_ldsc_format"
ldsctrait_filelist <- list.files(ldsc_sumstats_path, "*.tsv.bgz", full.names=TRUE, recursive = TRUE)

#################################################
#loop through sumstats file
#for ldsctrait_file in ldsctrait_filelist{do ...}
#################################################

#ldsctrait_file <- "systolic_bloodPressure.4080_irnt.ldsc.imputed_v3.both_sexes.tsv.bgz"
ldsctrait_file <- ldsctrait_filelist[grepl(traitname, ldsctrait_filelist)]
traitname <- paste(unlist(strsplit(basename(ldsctrait_file), "\\."))[1:2], collapse=".")

#read ldsc sumstats file(1.1M)
infile <- file.path(ldsctrait_file)
ldsctrait_df <- read.table(gzfile(infile), header=T, sep="\t") #read gzip file with read.table
N_samplesize <- ldsctrait_df %>% distinct(N) %>% as.numeric

#read full summarystats file (13M SNPs) derived by ldpredfunct_format.py script
#this is to derive the pvalues for all the 1KG annotation SNPs(MC)for later filters(i.e P+T thresholding)
sumstatsfile <- "/work-zfs/abattle4/surya/datasets/WHIMS/UKBB_sumstats/SBP_Unique_summarystats_baseline_ldpredFormat.txt"
ssf_df <- fread(sumstatsfile, sep=" ", fill=TRUE) #fill NA for empty elements
sumstats_df <- ssf_df %>% drop_na


for(annocoeff_file in annocoeff_list){

    #annocoeff_file <- annocoeff_list[grepl("WITHbaseONLYwithPPH4", annocoeff_list)]
    df_anno <- fread(annocoeff_file, sep="\t")
    
    #df_anno <- fread(file.path(annocoeff_dir, "combined_ALL_results_EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_WITHbaseONLYwithPPH0to4.txt"), sep="\t")
    prefix <-  unlist(strsplit(basename(annocoeff_file), "results_"))[2] %>% str_replace(".txt", "")

    #ldsc binary annotation 1KG files for regression Coefficient (multiply with MC annotation)
    refanno_dir <- file.path(annotpath, prefix)
    annotation_list <- list.files(refanno_dir, ".annot.gz", full.names=TRUE, recursive = TRUE)
    #print(paste("Annotation list:", annotation_list))


    #process sldsc log files for cis-genome heritability
    filelist <- list.files(file.path(annotpath, prefix), "*ldsc.log", full.names=TRUE, recursive = TRUE)
    #print(paste("Results_dirlist:", filelist))

    trait_list <- list()

    for(strfile in filelist) {

        filename <- unlist(strsplit(dirname(strfile), "results_"))[2]
        print(paste("Processing...", filename))

        #load heritabilitly log file for cisherit capture
        #result_file <- file.path(res_dir, paste0(prefix, ".ldsc.log"))
        readfile <- read_file(strfile)

        #String capture 
        #cisherit <- readfile %>% str_match("scale h2: (0\\.[0-9]+)")
        cisherit <- readfile %>% str_match("scale h2: (-?0\\.[0-9]+)") #-? pick if contains -ve val
        cisherit <- cisherit[2] %>% as.numeric
        trait_list[[filename]] <- cisherit

    }


    #concat or cbind to dataframe
    trait_herit_df <- cbind(trait_list) %>% data.frame
    names(trait_herit_df) <- "cish2g" #wholegenome
    print(head(trait_herit_df))

    trait_herit_df$trait <- rownames(trait_herit_df)
    cisH2_filter <- trait_herit_df %>% filter(grepl(traitname, trait_herit_df$trait))
    cisH2 <- cisH2_filter$cish2g[[1]]

    df_anno_merged <- left_join(df_anno, trait_herit_df, by = "trait")
    df_anno_h2snp <- df_anno_merged %>% mutate(h2snp=as.numeric(Coefficient)/as.numeric(cish2g))

    #compute C1 matrix -- used for muliplication with MC matrix(i.e. SnpsbyAnnotation from the *.annot.gz files)
    #filter annotation and construct the matrix (MC and C1)
    annocoef_matrix <- df_anno_h2snp %>% filter(grepl(traitname, trait))  %>% select(h2snp) %>% as.matrix
    traitname_df <- df_anno_h2snp %>% filter(grepl(traitname, trait))
    rownames(annocoef_matrix) <- traitname_df$Category
    trait_filename <- traitname_df %>% distinct(trait) %>% as.character


    #process ldsc binary annotation files for regression Coefficient (multiply with MC annotation)
    #annot_list <- list()

    #computes ldsc funct file persnp-heritability and merges with original ssfs to derive the pvalues for
    #1KG annotated SNPs (MCmatrix) for later filtering purpos(similar to P+T thresholding)
    for (anno_binaryfile in annotation_list) {

        annotname <- paste0((unlist(strsplit(anno_binaryfile, "\\."))[-1]), collapse=".")
        cat("\n")
        print(paste("Processing annot file ::", basename(anno_binaryfile)))
        cat("\n")
        print(paste("Chrom...", annotname))
        cat("\n")

        #anno_binaryfile <- "EBV_LCL_Binary_cisgenebased_500kb_colocthresh075_WITHbaseONLYwithPPH0to4.20.annot.gz"

        #install.packages('R.utils') #for fread compressed files
        df_snpanno <- fread(anno_binaryfile, sep="\t")
        df_snpanno <- df_snpanno %>% rename(RSID=SNP) %>% mutate(SNP=paste0(CHR, ":", BP))
        snpanno_matrix <- df_snpanno %>% select(-CHR,-BP,-CM,-RSID, -SNP) %>% as.matrix
        rownames(snpanno_matrix) <- df_snpanno$SNP

        #snpanno_matrix %*% annocoef_matrix %>% tail

        #LDmatrix multiplication(MC X C1)
        ldpredfunct_matrix <- snpanno_matrix %*% annocoef_matrix
        ldpredfunct_df <- ldpredfunct_matrix %>% data.frame
        ldpredfunct_df$SNP <- df_snpanno$SNP

        #merge with sumstats to filter based on pvals
        #merges with original ssfs to derive the pvalues for all 1KG annotated SNPS
        merged_df <- inner_join(ldpredfunct_df, sumstats_df, by="SNP")
        annot_list[[annotname]] <- merged_df

    }

    #concat to dataframe
    concat_df <- bind_rows(annot_list) #.id="" for annotname

    #filter the final cols for ldpredfunct h2snp file
    final_df  <- concat_df %>% select(SNP, h2snp, P, BETA, Z)
    percent_overlap <- (nrow(final_df)/ 9997231)*100 #1000KG bim files
    
    cat("\n")
    print(paste0("Percent Overlap: ", round(percent_overlap, 2), "%"))
    cat("\n")
    
    #create new dir funct_file
    output_dir <- file.path(refanno_dir, "funct_file")

    if (!dir.exists(output_dir)){
    dir.create(output_dir)
    }

    #out filename
    funct_filename <- file.path(output_dir, paste0(trait_filename, ".ldpredfunct_file.txt"))

    #write ldpred funct file
    fwrite(final_df, funct_filename, row.names=F, col.names=T, sep="\t")


    #################################################################
    #prepare vars for passing args to simple bash script
    #which will further call on to complex sbatch scripts
    scriptpath <- "/work-zfs/abattle4/surya/datasets/prs_project/scripts/ldpredfunct"
    ldpredfunct_AA <- "new_ldpredfunct_AA.sh"
    ldpredfunct_EUR <- "new_ldpredfunct_EUR.sh"

    script_AA <- file.path(scriptpath,ldpredfunct_AA)
    script_EUR <- file.path(scriptpath,ldpredfunct_EUR)
    funct_file <- funct_filename
    N <- N_samplesize
    h2 <- cisH2

    outdir <- annotpath
    outfilename <- paste0(prefix, "-", traitname) 

    CMD_AA <- paste(script_AA, funct_file, N, h2, outdir, outfilename)
    CMD_EUR <- paste(script_EUR, funct_file, N, h2, outdir, outfilename)

    #sanity check
    #cat(CMD_EUR)
    #cat(CMD_AA)

}



#run the script
system(CMD_AA)
system(CMD_EUR)


################################################
#Extract PRS scores and R2 from LDpredfunct file

library(data.table)
library(tidyverse)

eval "$(conda shell.bash hook)"
conda activate r4base

cohort <- "AA"
prefix <- file.path(paste0("HighMEM45_LDpredfunct_Output_", cohort), "corrected_pvalthresh_prs_baselineLDwithPPH4")
#file <- "ldpredfunct_PRS_systolic_bloodPressure.4080_irnt_WITHbaseONLYwithPPH4"

#process ldpredfunct log files for prs r2
logfile_list <- list.files(file.path(annotpath, aa_prefix), ".*.pvalthresh.*\\.log", full.names=TRUE)

pattern0 <- "Final raw effects PRS r2:"
pattern1 <- "Final LDpred-funct-inf  PRS r2:"
pattern2 <- "Final in-sample LDpredfunct (10 bins) PRS adjusted-R2:"
pattern3 <- "Total number of SNPs:"


parse_compute_persnpR2 <- function(logfile_list, cohort){

    annot_list <-  list()
    category_list <- list()
    raw_prs_list <- list()
    weighted_list <- list()
    predfunct_list <- list()
    snpcount_list <- list()

    for(logfile in logfile_list) {

        #logfile <- logfile_list[1]
        filename <- unlist(strsplit(basename(logfile), "\\.txt"))[1]
        filesplit1 <- unlist(strsplit(filename, ".pvalthresh"))[1]

        annot <- str_replace(filesplit1, "ldpredfunct_PRS_", "")
        category <- unlist(strsplit(filename, ".pvalthresh"))[2] %>% as.numeric
        
        print(paste("Processing...", str_replace(filename, "ldpredfunct_PRS_", "")))

        #load ldpredfunct log file for R2 capture
        readfile <- read_file(logfile)

        #String capture 
        ldpred_r2 <- readfile %>% str_match("Final raw effects PRS r2: (-?0\\.[0-9]+)") #-? pick if contains -ve val
        raw_prs <- ldpred_r2[2] %>% as.numeric

        ldpred_r2 <- readfile %>% str_match("Final LDpred-funct-inf  PRS r2: (-?0\\.[0-9]+)") #-? pick if contains -ve val
        weighted_prs <- ldpred_r2[2] %>% as.numeric

        ldpred_r2 <- readfile %>% str_match("Final in-sample LDpredfunct \\(10 bins\\) PRS adjusted-R2: (-?0\\.[0-9]+)") #-? pick if contains -ve val
        predfunct_prs <- ldpred_r2[2] %>% as.numeric

        ldpred_snpcount <- readfile %>% str_match("Total number of SNPs: ([0-9]+)") #-? pick if contains -ve val
        snpcount <- ldpred_snpcount[2] %>% as.numeric

        trait_id <- paste0(annot, "pval", category) 
        
        annot_list[[trait_id]] <- annot
        category_list[[trait_id]] <- category
        raw_prs_list[[trait_id]] <- raw_prs
        weighted_list[[trait_id]] <- weighted_prs
        predfunct_list[[trait_id]] <- predfunct_prs
        snpcount_list[[trait_id]] <- snpcount

    }

    bindlist <- list(category_list, raw_prs_list, weighted_list, predfunct_list, snpcount_list, annot_list) 
    bind_df <- do.call(cbind, bindlist) %>% as.data.frame
    names(bind_df) <- c("pvalcategory", "raw_prs", "weighted_prs", "predfunct_prs", "snpcount", "annotation")
    final_df <- bind_df %>% arrange(pvalcategory)
    final_df["cohort"] <- cohort

    #compute R2perSNP
    #replace nans with 1 to filter out later
    final_df <- final_df %>% replace(is.na(.), 1)
    predfunct_prs <- final_df %>% mutate(RawPRS_R2perSNP = as.numeric(raw_prs)/as.numeric(snpcount),
                                        WeightedPRS_R2perSNP = as.numeric(weighted_prs)/as.numeric(snpcount),
                                        AdjFunctPRS_R2perSNP = as.numeric(predfunct_prs)/as.numeric(snpcount)
                                        )

    return(predfunct_prs)

}

############
############
#Parse and compute R2perSNP for baselineLDwithPPH4
#maindir
annotpath <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"

#process ldpredfunct log files for prs r2
cohort <- "AA"
aa_prefix <- file.path(paste0("HighMEM45_LDpredfunct_Output_", cohort), "corrected_pvalthresh_prs_baselineLDwithPPH4")
aa_logfile_list <- list.files(file.path(annotpath, aa_prefix), ".*.pvalthresh.*\\.log", full.names=TRUE)

#call function
aa_df <- parse_compute_persnpR2(aa_logfile_list, cohort)

cohort <- "EUR"
eur_prefix <- file.path(paste0("HighMEM45_LDpredfunct_Output_", cohort), "corrected_pvalthresh_prs_baselineLDwithPPH4")
eur_logfile_list <- list.files(file.path(annotpath, eur_prefix), ".*.pvalthresh.*\\.log", full.names=TRUE)

#call function
eur_df <- parse_compute_persnpR2(eur_logfile_list, cohort)

#concat aa and eur dfs
predfunct_df1 <- bind_rows(aa_df, eur_df)
rownames(predfunct_df1) <- seq(1:nrow(predfunct_df1))

#quick check
baselineLDpph4 <- predfunct_df1
baselineLDpph4 %>% filter(pvalcategory == 1)

outfilename <- file.path(annotpath, "baselineLDwithPPH4_ldpredfunct_PRS_parsedR2.txt")
fwrite(baselineLDpph4,  outfilename, sep="\t", col.names=TRUE, row.names=F)




#############
#############
#Parse and compute R2perSNP for WITHbaseONLYwithPPH4
#maindir
annotpath <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"

#process ldpredfunct log files for prs r2
cohort <- "AA"
aa_prefix <- file.path(paste0("HighMEM45_LDpredfunct_Output_", cohort), "corrected_pvalthresh_prs")
aa_logfile_list <- list.files(file.path(annotpath, aa_prefix), ".*.pvalthresh.*\\.log", full.names=TRUE)

#call function
aa_df <- parse_compute_persnpR2(aa_logfile_list, cohort)

cohort <- "EUR"
eur_prefix <- file.path(paste0("HighMEM45_LDpredfunct_Output_", cohort), "corrected_pvalthresh_prs")
eur_logfile_list <- list.files(file.path(annotpath, eur_prefix), ".*.pvalthresh.*\\.log", full.names=TRUE)

#call function
eur_df <- parse_compute_persnpR2(eur_logfile_list, cohort)

#concat aa and eur dfs
predfunct_df2 <- bind_rows(aa_df, eur_df)
rownames(predfunct_df2) <- seq(1:nrow(predfunct_df2))

#quick check
baseOnlypph4 <- predfunct_df2
baseOnlypph4 %>% filter(pvalcategory == 1)

outfilename <- file.path(annotpath, "WITHbaseOnlyPPH4_ldpredfunct_PRS_parsedR2.txt")
fwrite(baseOnlypph4,  outfilename, sep="\t", col.names=TRUE, row.names=F)


###################################################################
#Compute PRS R2 from scratch using multiple linear regression model



compute_adjR2_covariates <- function(prsfile_list, cohort, phenotype_df){

    annot_list <-  list()
    category_list <- list()
    raw_prs_list <- list()
    weighted_list <- list()
    predfunct_list <- list()
    snpcount_list <- list()

    for(prs_file in prsfile_list){

        #prs_file <- aa_prsfile_list[8]
        filename <- unlist(strsplit(basename(prs_file), "\\.txt"))[1]
        filesplit1 <- unlist(strsplit(filename, ".pvalthresh"))[1]

        annot <- str_replace(filesplit1, "ldpredfunct_PRS_", "")
        category <- unlist(strsplit(filename, ".pvalthresh"))[2] %>% as.numeric
        trait_id <- paste0(annot, "pval", category) 
        
        print(paste("Processing...", str_replace(filename, "ldpredfunct_PRS_", "")))

        #load prs file
        readprsfile <- fread(prs_file, sep=",") %>% as.data.frame

        #merge with phenotype files and use as lm model ref dataframe
        merged_df <- inner_join(phenotype_df, readprsfile, by="IID")

        #build final model with PCs and Age Covars
        raw_final_df <- merged_df %>% select(SBP, raw_effects_prs, age, BMI, PC1:PC10)
        raw_final_model <- lm(SBP ~ ., raw_final_df)

        #without BMI i.e with "-BMI" into the model directly
        #final_model <- lm(SBP ~ . -BMI, final_df)

        #summaries
        summary(raw_final_model)
        raw_final_coefs <- summary(raw_final_model)$coef[2,1]
        raw_final_coefs_pval <- summary(raw_final_model)$coef[2,4]
        raw_final_r2 <- summary(raw_final_model)$r.squared
        raw_prs_adjr2 <- summary(raw_final_model)$adj.r.squared

        #build final model with PCs and Age Covars
        weighted_final_df <- merged_df %>% select(SBP, pval_derived_effects_prs, age, BMI, PC1:PC10)
        weighted_final_model <- lm(SBP ~ ., weighted_final_df)

        #summaries
        summary(weighted_final_model)
        weighted_final_coefs <- summary(weighted_final_model)$coef[2,1]
        weighted_final_coefs_pval <- summary(weighted_final_model)$coef[2,4]
        weighted_final_r2 <- summary(weighted_final_model)$r.squared
        weighted_prs_adjr2 <- summary(weighted_final_model)$adj.r.squared

        #build the final 10bins-ldpredfunct model
        predfunct_df <- merged_df %>% select(SBP, pval_derived_effects_prs, age, BMI, PC1:PC10, Bin_1:Bin_10)
        #automate selection for multiple vars linear regression
        predfunct_model <- lm(SBP ~ . , predfunct_df)

        #summaries
        summary(predfunct_model)
        predfunct_coefs <- summary(predfunct_model)$coef[2,1]
        predfunct_coefs_pval <- summary(predfunct_model)$coef[2,4]
        predfunct_r2 <- summary(predfunct_model)$r.squared
        predfunct_prs_adjr2 <- summary(predfunct_model)$adj.r.squared

        annot_list[[trait_id]] <- annot
        category_list[[trait_id]] <- category
        raw_prs_list[[trait_id]] <- raw_prs_adjr2
        weighted_list[[trait_id]] <- weighted_prs_adjr2
        predfunct_list[[trait_id]] <- predfunct_prs_adjr2
        #snpcount_list[[trait_id]] <- snpcount

    }

    #for snpcount merge with previously parsed ldpredfunct file logs
    bindlist <- list(category_list, raw_prs_list, weighted_list, predfunct_list, annot_list) 
    bind_df <- do.call(cbind, bindlist) %>% as.data.frame
    names(bind_df) <- c("pvalcategory", "raw_prs_adj", "weighted_prs_adj", "predfunct_prs_adj", "annotation")
    predfunct_df <- bind_df %>% arrange(pvalcategory)
    predfunct_df["cohort"] <- cohort

    #compute R2perSNP
    #replace nans with 1 to filter out later
    predfunct_df <- predfunct_df %>% replace(is.na(.), 1)

    return(predfunct_df)

}

############
############
#Parse and compute R2perSNP for baselineLDwithPPH4
#maindir
annotpath <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"
trait <- "SBP"

cohort <- "AA"
annotpath <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"
aa_prefix <- file.path(paste0("HighMEM45_LDpredfunct_Output_", cohort), "corrected_pvalthresh_prs_baselineLDwithPPH4")
aa_prsfile_list <- list.files(file.path(annotpath, aa_prefix), ".*.pvalthresh.*\\.*validation.*", full.names=TRUE)

#read phenotype_files
phenotype_aa <- "/work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/PAGE/race3/unrelated/pheno_mendrand.txt"
phenotype_aa_df <- fread(phenotype_aa, sep="\t")

#remove nan trait values
phenotype_aa_df <- phenotype_aa_df %>% filter(!is.na(!!as.name(trait))) #!!for passing argvars to dplyr

#call function
aa_df <- compute_adjR2_covariates(aa_prsfile_list, cohort, phenotype_aa_df)


cohort <- "EUR"
annotpath <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"
eur_prefix <- file.path(paste0("HighMEM45_LDpredfunct_Output_", cohort), "corrected_pvalthresh_prs_baselineLDwithPPH4")
eur_prsfile_list <- list.files(file.path(annotpath, eur_prefix), ".*.pvalthresh.*\\.*validation.*", full.names=TRUE)

#read phenotype files
phenotype_eur <- "/work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/WHIMS/race5/unrelated/pheno_mendrand.txt"
phenotype_eur_df <- fread(phenotype_eur, sep="\t")

#remove nan trait values
phenotype_eur_df <- phenotype_eur_df %>% filter(!is.na(!!as.name(trait))) #!!for passing argvars to dplyr

#call function
eur_df <- compute_adjR2_covariates(eur_prsfile_list, cohort, phenotype_eur_df)

#concat aa and eur dfs
predfunct_prs_df1 <- bind_rows(aa_df, eur_df)
rownames(predfunct_prs_df1) <- seq(1:nrow(predfunct_prs_df1))

#quick check
baselineLDpph4_prs <- predfunct_prs_df1
baselineLDpph4_prs %>% filter(pvalcategory == 1)

outfilename <- file.path(annotpath, "baselineLDwithPPH4_ldpredfunct_PRS_computedR2.txt")
fwrite(baselineLDpph4_prs,  outfilename, sep="\t", col.names=TRUE, row.names=F)





#############
#############
#Parse and compute R2perSNP for WITHbaseONLYwithPPH4
#maindir
annotpath <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"
trait <- "SBP"

cohort <- "AA"
annotpath <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"
aa_prefix <- file.path(paste0("HighMEM45_LDpredfunct_Output_", cohort), "corrected_pvalthresh_prs")
aa_prsfile_list <- list.files(file.path(annotpath, aa_prefix), ".*.pvalthresh.*\\.*validation.*", full.names=TRUE)

#read phenotype_files
phenotype_aa <- "/work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/PAGE/race3/unrelated/pheno_mendrand.txt"
phenotype_aa_df <- fread(phenotype_aa, sep="\t")

#remove nan trait values
phenotype_aa_df <- phenotype_aa_df %>% filter(!is.na(!!as.name(trait))) #!!for passing argvars to dplyr

#call function
aa_df <- compute_adjR2_covariates(aa_prsfile_list, cohort, phenotype_aa_df)


cohort <- "EUR"
annotpath <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"
eur_prefix <- file.path(paste0("HighMEM45_LDpredfunct_Output_", cohort), "corrected_pvalthresh_prs")
eur_prsfile_list <- list.files(file.path(annotpath, eur_prefix), ".*.pvalthresh.*\\.*validation.*", full.names=TRUE)

#read phenotype files
phenotype_eur <- "/work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/WHIMS/race5/unrelated/pheno_mendrand.txt"
phenotype_eur_df <- fread(phenotype_eur, sep="\t")

#remove nan trait values
phenotype_eur_df <- phenotype_eur_df %>% filter(!is.na(!!as.name(trait))) #!!for passing argvars to dplyr

#call function
eur_df <- compute_adjR2_covariates(eur_prsfile_list, cohort, phenotype_eur_df)

#concat aa and eur dfs
predfunct_prs_df2 <- bind_rows(aa_df, eur_df)
rownames(predfunct_prs_df2) <- seq(1:nrow(predfunct_prs_df2))

#quick check
baseOnlypph4_prs <- predfunct_prs_df2
baseOnlypph4_prs %>% filter(pvalcategory == 1)

outfilename <- file.path(annotpath, "WITHbaseOnlyPPH4_ldpredfunct_PRS_computedR2.txt")
fwrite(baseOnlypph4_prs,  outfilename, sep="\t", col.names=TRUE, row.names=F)



#################
#################
#Merge SNP count for persnpR2

#baselineLDwithPPH4
r2_snpcount_merged <- inner_join(baselineLDpph4_prs, baselineLDpph4, 
                                by=c("annotation" = "annotation", "cohort" = "cohort", "pvalcategory" = "pvalcategory"))

#compute R2perSNP
#replace nans with 1 to filter out later
#final_df <- final_df %>% replace(is.na(.), 1)
final_predfunct_prs <- r2_snpcount_merged %>% mutate(RawPRS_R2perSNP_adj = as.numeric(raw_prs_adj)/as.numeric(snpcount),
                                    WeightedPRS_R2perSNP_adj = as.numeric(weighted_prs_adj)/as.numeric(snpcount),
                                    AdjFunctPRS_R2perSNP_adj = as.numeric(predfunct_prs_adj)/as.numeric(snpcount)
                                    )
#quick check
final_predfunct_prs %>% filter(pvalcategory == 1)

outfilename <- file.path(annotpath, "baselineLDPPH4_Final_ldpredfunct_PRS_computedR2_plus_persnpR2.txt")
fwrite(final_predfunct_prs,  outfilename, sep="\t", col.names=TRUE, row.names=F)


# #select specific cols
# select_df <- final_predfunct_prs %>% select(contains("_adj"))
# #select_df <- final_predfunct_prs %>% select(ends_with("_adj"))
# #select_df <- final_predfunct_prs %>% select(starts_with("_adj"))


#WITHbaseOnlyPPH4
r2_snpcount_merged <- inner_join(baseOnlypph4_prs, baseOnlypph4, 
                                by=c("annotation" = "annotation", "cohort" = "cohort", "pvalcategory" = "pvalcategory"))

#compute R2perSNP
#replace nans with 1 to filter out later
#final_df <- final_df %>% replace(is.na(.), 1)
final_predfunct_prs <- r2_snpcount_merged %>% mutate(RawPRS_R2perSNP_adj = as.numeric(raw_prs_adj)/as.numeric(snpcount),
                                    WeightedPRS_R2perSNP_adj = as.numeric(weighted_prs_adj)/as.numeric(snpcount),
                                    AdjFunctPRS_R2perSNP_adj = as.numeric(predfunct_prs_adj)/as.numeric(snpcount)
                                    )
#quick check
final_predfunct_prs %>% filter(pvalcategory == 1)

outfilename <- file.path(annotpath, "WITHbaseOnlyPPH4_Final_ldpredfunct_PRS_computedR2_plus_persnpR2.txt")
fwrite(final_predfunct_prs,  outfilename, sep="\t", col.names=TRUE, row.names=F)





###############################################################################

#Plot PRS graphs


###############################################################################

#load previously written Final_ldpredfunct_PRS_computedR2_plus_persnpR2.txt file

annotpath <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"
outfilename <- file.path(annotpath, "Final_ldpredfunct_PRS_computedR2_plus_persnpR2.txt")

#read prs file
prs_df <- fread(outfilename, sep="\t")





















# ############
# #Parse and compute R2perSNP for baselineLDwithPPH4
# #maindir
# annotpath <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"

# #process ldpredfunct log files for prs r2
# cohort <- "AA"
# aa_prefix <- file.path(paste0("HighMEM45_LDpredfunct_Output_", cohort), "corrected_pvalthresh_prs_baselineLDwithPPH4")
# aa_logfile_list <- list.files(file.path(annotpath, aa_prefix), ".*.pvalthresh.*\\.log", full.names=TRUE)

# #call function
# aa_df <- parse_compute_persnpR2(aa_logfile_list, cohort)

# cohort <- "EUR"
# eur_prefix <- file.path(paste0("HighMEM45_LDpredfunct_Output_", cohort), "corrected_pvalthresh_prs_baselineLDwithPPH4")
# eur_logfile_list <- list.files(file.path(annotpath, eur_prefix), ".*.pvalthresh.*\\.log", full.names=TRUE)

# #call function
# eur_df <- parse_compute_persnpR2(eur_logfile_list, cohort)

# #concat aa and eur dfs
# predfunct_df1 <- bind_rows(aa_df, eur_df)
# rownames(predfunct_df1) <- seq(1:nrow(predfunct_df1))

# baselineLDpph4 <- predfunct_df1
# baselineLDpph4 %>% filter(pvalcategory == 1)


# Final raw effects PRS correlation: -0.0448
# Final raw effects PRS r2: 0.0020

# Final LDpred-funct-inf PRS correlation: -0.0787
# Final LDpred-funct-inf  PRS r2: 0.0062

# Final in-sample LDpredfunct (10 bins) PRS correlation: 0.1036
# Final in-sample LDpredfunct (10 bins) PRS R2: 0.0107
# Final in-sample LDpredfunct (10 bins) PRS adjusted-R2: 0.0092


sanity_check_compute_adjR2_covariates <- function(prsfile_list, cohort, phenotype_df){

    annot_list <-  list()
    category_list <- list()
    raw_prs_list <- list()
    weighted_list <- list()
    predfunct_list <- list()
    snpcount_list <- list()

    for(prs_file in prsfile_list){

        prs_file <- aa_prsfile_list[8]
        filename <- unlist(strsplit(basename(prs_file), "\\.txt"))[1]
        filesplit1 <- unlist(strsplit(filename, ".pvalthresh"))[1]

        annot <- str_replace(filesplit1, "ldpredfunct_PRS_", "")
        category <- unlist(strsplit(filename, ".pvalthresh"))[2] %>% as.numeric
        trait_id <- paste0(annot, "pval", category) 
        
        print(paste("Processing...", str_replace(filename, "ldpredfunct_PRS_", "")))

        #load prs file
        readprsfile <- fread(prs_file, sep=",") %>% as.data.frame

        #merge with phenotype files and use as lm model ref dataframe
        merged_df <- inner_join(phenotype_aa_df, readprsfile, by="IID")

        #build the raw model
        raw_model <- lm(SBP ~ raw_effects_prs, merged_df)

        #summaries
        summary(raw_model)
        raw_coefs <- summary(raw_model)$coef[2,1]
        raw_coefs_pval <- summary(raw_model)$coef[2,4]
        raw_r2 <- summary(raw_model)$r.squared
        raw_adjr2 <- summary(raw_model)$adj.r.squared

        #build the weighted model
        weighted_model <- lm(SBP ~ pval_derived_effects_prs, merged_df)

        #summaries
        summary(weighted_model)
        weighted_coefs <- summary(raw_model)$coef[2,1]
        weighted_coefs_pval <- summary(raw_model)$coef[2,4]
        weighted_r2 <- summary(raw_model)$r.squared
        weighted_adjr2 <- summary(raw_model)$adj.r.squared
        
        #build the 10bins-ldpredfunct model
        selected_df <- merged_df %>% select(SBP, pval_derived_effects_prs, Bin_1:Bin_10)
        #automate selection for multiple vars linear regression
        predfunct_model <- lm(SBP ~ . , selected_df)

        #summaries
        summary(predfunct_model)
        predfunct_coefs <- summary(predfunct_model)$coef[2,1]
        predfunct_coefs_pval <- summary(predfunct_model)$coef[2,4]
        predfunct_r2 <- summary(predfunct_model)$r.squared
        predfunct_adjr2 <- summary(predfunct_model)$adj.r.squared

        #or, build the 10bins alternatively with paste function
        bins <- paste0("Bin", "_", seq(1:10), collapse=" + ")
        expr <- paste0("SBP ~ pval_derived_effects_prs", " + ", bins)
        covars <- as.formula(expr)

        #run lm model
        bin_model <- lm(covars, merged_df)
        summary(bin_model)

        #or, run equivalent manual lm model
        test_model <- lm(SBP ~ pval_derived_effects_prs + Bin_1 + Bin_2 + Bin_3 + Bin_4 + Bin_5 + 
                        Bin_6 + Bin_7 + Bin_8 + Bin_9 + Bin_10, merged_df)
        summary(test_model)



        #build final model with PCs and Age Covars
        raw_final_df <- merged_df %>% select(SBP, raw_effects_prs, age, BMI, PC1:PC10)
        raw_final_model <- lm(SBP ~ ., raw_final_df)

        #without BMI i.e with "-BMI" into the model directly
        #raw_final_model <- lm(SBP ~ . -age, raw_final_df)

        #summaries
        summary(raw_final_model)
        raw_final_coefs <- summary(raw_final_model)$coef[2,1]
        raw_final_coefs_pval <- summary(raw_final_model)$coef[2,4]
        raw_final_r2 <- summary(raw_final_model)$r.squared
        raw_prs_adjr2 <- summary(raw_final_model)$adj.r.squared

        #build final model with PCs and Age Covars
        weighted_final_df <- merged_df %>% select(SBP, pval_derived_effects_prs, age, BMI, PC1:PC10)
        weighted_final_model <- lm(SBP ~ ., weighted_final_df)

        #summaries
        summary(weighted_final_model)
        weighted_final_coefs <- summary(weighted_final_model)$coef[2,1]
        weighted_final_coefs_pval <- summary(weighted_final_model)$coef[2,4]
        weighted_final_r2 <- summary(weighted_final_model)$r.squared
        weighted_prs_adjr2 <- summary(weighted_final_model)$adj.r.squared

        #build the final 10bins-ldpredfunct model
        predfunct_df <- merged_df %>% select(SBP, pval_derived_effects_prs, age, BMI, PC1:PC10, Bin_1:Bin_10)
        #automate selection for multiple vars linear regression
        predfunct_model <- lm(SBP ~ . , predfunct_df)

        #summaries
        summary(predfunct_model)
        predfunct_coefs <- summary(predfunct_model)$coef[2,1]
        predfunct_coefs_pval <- summary(predfunct_model)$coef[2,4]
        predfunct_r2 <- summary(predfunct_model)$r.squared
        predfunct_prs_adjr2 <- summary(predfunct_model)$adj.r.squared

        annot_list[[trait_id]] <- annot
        category_list[[trait_id]] <- category
        raw_prs_list[[trait_id]] <- raw_prs_adjr2
        weighted_list[[trait_id]] <- weighted_prs_adjr2
        predfunct_list[[trait_id]] <- predfunct_prs_adjr2
        #snpcount_list[[trait_id]] <- snpcount

    }

    #for snpcount merge with previously parsed ldpredfunct file logs
    bindlist <- list(category_list, raw_prs_list, weighted_list, predfunct_list, annot_list) 
    bind_df <- do.call(cbind, bindlist) %>% as.data.frame
    names(bind_df) <- c("pvalcategory", "raw_prs_adj", "weighted_prs_adj", "predfunct_prs_adj", "annotation")
    predfunct_df <- bind_df %>% arrange(pvalcategory)
    predfunct_df["cohort"] <- cohort

    #compute R2perSNP
    #replace nans with 1 to filter out later
    predfunct_df <- predfunct_df %>% replace(is.na(.), 1)

    return(predfunct_df)

}






# #############
# #############
# #Parse and compute R2perSNP for WITHbaseONLYwithPPH4
# #maindir
# annotpath <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/cisregion_based"

# #process ldpredfunct log files for prs r2
# cohort <- "AA"
# aa_prefix <- file.path(paste0("HighMEM45_LDpredfunct_Output_", cohort), "corrected_pvalthresh_prs")
# aa_logfile_list <- list.files(file.path(annotpath, aa_prefix), ".*.pvalthresh.*\\.log", full.names=TRUE)

# #call function
# aa_df <- parse_compute_persnpR2(aa_logfile_list, cohort)

# cohort <- "EUR"
# eur_prefix <- file.path(paste0("HighMEM45_LDpredfunct_Output_", cohort), "corrected_pvalthresh_prs")
# eur_logfile_list <- list.files(file.path(annotpath, eur_prefix), ".*.pvalthresh.*\\.log", full.names=TRUE)

# #call function
# eur_df <- parse_compute_persnpR2(eur_logfile_list, cohort)

# #concat aa and eur dfs
# predfunct_df2 <- bind_rows(aa_df, eur_df)
# rownames(predfunct_df2) <- seq(1:nrow(predfunct_df2))

# baseOnlypph4 <- predfunct_df2
# baseOnlypph4 %>% filter(pvalcategory == 1)


        # #build the raw model
        # raw_model <- lm(SBP ~ raw_effects_prs, merged_df)

        # #summaries
        # summary(raw_model)
        # raw_coefs <- summary(raw_model)$coef[2,1]
        # raw_coefs_pval <- summary(raw_model)$coef[2,4]
        # raw_r2 <- summary(raw_model)$r.squared
        # raw_adjr2 <- summary(raw_model)$adj.r.squared

        # #build the weighted model
        # weighted_model <- lm(SBP ~ pval_derived_effects_prs, merged_df)

        # #summaries
        # summary(weighted_model)
        # weighted_coefs <- summary(raw_model)$coef[2,1]
        # weighted_coefs_pval <- summary(raw_model)$coef[2,4]
        # weighted_r2 <- summary(raw_model)$r.squared
        # weighted_adjr2 <- summary(raw_model)$adj.r.squared
        
        # #build the 10bins-ldpredfunct model
        # selected_df <- merged_df %>% select(SBP, pval_derived_effects_prs, Bin_1:Bin_10)
        # #automate selection for multiple vars linear regression
        # predfunct_model <- lm(SBP ~ . , selected_df)

        # #summaries
        # summary(predfunct_model)
        # predfunct_coefs <- summary(predfunct_model)$coef[2,1]
        # predfunct_coefs_pval <- summary(predfunct_model)$coef[2,4]
        # predfunct_r2 <- summary(predfunct_model)$r.squared
        # predfunct_adjr2 <- summary(predfunct_model)$adj.r.squared

        # #or, build the 10bins alternatively with paste function
        # bins <- paste0("Bin", "_", seq(1:10), collapse=" + ")
        # expr <- paste0("SBP ~ pval_derived_effects_prs", " + ", bins)
        # covars <- as.formula(expr)

        # #run lm model
        # bin_model <- lm(covars, merged_df)
        # summary(bin_model)

        # #or, run equivalent manual lm model
        # test_model <- lm(SBP ~ pval_derived_effects_prs + Bin_1 + Bin_2 + Bin_3 + Bin_4 + Bin_5 + 
        #                 Bin_6 + Bin_7 + Bin_8 + Bin_9 + Bin_10, merged_df)
        # summary(test_model)


#######################################################

        # #recompute ldpredfunct r2 for validations
        # #build the raw model
        # raw_model <- lm(SBP ~ raw_effects_prs, merged_df)

        # #summaries
        # summary(raw_model)
        # raw_coefs <- summary(raw_model)$coef[2,1]
        # raw_coefs_pval <- summary(raw_model)$coef[2,4]
        # raw_r2 <- summary(raw_model)$r.squared
        # raw_adjr2 <- summary(raw_model)$adj.r.squared

        # #build the weighted model
        # weighted_model <- lm(SBP ~ pval_derived_effects_prs, merged_df)

        # #summaries
        # summary(weighted_model)
        # weighted_coefs <- summary(raw_model)$coef[2,1]
        # weighted_coefs_pval <- summary(raw_model)$coef[2,4]
        # weighted_r2 <- summary(raw_model)$r.squared
        # weighted_adjr2 <- summary(raw_model)$adj.r.squared
        
        # #build the 10bins-ldpredfunct model
        # selected_df <- merged_df %>% select(SBP, pval_derived_effects_prs, Bin_1:Bin_10)
        # #automate selection for multiple vars linear regression
        # predfunct_model <- lm(SBP ~ . , selected_df)

        # #summaries
        # summary(predfunct_model)
        # predfunct_coefs <- summary(predfunct_model)$coef[2,1]
        # predfunct_coefs_pval <- summary(predfunct_model)$coef[2,4]
        # predfunct_r2 <- summary(predfunct_model)$r.squared
        # predfunct_adjr2 <- summary(predfunct_model)$adj.r.squared

        # #or, build the 10bins alternatively with paste function
        # bins <- paste0("Bin", "_", seq(1:10), collapse=" + ")
        # expr <- paste0("SBP ~ pval_derived_effects_prs", " + ", bins)
        # covars <- as.formula(expr)

        # #run lm model
        # bin_model <- lm(covars, merged_df)
        # summary(bin_model)

        # #or, run equivalent manual lm model
        # test_model <- lm(SBP ~ pval_derived_effects_prs + Bin_1 + Bin_2 + Bin_3 + Bin_4 + Bin_5 + 
        #                 Bin_6 + Bin_7 + Bin_8 + Bin_9 + Bin_10, merged_df)
        # summary(test_model)

########################################################






# > snpanno_matrix %>% tail
#             base PP.H0.75 PP.H1.75 PP.H2.75 PP.H3.75 PP.H4.75
# rs73628077     1        0        0        0        1        1
# rs568825360    1        0        0        0        1        1
# rs11698187     1        0        0        0        1        1
# rs186809596    1        0        0        0        1        1
# rs140775622    1        0        0        0        1        1
# rs542224338    1        0        0        0        1        1
# > annocoef_matrix %>% tail
#              Enrichment
# baseL2_0       1.000000
# PP.H0.75L2_0   1.161670
# PP.H1.75L2_0   1.559310
# PP.H2.75L2_0   1.425258
# PP.H3.75L2_0   1.376064
# PP.H4.75L2_0   1.277765
# > snpanno_matrix %*% annocoef_matrix %>% tail
#             Enrichment
# rs73628077    3.653829
# rs568825360   3.653829
# rs11698187    3.653829
# rs186809596   3.653829
# rs140775622   3.653829
# rs542224338   3.653829
# > 1.000000 + 1.376064 + 1.277765
# [1] 3.653829


# snpanno_matrix %*% annocoef_matrix %>% head
#             Enrichment
# rs6078030     2.376064
# rs143291093   2.376064
# rs4814683     2.376064
# rs34147676    2.376064
# rs6076506     2.376064
# rs6139074     2.376064
# > annocoef_matrix
#              Enrichment
# baseL2_0       1.000000
# PP.H0.75L2_0   1.161670
# PP.H1.75L2_0   1.559310
# PP.H2.75L2_0   1.425258
# PP.H3.75L2_0   1.376064
# PP.H4.75L2_0   1.277765
# > snpanno_matrix %>% head
#             base PP.H0.75 PP.H1.75 PP.H2.75 PP.H3.75 PP.H4.75
# rs6078030      1        0        0        0        1        0
# rs143291093    1        0        0        0        1        0
# rs4814683      1        0        0        0        1        0
# rs34147676     1        0        0        0        1        0
# rs6076506      1        0        0        0        1        0
# rs6139074      1        0        0        0        1        0

