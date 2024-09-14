
#LDPREDFunct format
import pandas as pd
df = pd.read_csv("4080_irnt.sbp.gwas.imputed_v3.both_sexes.tsv", sep ="\t")

df[["CHR", "BP", "A2", "A1"]] = df["variant"].str.split(":", expand=True)
df["Z"] = df["beta"]/df["se"]
df.rename(columns={"pval":"P", "beta": "BETA"}, inplace=True)

select_cols = ["CHR", "SNP", "BP", "A1", "A2", "P", "BETA", "Z"]
df["SNP"] = df["CHR"] + ":" + df["BP"]
df_funct = df.loc[:,select_cols]
df_funct.to_csv("SBP_summarystats_baseline_ldpredFormat.txt", sep=" ", index=False)

#drop duplicates of summary stats:
df_uniq = df_funct.drop_duplicates(["SNP"])
df_uniq.to_csv("SBP_Unique_summarystats_baseline_ldpredFormat.txt", sep=" ", index=False)


#eQTL filtered PRS:
df_sbp = pd.read_csv("Cells_EBV-transformed_lymphocytes_prs_format_0.1.txt", sep="\t")
df_merge = pd.merge(df_sbp, df, left_on="ID", right_on="variant")
select_cols = ["CHR", "SNP", "BP", "A1", "A2", "P", "BETA", "Z"]

df_merge["SNP"] = df_merge["CHR"] + ":" + df_merge["BP"]
df_mergefunct = df_merge.loc[:,select_cols]
df_mergefunct.to_csv("SBP_eQTLprsformat_summarystats_baseline_ldpredFormat.txt", sep=" ", index=False)

#dummy per snp heritability baseline funct file:
df["h2SNP"] = 0
df_baseline = df.loc[:, ["SNP", "h2SNP"]]
df_baseline.to_csv("functfile_allsnp_baseline.txt", sep="\t", index=False)

#dummy per snp heritability eQTL funct file:
df_merge["h2SNP"] = 0
df_eQTLbaseline = df_merge.loc[:, ["SNP", "h2SNP"]]
df_eQTLbaseline.to_csv("functfile_eQTL_baseline.txt", sep="\t", index=False)


# slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y2)

#phenotype file format:
AA_file = "/work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/PAGE/race3/unrelated/pheno_mendrand.txt"
df_pheno = pd.read_csv(AA_file, sep="\t") 
df_phenoformat = df_pheno.loc[:,["IID", "SBP"]]
df_phenoformat.to_csv("phenotype_AA.txt", sep="\t", header=False, index=False)

EUR_file = "/work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/WHIMS/race5/unrelated/pheno_mendrand.txt"
df_pheno = pd.read_csv(EUR_file, sep="\t") 
df_phenoformat = df_pheno.loc[:,["IID", "SBP"]]
df_phenoformat.to_csv("phenotype_EUR.txt", sep="\t", header=False, index=False)


#remove null values from phenotype:
df_pheno = pd.read_csv("/work-zfs/abattle4/surya/datasets/WHIMS/phenotype_AA.txt", header=None, sep="\t")
df_pheno_final = df_pheno.dropna()
df_pheno_final.to_csv("phenotype_AA.txt", sep="\t", header=False, index=False)


