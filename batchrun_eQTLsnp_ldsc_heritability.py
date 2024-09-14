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
#work_dir = "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps"
#work_dir = "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/high_ldsnps"

annotpath = join(work_dir, "ldsc_annot")
prefix="EBV_LCL_Binary_100kb_colocthresh075_nonull"
prefix="EBV_LCL_Binary_100kb_colocthresh075_withBase"


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
    # annotation_list = ["PP.H0.7", "PP.H1.7", "PP.H2.7", "PP.H3.7", "PP.H4.7"]
    annotation_list = ["PP.H1.75", "PP.H2.75", "PP.H3.75", "PP.H4.75"]
    annotation_list = ["PP.H0.75", "PP.H1.75", "PP.H2.75", "PP.H3.75", "PP.H4.75"]
    
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

    #dummy intercept base variable
    df_merged["Base"] = 1

    df_merged.to_csv(args_annot_file, sep = "\t", index = False, compression="gzip")

print("Task completed...")


####################
##run s-LDSC program:
####################

scriptpath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc"
#scriptpath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc"

bfile_path="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_EUR_Phase3_plink"
#bfile_path="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_EUR_Phase3_plink"

hapmap3_snpspath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/hapmap3_snps"
#hapmap3_snpspath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/hapmap3_snps"

freqpath="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_Phase3_frq"
#freqpath="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_Phase3_frq"

regweights_path="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/1000G_Phase3_weights_hm3_no_MHC"
#regweights_path="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/ldscDB/1000G_Phase3_weights_hm3_no_MHC"

##################
#parameter changes

#work_dir="/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/high_ldsnps"
#work_dir="/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps"
work_dir="/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps"

annotpath="${work_dir}/ldsc_annot"

prefix="EBV_LCL_Binary_100kb_colocthresh075_nonull"
prefix="EBV_LCL_Binary_100kb_colocthresh075_withBase"

annotation_dir="${annotpath}/${prefix}"

genome_1kg_prefix="1000G.EUR.QC"

hm3snp_prefix="weights.hm3_noMHC"


#compute ldscore estimates
for chr in {1..22}; do \
    python ${scriptpath}/ldsc.py \
      --l2 \
      --bfile ${bfile_path}/1000G.EUR.QC.${chr} \
      --ld-wind-cm 1 \
      --print-snps ${hapmap3_snpspath}/hm.${chr}.snp \
      --annot ${annotation_dir}/${prefix}.${chr}.annot.gz \
      --out ${annotation_dir}/${prefix}.${chr}
done

sumstats_path="/Users/suryachhetri/datasets/prs_project/ldsc/ldsc/ldscDB/my_annot/EBV_LCL_Binary"
sumstats_path="/work-zfs/abattle4/surya/datasets/prs_project/ldsc/summary_stats_ldsc_format"
sumstats_path="/work-zfs/abattle4/surya/datasets/prs_project/UKBB/highly_heritable_traits/ldsc_sumstats"

for sumstatsfile in ${sumstats_path}/*.bgz; do \

    filename=$(basename ${sumstatsfile} | cut -f1 -d "."); \
    filename=$(basename ${sumstatsfile} | cut -f1-2 -d "."); \
    outputdir="${annotation_dir}/results_${filename}"; \
    outputdir="${annotation_dir}/highly_heritableTrait_results_${filename}"; \
    if [[ ! -d ${outputdir} ]];then mkdir -p ${outputdir}; fi

    #compute proportion of heritability explained by snps
    python ${scriptpath}/ldsc.py \
        --h2 ${sumstatsfile} \
        --ref-ld-chr ${annotation_dir}/${prefix}. \
        --frqfile-chr ${freqpath}/${genome_1kg_prefix}. \
        --w-ld-chr ${regweights_path}/${hm3snp_prefix}. \
        --overlap-annot \
        --print-coefficients \
        --print-delete-vals \
        --out ${outputdir}/${prefix}.ldsc; \

done



#####################
# R Concat dataframe

library(tidyverse)
library(data.table)

prefix <- "EBV_LCL_Binary_100kb_colocthresh075_nonull"

main_dir <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot/EBV_LCL_Binary_100kb_colocthresh075_nonull" 
output_dir <- "/work-zfs/abattle4/surya/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldsc_annot"
results_dirlist <- list.files(main_dir, "^results*", include.dirs = TRUE, full.names=TRUE)
results_dirlist <- list.files(main_dir, "^highly_heritableTrait_results*", include.dirs = TRUE, full.names=TRUE)
results_dirlist <- list.files(main_dir, "*results*", include.dirs = TRUE, full.names=TRUE)

results_list <- list() 

for(res_dir in results_dirlist) {
    filename <- unlist(strsplit(basename(res_dir), "results_"))[2]
    print(paste("Processing...", filename))
    result_file <- file.path(res_dir, paste0(prefix, ".ldsc.results"))
    result_df <- fread(result_file, sep="\t")
    result_df$trait <- filename
    results_list[[filename]] <- result_df
}

#concat dataframes
result_concat_df <- bind_rows(results_list) %>% data.frame #.id="id" for giving file_ids

#file names
outfile <- file.path(output_dir, paste0("combined_results_", prefix, ".txt"))
outfile1 <- file.path(output_dir, paste0("combined_highly_heritableTrait_results_", prefix, ".txt"))
outfile2 <- file.path(output_dir, paste0("combined_ALL_results_", prefix, ".txt"))

#write final df
fwrite(result_concat_df, outfile, sep="\t", row.names=F, quote=F, col.names=T)
fwrite(result_concat_df, outfile1, sep="\t", row.names=F, quote=F, col.names=T)
fwrite(result_concat_df, outfile2, sep="\t", row.names=F, quote=F, col.names=T)



############################################
#barplots of trait-SNP heritatbility results

library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(pheatmap)
library(ggh4x) #rbase like plots with truncated x-y axis
library(ggthemes) #rbase like boxed-background panel plots
library(ggdendro) #for clustering
library(scales)
library(RColorBrewer)

main_dir <- "/Users/suryachhetri/datasets/prs_project/ldsc"

#mark input files
all_result_file <- file.path(main_dir, "combined_ALL_results_EBV_LCL_Binary_100kb_colocthresh075_nonull.txt")
all_result_file <- file.path(main_dir, "combined_ALL_results_EBV_LCL_Binary_100kb_colocthresh075.txt")
herit_result_file <- file.path(main_dir, "combined_highly_heritableTrait_results_EBV_LCL_Binary_100kb_colocthresh075_nonull.txt")
herit_result_file <- file.path(main_dir, "combined_highly_heritableTrait_results_EBV_LCL_Binary_100kb_colocthresh075.txt")
result_file <- file.path(main_dir, "combined_results_EBV_LCL_Binary_100kb_colocthresh075_nonull.txt")
result_file <- file.path(main_dir, "combined_results_EBV_LCL_Binary_100kb_colocthresh075.txt")

#load input files
all_result_df <- fread(all_result_file, sep="\t")
all_result_df$Category <- gsub(".75L2_0", "", all_result_df$Category)
all_result_df$log10PValue <- -log10(all_result_df$Enrichment_p)

herit_result_df <- fread(herit_result_file, sep="\t")
herit_result_df$Category <- gsub(".75L2_0", "", herit_result_df$Category)
herit_result_df$log10PValue <- -log10(herit_result_df$Enrichment_p)

result_df <- fread(result_file, sep="\t")
result_df$Category <- gsub(".75L2_0", "", result_df$Category)
result_df$log10PValue <- -log10(result_df$Enrichment_p)

#filter based on traits
rbc_df <- result_df[grepl("rbc_bloodtrait", result_df$trait)]
wbc_df <- result_df[grepl("wbc_bloodtrait", result_df$trait)]
plt_df <- result_df[grepl("plt_bloodtrait", result_df$trait)]

bloodtrait_df <- result_df[grepl("bloodtrait", result_df$trait)] %>% filter(trait!="nucleated_erythrocyte_count_rbc_bloodtrait")
immunetrait_df <- result_df[grepl("immune", result_df$trait)]
cancertrait_df <- result_df[grepl("cancer", result_df$trait)]


#Bloodtrait: SNP-HERITABILITY PROPORTION
pdf(file.path(main_dir, "colocSnpHeritability_proportion_bloodtraits.pdf"))

ggplot(bloodtrait_df, aes(x = trait, y = Prop._h2, fill = Category)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +
    coord_flip()

dev.off()

#Bloodtrait: SNP-ENRICHMENT
pdf(file.path(main_dir, "colocSnpHeritability_enrichment_bloodtraits.pdf"), width=8, height=7)

ggplot(bloodtrait_df, aes(x = trait, y = Enrichment, fill = Category)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +
    coord_flip()

dev.off()

#Immunetrait: SNP-HERITABILITY PROPORTION
pdf(file.path(main_dir, "colocSnpHeritability_proportion_immunetraits.pdf"), width=8, height=7)

ggplot(immunetrait_df, aes(x = trait, y = Prop._h2, fill = Category)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +
    coord_flip()

dev.off()

#Immunetrait: SNP-ENRICHMENT
pdf(file.path(main_dir, "colocSnpHeritability_enrichment_immunetraits.pdf"), width=8, height=7)

ggplot(immunetrait_df, aes(x = trait, y = Enrichment, fill = Category)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +
    coord_flip()

dev.off()

#Cancertrait: SNP-HERITABILITY PROPORTION
pdf(file.path(main_dir, "colocSnpHeritability_proportion_Cancertraits.pdf"), width=8, height=7)

ggplot(cancertrait_df, aes(x = trait, y = Prop._h2, fill = Category)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +
    coord_flip()

dev.off()

#Cancertrait: SNP-ENRICHMENT
pdf(file.path(main_dir, "colocSnpHeritability_enrichment_Cancertraits.pdf"), width=8, height=7)

ggplot(cancertrait_df, aes(x = trait, y = Enrichment, fill = Category)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +
    coord_flip()

dev.off()

#Highly Heritable Trait: SNP-HERITABILITY PROPORTION
pdf(file.path(main_dir, "colocSnpHeritability_proportion_highlyHeritabletraits.pdf"), width=8, height=7)

herit_result_df1 <- herit_result_df %>% separate(trait, c("Trait", "phenoid"), sep="\\.")
ggplot(herit_result_df1, aes(x = Trait, y = Prop._h2, fill = Category)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() + 
    theme(text = element_text(size=7, face="bold")) +
    coord_flip() 
    # theme(
    #     axis.text.x = element_text(size=5, angle = 90),
    #     axis.text.y = element_text(size=9, angle = 90),
    #     plot.title = element_text(hjust =0.5)
    #     )
    

dev.off()


#Highly Heritable Trait: SNP-ENRICHMENT
pdf(file.path(main_dir, "colocSnpHeritability_enrichment_highlyHeritabletraits.pdf"), width=8, height=7)

herit_result_df2 <- herit_result_df %>% separate(trait, c("Trait", "phenoid"), sep="\\.")
ggplot(herit_result_df2, aes(x = Trait, y = Enrichment, fill = Category)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +
    theme(text = element_text(size=7, face="bold")) +
    coord_flip()

dev.off()


#########################################
#TOTAL result - all results combined
#HEATMAP OF SNP-Enrichment HERITABILITY

#cast data to wide format
read_df <- dcast(all_result_df, trait~Category, value.var="Enrichment")
data <- read_df[,2:ncol(read_df)]

#generate matrix
mat_data <- as.matrix(data)

#scale rowwise groupby traits:
#mat_data <- t(apply(result_mat, 1, scale))

#naming of cols and rows
colnames(mat_data) <- colnames(data)
rownames(mat_data) <- read_df$trait

#annotate traits
read_df$annotation <- "Other trait"
read_df[grepl("rbc_bloodtrait", read_df$trait),"annotation"] = "RBC Bloodtrait"
read_df[grepl("wbc_bloodtrait", read_df$trait),"annotation"] = "WBC Bloodtrait"
read_df[grepl("plt_bloodtrait", read_df$trait),"annotation"] = "PLT Bloodtrait"
read_df[grepl("immune", read_df$trait),"annotation"] = "Immune trait"
read_df[grepl("cancer", read_df$trait),"annotation"] = "Cancer trait"

#annotation prepare for heritable result
heritable_df <- herit_result_df %>% select(trait)
heritable_df$annotation <- "Heritable Trait"

#matched and unmatched values
inner_df <- inner_join(read_df, heritable_df, by="trait") %>% distinct(trait, .keep_all=TRUE)
anti_df <- anti_join(read_df, heritable_df, by="trait")

#bind rows and assign values
combined_df <- bind_rows(inner_df, anti_df)
combined_df[grepl("Heritable", combined_df$annotation.y), "annotation"] = "Heritable Trait"
final_df <- combined_df %>% select(-c(annotation.x, annotation.y))

#for complex cases using mutate dplyr
# regtags1 <- c("eQTL", "binding", "Minimal", "Other")
# final_df %>% 
#   mutate(
#     regulomtag1 = case_when(
#     rankscore=="1a" ~ regtags1[1],
#     rankscore=="2a" ~ regtags1[2],
#     rankscore=="3a" ~ regtags1[3])
#     )


# #file save
# pdf(file.path(main_dir, "colocSnpHeritability_enrichmentHeatmap.pdf"))

# #pheatmap plot:
# hmap1 <- pheatmap(mat_data, fontsize_col=9, fontsize_row=2.5)

# dev.off()


#####################
#rowannotated heatmap
metadata <- combined_df %>% select(annotation) %>% data.frame
rownames(metadata) <- combined_df$trait

#rearrange of levels
re_annotation <- c("Cancer trait", "Heritable Trait", "Immune trait", 
                "RBC Bloodtrait", "WBC Bloodtrait", "PLT Bloodtrait", 
                "Other trait")

#levels refactor
metadata$annotation <- with(metadata, factor(annotation, re_annotation))

#custom color
colornew = colorRampPalette(c("yellow", "red"))(100)
color = colorRampPalette(
        rev(brewer.pal(n = 7, name ="RdYlBu"))
        )(100)

#file save
pdf(file.path(main_dir, "colocSnpHeritability_enrichmentAnnotHeatmap3.pdf"))

#pheatmap plot:
#colornew = colorRampPalette(c("yellow", "white", "red"))(100)
hmap2 <- pheatmap(mat_data, color=colornew, annotation_row=metadata, fontsize_col=9, fontsize_row=2, annotation_legend=TRUE)

dev.off()


#########################################
#TOTAL result - all results combined
#HEATMAP OF SNP-PValue HERITABILITY

#cast data to wide format
read_df <- dcast(all_result_df, trait~Category, value.var="log10PValue")
data <- read_df[,2:ncol(read_df)]

#generate matrix
mat_data <- as.matrix(data)

#scale rowwise groupby traits:
#mat_data <- t(apply(result_mat, 1, scale))

#naming of cols and rows
colnames(mat_data) <- colnames(data)
rownames(mat_data) <- read_df$trait

#annotate traits
read_df$annotation <- "Other trait"
read_df[grepl("rbc_bloodtrait", read_df$trait),"annotation"] = "RBC Bloodtrait"
read_df[grepl("wbc_bloodtrait", read_df$trait),"annotation"] = "WBC Bloodtrait"
read_df[grepl("plt_bloodtrait", read_df$trait),"annotation"] = "PLT Bloodtrait"
read_df[grepl("immune", read_df$trait),"annotation"] = "Immune trait"
read_df[grepl("cancer", read_df$trait),"annotation"] = "Cancer trait"

#annotation prepare for heritable result
heritable_df <- herit_result_df %>% select(trait)
heritable_df$annotation <- "Heritable Trait"

#matched and unmatched values
inner_df <- inner_join(read_df, heritable_df, by="trait") %>% distinct(trait, .keep_all=TRUE)
anti_df <- anti_join(read_df, heritable_df, by="trait")

#bind rows and assign values
combined_df <- bind_rows(inner_df, anti_df)
combined_df[grepl("Heritable", combined_df$annotation.y), "annotation"] = "Heritable Trait"
final_df <- combined_df %>% select(-c(annotation.x, annotation.y))


#####################
#rowannotated heatmap
metadata <- combined_df %>% select(annotation) %>% data.frame
rownames(metadata) <- combined_df$trait

#rearrange of levels
re_annotation <- c("Cancer trait", "Heritable Trait", "Immune trait", 
                "RBC Bloodtrait", "WBC Bloodtrait", "PLT Bloodtrait", 
                "Other trait")

#levels refactor
metadata$annotation <- with(metadata, factor(annotation, re_annotation))


#breaks at the quantiles of the data, 
#then each color will represent an equal proportion of the data:

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(mat_data, n = 100)

#file save
pdf(file.path(main_dir, "colocSnpHeritability_PValueAnnotHeatmap.pdf"))

#pheatmap plot:
hmap2 <- pheatmap(mat_data,annotation_row=metadata, 
          breaks = mat_breaks, fontsize_col=9, fontsize_row=2, annotation_legend=TRUE)

dev.off()

#########################################\
#HEATMAP OF SNP-PValue HERITABILITY for subset of dataframe

#cast data to wide format
read_df <- dcast(all_result_df, trait~Category, value.var="log10PValue")


#scale rowwise groupby traits:
#mat_data <- t(apply(result_mat, 1, scale))

#annotate traits
read_df$annotation <- "Other trait"
read_df[grepl("rbc_bloodtrait", read_df$trait),"annotation"] = "RBC Bloodtrait"
read_df[grepl("wbc_bloodtrait", read_df$trait),"annotation"] = "WBC Bloodtrait"
read_df[grepl("plt_bloodtrait", read_df$trait),"annotation"] = "PLT Bloodtrait"
read_df[grepl("immune", read_df$trait),"annotation"] = "Immune trait"
read_df[grepl("cancer", read_df$trait),"annotation"] = "Cancer trait"

#subset only annotated data
final_df <- read_df[read_df$annotation != "Other trait",]

#trim numeric data
data <- final_df[,2:(ncol(final_df)-1)]

#generate matrix
mat_data <- as.matrix(data)

#naming of cols and rows
colnames(mat_data) <- colnames(data)
rownames(mat_data) <- final_df$trait

#####################
#rowannotated heatmap
metadata <- final_df %>% select(annotation) %>% data.frame
rownames(metadata) <- final_df$trait

#rearrange of levels
re_annotation <- c("Cancer trait", "Immune trait", 
                "RBC Bloodtrait", "WBC Bloodtrait", "PLT Bloodtrait")

#levels refactor
metadata$annotation <- with(metadata, factor(annotation, re_annotation))


#breaks at the quantiles of the data, 
#then each color will represent an equal proportion of the data:

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(mat_data, n = 100)

#file save
pdf(file.path(main_dir, "colocSnpHeritability_PValueAnnotSubsetHeatmap1.pdf"))

#pheatmap plot:
hmap2 <- pheatmap(mat_data,annotation_row=metadata, 
          breaks = mat_breaks, fontsize_col=9, fontsize_row=5, annotation_legend=TRUE)

dev.off()

#########################################\
#HEATMAP OF SNP-Enrichment HERITABILITY for subset of dataframe

#cast data to wide format
read_df <- dcast(all_result_df, trait~Category, value.var="Enrichment")


#annotate traits
read_df$annotation <- "Other trait"
read_df[grepl("rbc_bloodtrait", read_df$trait),"annotation"] = "RBC Bloodtrait"
read_df[grepl("wbc_bloodtrait", read_df$trait),"annotation"] = "WBC Bloodtrait"
read_df[grepl("plt_bloodtrait", read_df$trait),"annotation"] = "PLT Bloodtrait"
read_df[grepl("immune", read_df$trait),"annotation"] = "Immune trait"
read_df[grepl("cancer", read_df$trait),"annotation"] = "Cancer trait"

#subset only annotated data
final_df <- read_df[read_df$annotation != "Other trait",]

#trim numeric data
data <- final_df[,2:(ncol(final_df)-1)]

#generate matrix
mat_data <- as.matrix(data)

#scale rowwise groupby traits:
mat_data <- t(apply(mat_data, 1, scale))

#naming of cols and rows
colnames(mat_data) <- colnames(data)
rownames(mat_data) <- final_df$trait

#####################
#rowannotated heatmap
metadata <- final_df %>% select(annotation) %>% data.frame
rownames(metadata) <- final_df$trait

#rearrange of levels
re_annotation <- c("Cancer trait", "Immune trait", 
                "RBC Bloodtrait", "WBC Bloodtrait", "PLT Bloodtrait")

#levels refactor
metadata$annotation <- with(metadata, factor(annotation, re_annotation))

#custom color
colornew = colorRampPalette(c("yellow", "red"))(100)
color = colorRampPalette(
        rev(brewer.pal(n = 7, name ="RdYlBu"))
        )(100)

#file save
pdf(file.path(main_dir, "colocSnpHeritability_EnrichmentAnnotSubsetHeatmap2.pdf"))

#pheatmap plot:
hmap2 <- pheatmap(mat_data, color=color, annotation_row=metadata, 
          fontsize_col=9, fontsize_row=5, annotation_legend=FALSE)

dev.off()



#########################################
#########################################

# New BASELINE LD based analysis

#########################################
#########################################

library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(pheatmap)
library(ggh4x) #rbase like plots with truncated x-y axis
library(ggthemes) #rbase like boxed-background panel plots
library(ggdendro) #for clustering
library(scales)
library(RColorBrewer)

main_dir <- "/Users/suryachhetri/datasets/prs_project/ldsc"
main_dir <- "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldscrun_combined_results"

#mark input files
all_result_file <- file.path(main_dir, "combined_ALL_results_EBV_LCL_Binary_100kb_colocthresh075.txt")
all_result_file <- file.path(main_dir, "combined_ALL_results_EBV_LCL_Binary_100kb_colocthresh075_baselineLD.txt")

#load input files
all_result_df <- fread(all_result_file, sep="\t")
all_result_df$Category <- gsub(".75L2_0", "", all_result_df$Category)
all_result_df$Category <- gsub("L2_0", "", all_result_df$Category)
all_result_df$log10PValue <- -log10(all_result_df$Enrichment_p)

#annotate traits
annotation_file <- file.path(main_dir, "Combined_trait_results.annotation.txt")
annot_df <- fread(annotation_file, sep="\t")

#merge the annotation
final_df_left <- left_join(all_result_df, annot_df, by="trait") %>% distinct(trait, .keep_all=TRUE)
final_df <- inner_join(all_result_df, annot_df, by="trait")

# #filter the traits
# annot_df %>% filter(grepl("rbc_bloodtrait", trait))
# annot_df %>% filter(grepl("wbc_bloodtrait", trait))
# annot_df %>% filter(grepl("plt_bloodtrait", trait))

annot_df %>% filter(annotation == "Anthropological trait")


rbc_trait <- final_df %>% filter(grepl("rbc_bloodtrait", trait))
wbc_trait <- final_df %>% filter(grepl("wbc_bloodtrait", trait))
plt_trait <- final_df %>% filter(grepl("plt_bloodtrait", trait))
immune_trait <- final_df %>% filter(grepl("ulcerative_collitis_ICDK51.K51", trait))
immune_trait <- final_df %>% filter(grepl("ulcerative_collitis.ULCERNAS", trait))

anthro_trait <- final_df %>% filter(grepl("Heel_BMD_left.4106_irnt", trait))
anthro_trait <- final_df %>% filter(grepl("Heel_BMD_right.4125_irnt", trait))

#find distinct traits
rbc_trait %>% distinct(trait)
wbc_trait %>% distinct(trait)
plt_trait %>% distinct(trait)
immune_trait %>% distinct(trait)

erythro <- rbc_trait %>% filter(trait == "erythrocyte_count_rbc_bloodtrait.30010_irnt")
erythro %>% arrange(desc(Enrichment))

################
#rbc trait
erythro1 <- rbc_trait %>% filter(trait == "haemoglobin_conc_rbc_bloodtrait.30020_irnt")
erythro_sort <- erythro1 %>% arrange(desc(Enrichment))
erythro_filt <- erythro_sort %>% filter(Enrichment > 0)

erythro2 <- rbc_trait %>% filter(trait == "mchc_rbc_bloodtrait.30060_irnt")
erythro2 %>% arrange(desc(Enrichment))

erythro3 <- rbc_trait %>% filter(trait == "mcv_rbc_bloodtrait.30040_irnt")
erythro3 %>% arrange(desc(Enrichment))

################
#plt trait
plt_df <- plt_trait %>% filter(trait == "pdw_plt_bloodtrait.30110_irnt")
plt_df %>% arrange(desc(Enrichment))

###############
#wbc trait
wbc_df <- wbc_trait %>% filter(trait == "baso_count_wbc_bloodtrait.30160")
wbc_df %>% arrange(desc(Enrichment))
> 1.295233e-06/2.061816e-07
[1] 6.282001

immune_df <- immune_trait %>% filter(trait == "ulcerative_collitis_ICDK51.K51")

##################
#immune trait
immune_df <- immune_trait %>% filter(trait == "ulcerative_collitis.ULCERNAS")
immune_df %>% arrange(desc(Enrichment))
#highlight PP.H0
4.948604e-07/6.021395e-08
8.218

#anthro trait
anthro_df <- anthro_trait %>% filter(trait == "Heel_BMD_left.4106_irnt")
anthro_df <- anthro_trait %>% filter(trait == "Heel_BMD_right.4125_irnt")
anthro_df %>% arrange(desc(Enrichment))


#annotate traits
read_df$annotation <- "Other trait"
read_df[grepl("rbc_bloodtrait", read_df$trait),"annotation"] = "RBC Bloodtrait"
read_df[grepl("wbc_bloodtrait", read_df$trait),"annotation"] = "WBC Bloodtrait"
read_df[grepl("plt_bloodtrait", read_df$trait),"annotation"] = "PLT Bloodtrait"
read_df[grepl("immune", read_df$trait),"annotation"] = "Immune trait"
read_df[grepl("cancer", read_df$trait),"annotation"] = "Cancer trait"

#Barplot Trait: SNP-ENRICHMENT
pdf(file.path(main_dir, "colocSnpHeritability_enrichment_102trait-bar.pdf"), width=8, height=7)




################
#################
#################


#rbc trait
erythro1 <- rbc_trait %>% filter(trait == "haemoglobin_conc_rbc_bloodtrait.30020_irnt")
erythro_sort <- erythro1 %>% arrange(desc(Enrichment))
erythro_filt <- erythro_sort %>% filter(Enrichment > 0)

#herit_result_df2 <- herit_result_df %>% separate(trait, c("Trait", "phenoid"), sep="\\.")
ggplot(erythro_filt, aes(x = reorder(Category, -Enrichment), y = Enrichment)) + 
    geom_bar(stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    #scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +  coord_flip() +
    scale_x_discrete(limits=rev) +

    # +
    theme(#text = element_text(size=6, face="bold"),
        #axis.text.x = element_text(size=8, face="bold"),
        axis.text.y = element_text(size=4, face="bold")
        ) +
    labs(x="Annotation Category")

#plt trait
wbc1 <- wbc_trait %>% filter(trait == "pdw_wbc_bloodtrait.30110_irnt")
wbc_sort <- wbc1 %>% arrange(desc(Enrichment))
wbc_filt <- wbc_sort %>% filter(Enrichment > 0)

#herit_result_df2 <- herit_result_df %>% separate(trait, c("Trait", "phenoid"), sep="\\.")
ggplot(plt_filt, aes(x = reorder(Category, -Enrichment), y = Enrichment)) + 
    geom_bar(stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    #scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +  coord_flip() +
    scale_x_discrete(limits=rev) +

    # +
    theme(#text = element_text(size=6, face="bold"),
        #axis.text.x = element_text(size=8, face="bold"),
        axis.text.y = element_text(size=4, face="bold")
        ) +
    labs(x="Annotation Category")

#wbc trait
wbc1 <- wbc_trait %>% filter(trait == "baso_count_wbc_bloodtrait.30160")
wbc_sort <- wbc1 %>% arrange(desc(Enrichment))
wbc_filt <- wbc_sort %>% filter(Enrichment > 0)

#herit_result_df2 <- herit_result_df %>% separate(trait, c("Trait", "phenoid"), sep="\\.")
ggplot(wbc_filt, aes(x = reorder(Category, -Enrichment), y = Enrichment)) + 
    geom_bar(stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    #scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +  coord_flip() +
    scale_x_discrete(limits=rev) +

    # +
    theme(#text = element_text(size=6, face="bold"),
        #axis.text.x = element_text(size=8, face="bold"),
        axis.text.y = element_text(size=4, face="bold")
        ) +
    labs(x="Annotation Category")


#immune trait
immune_trait <- final_df %>% filter(grepl("ulcerative_collitis.ULCERNAS", trait))
immune1 <- immune_trait %>% filter(trait == "ulcerative_collitis.ULCERNAS")
immune_sort <- immune1 %>% arrange(desc(Enrichment))
immune_filt <- immune_sort %>% filter(Enrichment > 0)

#herit_result_df2 <- herit_result_df %>% separate(trait, c("Trait", "phenoid"), sep="\\.")
ggplot(immune_filt, aes(x = reorder(Category, -Enrichment), y = Enrichment)) + 
    geom_bar(stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    #scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +  coord_flip() +
    scale_x_discrete(limits=rev) +

    # +
    theme(#text = element_text(size=6, face="bold"),
        #axis.text.x = element_text(size=8, face="bold"),
        axis.text.y = element_text(size=4, face="bold")
        ) +
    labs(x="Annotation Category")

#anthro trait
anthro_trait <- final_df %>% filter(grepl("Heel_BMD_left.4106_irnt", trait))
anthro_trait <- final_df %>% filter(grepl("Heel_BMD_right.4125_irnt", trait))

anthro1 <- anthro_trait %>% filter(trait == "Heel_BMD_left.4106_irnt")
anthro1 <- anthro_trait %>% filter(trait == "Heel_BMD_right.4125_irnt")
anthro_sort <- anthro1 %>% arrange(desc(Enrichment))
anthro_filt <- anthro_sort %>% filter(Enrichment > 0)

#herit_result_df2 <- herit_result_df %>% separate(trait, c("Trait", "phenoid"), sep="\\.")
ggplot(anthro_filt, aes(x = reorder(Category, -Enrichment), y = Enrichment)) + 
    geom_bar(stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    #scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +  coord_flip() +
    scale_x_discrete(limits=rev) +

    # +
    theme(#text = element_text(size=6, face="bold"),
        #axis.text.x = element_text(size=8, face="bold"),
        axis.text.y = element_text(size=4, face="bold")
        ) +
    labs(x="Annotation Category")



    # #load file
    # read_df <- fread(input_file, header=T) %>% 
    #   as.data.frame() %>%
    #   separate(Bed_File, c("TF", "txt"),  sep="[.]")

    # #indexing by TFs itself
    # rownames(read_df) <- read_df$TF

    # #fold enrich i.e inbed_hits/expected_hits :
    # #compute enrichment values (observed/expected TF-SNP matches)
    # enrich <- read_df[,3]/read_df[,4]
    # read_df$enrich <- enrich
    # read_df$log10PValue <- -log10(read_df$PValue)

    # #order by highest -log10 Pvalue enrichment first
    # read_df <- read_df[order(read_df[,7], decreasing=T),] 
    # read_df <- read_df[!is.na(read_df$PValue),]     #drop na

plt1 <- plt_trait %>% filter(trait == "pdw_plt_bloodtrait.30110_irnt")
plt_sort <- plt1 %>% arrange(desc(Coefficient))
plt_filt <- plt_sort %>% filter(Coefficient > 0)

read_df <- plt_sort %>% mutate(Tau = Coefficient*1000000)
read_df <- read_df %>% filter(Enrichment >0)

label_data <- read_df %>% arrange(desc(Enrichment)) %>% head(5)

wbc1 <- wbc_trait %>% filter(trait == "baso_count_wbc_bloodtrait.30160")
wbc_sort <- wbc1 %>% arrange(desc(Coefficient))
wbc_filt <- wbc_sort %>% filter(Coefficient > 0)

read_df <- wbc_sort %>% mutate(Tau = Coefficient*1000000)
read_df <- read_df %>% filter(Enrichment >0)

rbc1 <- rbc_trait %>% filter(trait == "haemoglobin_conc_rbc_bloodtrait.30020_irnt")
rbc_sort <- rbc1 %>% arrange(desc(Coefficient))
rbc_filt <- rbc_sort %>% filter(Coefficient > 0)

read_df <- rbc_sort %>% mutate(Tau = Coefficient*1000000)
read_df <- read_df %>% filter(Enrichment >0)


immune1 <- immune_trait %>% filter(trait == "ulcerative_collitis.ULCERNAS")
immune_sort <- immune1 %>% arrange(desc(Coefficient))
immune_filt <- immune_sort %>% filter(Coefficient > 0)

read_df <- immune_sort %>% mutate(Tau = Coefficient*1000000)
read_df <- read_df %>% filter(Enrichment >0)

label_data <- read_df %>% arrange(desc(Enrichment)) %>% head(8)

#library(ggrepel)
options(ggrepel.max.overlaps = 100)

# #label data points
# label_data1 <- read_final_df[rev(order(read_final_df$log10BonfPValue)),][c(1:30),]
# label_data2 <- read_final_df[rev(order(read_final_df$enrich)),][c(1:10),]

# label_data <- bind_rows(label_data1, label_data2)
# circle_point_data <- bind_rows(label_data1, label_data2)

# #remove all duplicate circles
# label_data <- label_data %>% distinct(TF, .keep_all=T)
# circle_point_data <- circle_point_data %>% distinct(TF, .keep_all=T)

# filename <- gsub("_Sumstats.txt", "_TF_enrichplot.pdf", infile)
# outfile <- file.path(output_dir, filename)

# #pdf(outfile, width=7, height=7)
# print(paste("processing", filename))
# cat("\n")

read_df %>% 
    ggplot(aes(x=Tau, y=Enrichment)) +
    geom_point(color="red") + #geom_line() +
    #geom_vline(colour="dodgerblue", xintercept=1, linetype=2) +
    #geom_hline(colour="dodgerblue", yintercept=2, linetype=2) + 
    #labs(y= "-Log10(Adj PValue)", x="Observed/Expected Fold-Enrichment", title="Transcription Factor Binding Sites") +
    theme_bw() + 

    theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)) +

    geom_text_repel(data = label_data,
                aes(label = Category), 
                box.padding = 1,
                #max.overlaps =10,
                show.legend = FALSE, #this removes the 'a' from the legend
                size=2.1
                ) 
    #+

    #geom_point(data=circle_point_data,
    #         pch=21, fill=NA, size=3, colour="dodgerblue", stroke=0.5)


dev.off()


#########################################
#########################################
#########################################
#########################################

#Cleaned up total results (Selective)
#filter out traits with negative Prop._h2
#Full 95 traits
#########################################
#########################################
#########################################
#########################################

library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(pheatmap)
library(ggh4x) #rbase like plots with truncated x-y axis
library(ggthemes) #rbase like boxed-background panel plots
library(ggdendro) #for clustering
library(scales)
library(RColorBrewer)

main_dir <- "/Users/suryachhetri/datasets/prs_project/ldsc"
main_dir <- "/Users/suryachhetri/datasets/prs_project/genoa/coloc_genoa_eQTL/coloc_topsnps/ldscrun_combined_results"

#mark input files
#all_result_file <- file.path(main_dir, "combined_ALL_results_EBV_LCL_Binary_100kb_colocthresh075_nonull.txt")
all_result_file <- file.path(main_dir, "combined_ALL_results_EBV_LCL_Binary_100kb_colocthresh075.txt")
all_result_file <- file.path(main_dir, "combined_ALL_results_EBV_LCL_Binary_100kb_colocthresh075_baselineLD.txt")

#herit_result_file <- file.path(main_dir, "combined_highly_heritableTrait_results_EBV_LCL_Binary_100kb_colocthresh075_nonull.txt")
herit_result_file <- file.path(main_dir, "combined_highly_heritableTrait_results_EBV_LCL_Binary_100kb_colocthresh075.txt")

#result_file <- file.path(main_dir, "combined_results_EBV_LCL_Binary_100kb_colocthresh075_nonull.txt")
result_file <- file.path(main_dir, "combined_results_EBV_LCL_Binary_100kb_colocthresh075.txt")

#load input files
all_result_df <- fread(all_result_file, sep="\t")
all_result_df$Category <- gsub(".75L2_0", "", all_result_df$Category)
all_result_df$Category <- gsub("L2_0", "", all_result_df$Category)
all_result_df$log10PValue <- -log10(all_result_df$Enrichment_p)

herit_result_df <- fread(herit_result_file, sep="\t")
herit_result_df$Category <- gsub(".75L2_0", "", herit_result_df$Category)
herit_result_df$log10PValue <- -log10(herit_result_df$Enrichment_p)

result_df <- fread(result_file, sep="\t")
result_df$Category <- gsub(".75L2_0", "", result_df$Category)
result_df$log10PValue <- -log10(result_df$Enrichment_p)


#########################################
#HEATMAP OF SNP-PValue HERITABILITY for subset of dataframe

#gather trait list with negative Prop._h2 
filter_trait_list <- all_result_df %>% filter(Prop._h2 < 0) %>% distinct(trait) %>% unlist(use.names=F)

#filter out traits with negative Prop._h2
all_clean_result_df <- all_result_df %>% filter(!trait %in% filter_trait_list) %>% data.table

#cast data to wide format
read_df <- dcast(all_clean_result_df, trait~Category, value.var="log10PValue")

#annotate traits
annotation_file <- file.path(main_dir, "Combined_trait_results.annotation.txt")
annot_df <- fread(annotation_file, sep="\t")

#matched and unmatched values
final_df <- left_join(read_df, annot_df, by="trait") %>% distinct(trait, .keep_all=TRUE)
#final_df <- inner_df %>% select(-c(annotation.x)) %>% rename("annotation"="annotation.y")

#trim numeric data
data <- final_df[,2:(ncol(final_df)-1)]

#generate matrix
mat_data <- as.matrix(data)

#scale rowwise groupby traits:
#mat_data <- t(apply(mat_data, 1, scale))

#naming of cols and rows
colnames(mat_data) <- colnames(data)
rownames(mat_data) <- final_df$trait

#####################
#rowannotated heatmap
metadata <- final_df %>% select(annotation) %>% data.frame
rownames(metadata) <- final_df$trait

#rearrange of levels
re_annotation <- c("Anthropological trait", "Immune trait", "Cancer trait", "Neurological trait",
                "RBC Bloodtrait", "WBC Bloodtrait", "PLT Bloodtrait", "BP Bloodtrait",
                "Cardio trait", "Lipid trait", "Skin trait", "Other trait")

#levels refactor
metadata$annotation <- with(metadata, factor(annotation, re_annotation))

#breaks at the quantiles of the data, 
#then each color will represent an equal proportion of the data:

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(mat_data, n = 100)

#file save
pdf(file.path(main_dir, "Cleaned_colocSnpHeritability_PValueAnnotSubsetHeatmap1.pdf"), height=8, width=9)

#pheatmap plot:
hmap2 <- pheatmap(mat_data,annotation_row=metadata, cluster_cols=FALSE,
          breaks = mat_breaks, fontsize_col=9, fontsize_row=4, annotation_legend=TRUE)

dev.off()


#########################################\
#HEATMAP OF SNP-Enrichment HERITABILITY for subset of dataframe

#gather trait list with negative Prop._h2 
filter_trait_list <- all_result_df %>% filter(Prop._h2 < 0) %>% distinct(trait) %>% unlist(use.names=F)

#filter out traits with negative Prop._h2
all_clean_result_df <- all_result_df %>% filter(!trait %in% filter_trait_list) %>% data.table

#cast data to wide format
read_df <- dcast(all_clean_result_df, trait~Category, value.var="Enrichment")

#annotate traits
annotation_file <- file.path(main_dir, "Combined_trait_results.annotation.txt")
annot_df <- fread(annotation_file, sep="\t")

#matched and unmatched values
final_df <- inner_join(read_df, annot_df, by="trait") %>% distinct(trait, .keep_all=TRUE)
#final_df <- inner_df %>% select(-c(annotation.x)) %>% rename("annotation"="annotation.y")


#trim numeric data
data <- final_df[,2:(ncol(final_df)-1)]

#generate matrix
mat_data <- as.matrix(data)

#scale rowwise groupby traits:
mat_data <- t(apply(mat_data, 1, scale))

#naming of cols and rows
colnames(mat_data) <- colnames(data)
rownames(mat_data) <- final_df$trait

#####################
#rowannotated heatmap
metadata <- final_df %>% select(annotation) %>% data.frame
rownames(metadata) <- final_df$trait

#rearrange of levels
re_annotation <- c("Anthropological trait", 
                "RBC Bloodtrait", "WBC Bloodtrait", "Immune trait", "PLT Bloodtrait",
                "Cancer trait", "Lipid trait", "Other trait")

#levels refactor
metadata$annotation <- with(metadata, factor(annotation, re_annotation))

#breaks at the quantiles of the data, 
#then each color will represent an equal proportion of the data:

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(mat_data, n = 100)

#custom color
colornew = colorRampPalette(c("yellow", "red"))(100)
color = colorRampPalette(
        rev(brewer.pal(n = 7, name ="RdYlBu"))
        )(100)

#file save
pdf(file.path(main_dir, "Cleaned_colocSnpHeritability_EnrichmentAnnotSubsetHeatmap_94trait.pdf"), height=8, width=9)

#pheatmap plot:
hmap2 <- pheatmap(mat_data, cluster_col=FALSE, color=colornew, annotation_row=metadata, 
          fontsize_col=9, fontsize_row=3.5, annotation_legend=TRUE)

dev.off()

#Barplot Trait: SNP-ENRICHMENT
pdf(file.path(main_dir, "colocSnpHeritability_enrichment_94trait-bar.pdf"), width=8, height=7)

#herit_result_df2 <- herit_result_df %>% separate(trait, c("Trait", "phenoid"), sep="\\.")
ggplot(all_clean_result_df, aes(x = trait, y = Enrichment, fill = Category)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +
    theme(text = element_text(size=6, face="bold"),
        axis.text.x = element_text(size=8, face="bold"),
        axis.text.y = element_text(size=4, face="bold")
        ) +
    coord_flip()

dev.off()


#######################################
#Trial = Full 102 traits
#######################################

main_dir <- "/Users/suryachhetri/datasets/prs_project/ldsc"

#mark input files
all_result_file <- file.path(main_dir, "combined_ALL_results_EBV_LCL_Binary_100kb_colocthresh075_nonull.txt")
#all_result_file <- file.path(main_dir, "combined_ALL_results_EBV_LCL_Binary_100kb_colocthresh075.txt")

herit_result_file <- file.path(main_dir, "combined_highly_heritableTrait_results_EBV_LCL_Binary_100kb_colocthresh075_nonull.txt")
#herit_result_file <- file.path(main_dir, "combined_highly_heritableTrait_results_EBV_LCL_Binary_100kb_colocthresh075.txt")

result_file <- file.path(main_dir, "combined_results_EBV_LCL_Binary_100kb_colocthresh075_nonull.txt")
#result_file <- file.path(main_dir, "combined_results_EBV_LCL_Binary_100kb_colocthresh075.txt")

#load input files
all_result_df <- fread(all_result_file, sep="\t")
all_result_df$Category <- gsub(".75L2_0", "", all_result_df$Category)
all_result_df$log10PValue <- -log10(all_result_df$Enrichment_p)

herit_result_df <- fread(herit_result_file, sep="\t")
herit_result_df$Category <- gsub(".75L2_0", "", herit_result_df$Category)
herit_result_df$log10PValue <- -log10(herit_result_df$Enrichment_p)

result_df <- fread(result_file, sep="\t")
result_df$Category <- gsub(".75L2_0", "", result_df$Category)
result_df$log10PValue <- -log10(result_df$Enrichment_p)


#########################################
#HEATMAP OF SNP-PValue HERITABILITY for subset of dataframe

#gather trait list with negative Prop._h2 
filter_trait_list <- all_result_df %>% filter(Prop._h2 < 0) %>% distinct(trait) %>% unlist(use.names=F)

#filter out traits with negative Prop._h2
all_clean_result_df <- all_result_df %>% filter(!trait %in% filter_trait_list) %>% data.table


#cast data to wide format
read_df <- dcast(all_clean_result_df, trait~Category, value.var="Enrichment")
read_df <- dcast(all_clean_result_df, trait~Category, value.var="Coefficient_z.score")


# #annotate traits
# read_df$annotation <- "Other trait"
# read_df[grepl("rbc_bloodtrait", read_df$trait),"annotation"] = "RBC Bloodtrait"
# read_df[grepl("wbc_bloodtrait", read_df$trait),"annotation"] = "WBC Bloodtrait"
# read_df[grepl("plt_bloodtrait", read_df$trait),"annotation"] = "PLT Bloodtrait"
# read_df[grepl("immune", read_df$trait),"annotation"] = "Immune trait"
# read_df[grepl("cancer", read_df$trait),"annotation"] = "Cancer trait"


# #annotation prepare for heritable result
# heritable_df <- herit_result_df %>% select(trait)
# heritable_df$annotation <- "Heritable Trait"

# #matched and unmatched values
# inner_df <- inner_join(read_df, heritable_df, by="trait") %>% distinct(trait, .keep_all=TRUE)
# anti_df <- anti_join(read_df, heritable_df, by="trait")

# #bind rows and assign values
# combined_df <- bind_rows(inner_df, anti_df)
# combined_df[grepl("Heritable", combined_df$annotation.y), "annotation"] = "Heritable Trait"
# final_df <- combined_df %>% select(-c(annotation.x, annotation.y))

# final_df %>% select(trait, annotation) %>% fwrite(.,file.path(main_dir, "Combined_trait_results.fullannotation.txt"), sep="\t", col.names=TRUE, row.names=FALSE)

#edit the #"Combined_trait_results.fullannotation.txt" files manually
#annotate traits
annotation_file <- file.path(main_dir, "Combined_trait_results.fullannotation.txt")
annot_df <- fread(annotation_file, sep="\t")

#matched and unmatched values
final_df <- inner_join(read_df, annot_df, by="trait") %>% distinct(trait, .keep_all=TRUE)
#final_df <- inner_df %>% select(-c(annotation.x)) %>% rename("annotation"="annotation.y")

#trim numeric data
data <- final_df[,2:(ncol(final_df)-1)]

#generate matrix
mat_data <- as.matrix(data)

#scale rowwise groupby traits:
mat_data <- t(apply(mat_data, 1, scale))

#naming of cols and rows
colnames(mat_data) <- colnames(data)
rownames(mat_data) <- final_df$trait

#####################
#rowannotated heatmap
metadata <- final_df %>% select(annotation) %>% data.frame
rownames(metadata) <- final_df$trait

#rearrange of levels
re_annotation <- c("Anthropological trait", 
                "RBC Bloodtrait", "WBC Bloodtrait", "Immune trait", "PLT Bloodtrait",
                "Cancer trait", "Lipid trait", "Other trait")

#levels refactor
metadata$annotation <- with(metadata, factor(annotation, re_annotation))

#breaks at the quantiles of the data, 
#then each color will represent an equal proportion of the data:

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(mat_data, n = 100)

#custom color
colornew = colorRampPalette(c("yellow", "red"))(100)
color = colorRampPalette(
        rev(brewer.pal(n = 7, name ="RdYlBu"))
        )(100)

#file save
pdf(file.path(main_dir, "Cleaned_colocSnpHeritability_EnrichmentAnnotSubsetHeatmap_102trait-Z.pdf"), height=8, width=9)

#pheatmap plot:
hmap2 <- pheatmap(mat_data, cluster_col=TRUE, color=colornew, annotation_row=metadata, 
          breaks=mat_breaks, fontsize_col=9, fontsize_row=3.5, annotation_legend=TRUE)

dev.off()


#Barplot Trait: SNP-ENRICHMENT
pdf(file.path(main_dir, "colocSnpHeritability_enrichment_102trait-bar.pdf"), width=8, height=7)

#herit_result_df2 <- herit_result_df %>% separate(trait, c("Trait", "phenoid"), sep="\\.")
ggplot(all_clean_result_df, aes(x = trait, y = Enrichment, fill = Category)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +
    theme(text = element_text(size=6, face="bold"),
        axis.text.x = element_text(size=8, face="bold"),
        axis.text.y = element_text(size=4, face="bold")
        ) +
    coord_flip()

dev.off()


result_melt_df <- melt(final_df, id.vars=c("annotation", "trait"))
#Barplot Trait: SNP-ENRICHMENT
pdf(file.path(main_dir, "bloodtrait_colocSnpHeritability_enrichment_102trait-bar.pdf"), width=8, height=7)

#herit_result_df2 <- herit_result_df %>% separate(trait, c("Trait", "phenoid"), sep="\\.")
ggplot(all_clean_result_df, aes(x = trait, y = Enrichment, fill = Category)) + 
    geom_bar(position = "fill",stat = "identity") +
    # or:
    # geom_bar(position = position_fill(), stat = "identity") 
    scale_y_continuous(labels = scales::percent_format()) +
    theme_bw() +
    theme(text = element_text(size=6, face="bold"),
        axis.text.x = element_text(size=8, face="bold"),
        axis.text.y = element_text(size=4, face="bold")
        ) +
    coord_flip()

dev.off()




#######################################
#if required:
#combine pheatmap with library(gridExtra)
a <- list(hmap1[[4]])
a[[2]] <- hmap2[[4]]
do.call(grid.arrange,a)
#plot(do.call(grid.arrange,a))

##########################################
#rbase like-plots generating ggplot codes
library(ggh4x)

ggplot(mpg, aes(model, hwy, fill=factor(cyl))) + 
    theme_classic() +
    geom_boxplot() + 
    guides(x = "axis_truncated", y= "axis_truncated") + #library(ggh4x)
    

    scale_y_continuous(name = "COOL", breaks=seq(10,50,10), limits=c(0,50), expand=c(0,0)) +
    labs(x="X-axis", y="COOL", fill="NEW LEGEND") + 

    theme(
        axis.ticks.length = unit(4, "pt"), #length of ticks

        axis.ticks = element_line(size = 3, color=c("red","blue","green")), #width of ticks
        #axis.ticks.x = element_line(size = 0.25, color="red"),
        #axis.ticks.y = element_line(size = 0.3),
        #axis.ticks.x = element_blank(), #rbase histogram style
        
        axis.line = element_line(color="black", size=0.25), #width of lines
        #axis.line.x = element_line(color="black", size=0.25),
        #axis.line.y = element_line(color="black", size=0.25),
        #axis.line.x = element_blank(),

        axis.text = element_text(angle=45, vjust = 0.5, hjust=1), #position of texts
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        #axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),

        axis.title = element_text(size=12, face="bold", color="black"),
        #axis.title.x=element_text(size=10),
        #axis.title.y=element_text(size=20, color="black"),

        legend.position="top", #bottom #right
        legend.title = element_text(color = "dodgerblue", size = 20),
        legend.text = element_text(color = "darkred"),
        legend.key.size = unit(1.2, "line")
        #legend.key = element_blank(),
        )


###########################################

# # # ## Read DNA methylation file
# input_file <- "~/Dropbox/DNAme_project/DNAme_ideas_output/Ideas_methpercent_table_states.txt"
# read_df <- fread(input_file, sep="\t", header=T) %>% as.data.frame()
# data <- read_df[,2:ncol(read_df)]
# mat_data <- as.matrix(data)

# # ## Rowwise scaling - interested in meth variance across cis-reg regions for any samples
# mat_data0 <- t(apply(mat_data, 1, scale))
# mat_data3 <- t(apply(mat_data, 1, function(x) (x-mean(x))/sd(x) ))


# ## Naming of cols and rows
# colnames(mat_data) <- colnames(data)
# rownames(mat_data) <- read_df$sample
# rownames(mat_data)
# colnames(mat_data)

# output_dir <- "~/Dropbox/DNAme_project/DNAme_ideas_output/plots"
# output_file_name <- file.path(output_dir, "Ideas_cisreg_methpercent_heatmap_figure_revised.pdf")       
# pdf(output_file_name)

# Heatmap(mat_data, name="DNAme Z-score", 
#     column_title="IDEAS Genomic States",
#     row_title="Cell-types",
#     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
#     column_title_gp = gpar(fontsize = 10, fontface = "bold"),
#     row_names_gp = gpar(fontsize = 8),
#     column_names_gp = gpar(fontsize = 6)
#     # km=4,
#     ## For heatmap and annotation legends (add new-legends)
#     # heatmap_legend_param = list(
#     # title_gp = gpar(fontsize = 10), 
#     # labels_gp = gpar(fontsize = 6),
#     # legend_direction = "horizontal",
#     # legend_width = unit(3, "cm"), title_position = "topcenter"
#     # )
# ) 

# dev.off()

######################################




# #Just adjust the widths:

# ggplot(df, aes(x=timepoint, y=mean, fill=group)) +
#   geom_bar(position=position_dodge(0.9), colour="black", stat="identity", width=0.9, , binwidth=0) +
#   geom_errorbar(position=position_dodge(0.9), width=0.85, aes(ymin=mean, ymax=mean+sem)) +
#   theme_bw() + coord_flip()


# ########################
# #citing example for plot
# #roadmap annotation
# roadmap_df <- fread(annotation_file, sep="\t")

# #first sync up with the matrix data
# anno_df <- left_join(read_df, roadmap_df, by=c("celltype"="EpigenomeID"))
# anno_df$log10PValue <- -log10(anno_df$PValue)


# df1 <- data.frame(Group_level=levels(factor(anno_df$Group)))
# df2 <- anno_df[,c("Group","Color")]
# df3 <- inner_join(df1, df2, by=c("Group_level"="Group"))
# color_df <- df3 %>% distinct(Group_level, .keep_all=T)
# color_df[color_df$Group_level == "ENCODE2012", "Color"] = "dodgerblue"

# #anno_df$Group <- factor(anno_df$Group, levels=anno_df$oldGroup), ordered=TRUE)

# #reverse the factor order if required
# #factor_reordered <- reorder(anno_df$celltype,anno_df$log10PValue)

# #library(ggrepel)
# options(ggrepel.max.overlaps = 10)
# label_data <- anno_df[rev(order(anno_df$log10PValue)),][c(1:5),]

# pdf(file.path(local_dir, "histone_H3K27Ac_marks_tl_topsnps.pdf"), width=9, height=7)
# anno_df %>%
#   ggplot(aes(x=reorder(celltype,log10PValue), y=log10PValue)) +
#   geom_point(size = 2, aes(color=Group)) + 
#   scale_color_manual(values=color_df$Color) +
#   #scale_y_continuous(limits = c(0, 10)) +
#   #scale_x_discrete(limits = rev(levels(factor_reordered)))+
#   geom_segment(aes(xend = celltype, yend = 0), size = 0.2, linetype=1) +
#   labs(y= "-log10(PValue)", x="CellType", title="H3K27ac Enrichment (n=98)") +
#   theme_minimal() + 
#   theme(axis.text.x = element_text(size=3, angle = 90, hjust = 1),
#     plot.title = element_text(hjust = 0.5)) +
#   geom_text_repel(data = label_data,
#                   aes(label = EpigenomeMnemonic), 
#                   box.padding = 1,
#                   #max.overlaps =10,
#                   show.legend = FALSE, #this removes the 'a' from the legend
#                   size=2
#                   )

# dev.off()


# #label_data <- bind_rows(label_data1, label_data2, data1, data2, data3, data4, data5, data6, data7, data8)
# #circle_point_data <- bind_rows(data1, data2, data3, data4, data5, 
# #  data6, data7, data8, data9, data12, data13, data14, data15, data16, data17)

# pdf(file.path(local_dir, "RemapNonRedundantTfbs_tl_topsnpsScatterplot.pdf"), width=10, height=8)

# read_df %>% 
#   ggplot(aes(x=enrich, y=log10PValue)) +
#   geom_point(color="red", size=1.3) + #geom_line() +
#   geom_vline(colour="dodgerblue", xintercept=1, linetype=2) +
#   geom_hline(colour="dodgerblue", yintercept=1, linetype=2) + 
#   #scale_color_manual(values=color_df$Color) +
#   #scale_x_discrete(limits = rev(levels(factor_reordered)))+
#   #geom_segment(aes(xend = celltype, yend = 0), size = 0.2, linetype=1) +
#   labs(y= "-log10P(Value)", x="Observed/Expected Fold-enrichment", title="Transcription Factor Binding Sites") +
#   theme_bw() + 
#   theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1),
#     plot.title = element_text(hjust = 0.5)) +
#   geom_text_repel(data = label_data,
#                 aes(label = TF), 
#                 box.padding = 1,
#                 #max.overlaps =10,
#                 show.legend = FALSE, #this removes the 'a' from the legend
#                 size=2.1
#                 ) +

#   geom_point(data=circle_point_data,
#              pch=21, fill=NA, size=2.5, colour="dodgerblue", stroke=1)

# dev.off()


# #########################################
# #example heatmap for reference
# input_file <- file.path(local_dir, infile)
# read_df <- fread(input_file, header=T) %>% 
#   as.data.frame() %>%
#   separate(Bed_File, c("celltype", "chromState", "txt"),  sep="[.]")


# #read_df %>% dcast(celltype+txt~chromState, value.var="PValue")
# read_tf <- read_df %>% dcast(celltype ~ chromState, value.var="PValue")
# rownames(read_tf) <- read_tf$celltype
# read_tf <- read_tf %>% select(-celltype) 

# mat_data <- as.matrix(read_tf)
# mat_data <- mat_data %>% head(111) #excluding encode cell-lines
# mat_data_log <- -log10(mat_data)
# #mat_data_log[is.na(mat_data_log)] <- 0

# ###########################################
# #rowname annotated heatmap
# roadmap_df <- fread(annotation_file, sep="\t")

# #first sync up with the matrix data
# anno_df <- left_join(data.frame(EpigenomeID=rownames(mat_data)), roadmap_df, by="EpigenomeID")

# #tabulate list of row annoations
# annotation_row <- data.frame(
#   CellType = anno_df$Group
#   #roadmap2 = anno_df$Group
#   )

# #sync tabulated list with matrix data rownames
# rownames(annotation_row) = rownames(mat_data) #or, anno_df$EpigenomeID

# #assign colors to thos tabulated lists
# mat_colors <- list(
#   CellType = anno_df$Color
#   #roadmap2 = anno_df$Color
#   )

# #sync colors with their initial group names
# names(mat_colors$CellType) <- anno_df$Group
# #names(mat_colors$roadmap2) <- anno_df$Group

# #mat_breaks <- quantile_breaks(mat_data_log, n = 100)

# pdf(file.path(local_dir, "chromHMM_15state_tl_topsnps_w-o.encode.pdf"), width = 7, height = 8)
# pdf(file.path(local_dir, "chromHMM_25state_tl_topsnps_w-o.encode.pdf"), width = 7, height = 8)

# pheat <- pheatmap(mat_data_log,
#   #fontsize=4, 
#   fontsize_row=3,
#   fontsize_col=9,
#   annotation_row = annotation_row,
#   #annotation_col = annotation_col,
#   annotation_colors = mat_colors,
#   annotation_legend = FALSE,
#   legend_labels = "-log10(PValue)",
#   #cellwidth=25, cellheight=3,
#   #cutree_rows = 2,
#   #cutree_cols = 5,
#   #breaks= mat_breaks,
#   main="CellType-Enrichment"
#   #color = colorRampPalette(brewer.pal(9,"Reds"))(400)
#   #colorRampPalette(rev(brewer.pal(n=5, name="RdBu")))(100)
# )

# dev.off()


# ## legends arrange in alphabetical order:
# anno_df <- anno_df %>% 
#     arrange(Group)

# #equivalent to drop_duplicates based on specific columns
# anno_df_uniq <- anno_df %>% distinct(Group, .keep_all = TRUE)

# pdf(file.path(local_dir, "roadmap_legend_point_nco1.pdf"), height=9)
# lgd = Legend(at = anno_df_uniq$Group, 
#     title = "CellType & Tissues", 
#     type = "points",
#     size=unit(4, "mm"),
#     #pch=rep(16, 49), 
#     background = "white",
#     legend_gp = gpar(col = anno_df_uniq$Color, cex = 14),
#     labels_gp = gpar(col = "black", fontsize=8), #fontface="bold")
#     ncol=1,
#     grid_height = unit(4, "mm")
#     #grid_width = unit(5, "mm")
#     #title_position = "topcenter")
#     )
# draw(lgd)
# dev.off()





#####################
#python plotting
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
