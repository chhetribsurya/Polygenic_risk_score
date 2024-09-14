#!/usr/bin/env Rscript

suppressMessages({

    #load libraries:
    library(data.table)
    library(tidyverse)
    library(coloc)
    library(optparse, quietly=TRUE)

})

#timeit
start_time <- Sys.time()

#pass arguments:
#args <- commandArgs(TRUE)
#srcFile <- args[1]

option_list <- list(

    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
    help="Print extra output [default]"),

    make_option(c("-q", "--quietly"), action="store_false",
    dest="verbose", help="Print little output"),

    make_option(c("-i", "--infile"), type="character",
    help="Input file with snpid|betas|ses|mafs (w/ Full path) [default %default]",
    metavar="file"),

    make_option(c("-d", "--dirTarget"), default=getwd(),
    help = "Output file(s) target directory [default \"%default\"]",
    metavar="dir"),

    make_option(c("-o", "--outputName"), type="character", default="trait_coloc",
    help="Output file(s) base name for coloc output [default %default]",
    metavar="char"),

    make_option(c("-t", "--tissue_celltype"), default="Adipose_Subcutaneous",
    help="Provide trait name or cell-tissue type [default %default]",
    metavar="char"),

    make_option(c("-c", "--count"), type="integer", default=10,
    help="Min no. of eQTL signals per gene [default %default]",
    metavar="number"), #Max=10000

    make_option("--sd", type="integer", default=1,
    help="Standard deviation if generator == \"rnorm\" [default %default]",
    metavar="standard deviation")

)

#set defaults if options not provided:
#opt <- parse_args(OptionParser(option_list=option_list)) #usage = "%prog [options] file"
parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options


#print some progress messages to stderr if "quietly" wasn't requested
if(opt$verbose ) {
    write("writing some verbose output to standard error...\n", stderr())
}


############################
                           #
                           #
                           #
#main args:
colocfile <- opt$infile
output_dir <- opt$dirTarget
filebase <- opt$outputName

#additional args:
tissue <- opt$tissue_celltype
snpcount <- opt$count      
stdev <- opt$sd                
                           #
                           #
                           #
                           #
############################

print("Provided parameters...")
print(paste0("input file: ", colocfile))
print(paste0("target dir: ", output_dir))
print(paste0("filebase name: ", filebase))
print(paste0("tissue : ", tissue))
print(paste0("snp count: ", snpcount))
print(paste0("stdDev (pheno): ", stdev)); cat("\n\n")

#target output files:
#output_dir <- "/home-4/schhetr1@jhu.edu/surya/datasets/prs_project/coloc_prs_analysis/coloc_output"

if (!dir.exists(output_dir)){
dir.create(output_dir)
} else {
    print("Target dir exists...")
}


output1 <- file.path(output_dir, paste0(tissue, "_genecoloc_", filebase))
# output2 <- file.path(output_dir, paste0(tissue, "_snpcoloc_", filebase))
# output3 <- file.path(output_dir, paste0(tissue, "_finemap_ea_", filebase))
# output4 <- file.path(output_dir, paste0(tissue, "_finemap_aa_", filebase))


#read input file:
snpcount <- 10000
#colocfile <- "/Users/suryachhetri/datasets/prs_project/final_hg19/Adipose_subcutaneous_coloc_betas_ses_eqtlgwas_1e5.txt"
#colocfile <- "/work-zfs/abattle4/surya/datasets/prs_project/coloc_prs_analysis/coloc/Adipose_Subcutaneous_coloc_betas_ses_eqtlgwas_1e4.txt"
colocfile <- "/Users/suryachhetri/datasets/prs_project/final_hg19/Cells_EBV-transformed_lymphocytes_coloc_betas_ses_eqtlgwas_1e4.txt"
colocfile <- "/Users/suryachhetri/datasets/prs_project/final_hg19/coloc_geuvadis_montgo_eQTL/Cells_EBV-transformed_lymphocytes_coloc_betas_ses_eqtlgwas_1e4.txt"
colocfile <- "/Users/suryachhetri/datasets/prs_project/final_hg19/coloc_geuvadis_montgo_eQTL/Cells_EBV-transformed_lymphocytes_coloc_betas_ses_ONLYshared_eQTLnogwas_1e4.txt"
coloc_df <- fread(colocfile, sep="\t")
coloc_df$maf_AA <- 0.1

#filter for genes with > 1000 snps
# coloc_filt_df <- coloc_df %>%
#     count(phenotype_id, sort=T, name="count") %>%
#     filter(count < snpcount) %>%
#     as_tibble()

# print("Dimension of Coloc dataframe...")
# print(coloc_filt_df); cat("\n")

#create unique genelist:
#genevect <- coloc_df %>% distinct(phenotype_id) %>% unlist(use.name=F) %>% head
genevect <- coloc_df %>% distinct(phenotype_id) %>% unlist(use.name=F)

#distinct gene count:
coloc_df %>% 
    distinct(phenotype_id) %>% nrow

#create empty list for rbind:
gene_coloc <- list()
snp_coloc <- list()
snp_ids <- list()

fmap_ea_snps <- list()
fmap_aa_snps <- list()

idx_vect = c(1:length(genevect))
#masterlist <- list(genevect, idx_vect)

for(i in idx_vect) {
    cat("\n")
    print(paste("Processing gene number...:", i))
    cat("\n")
    #gene <- masterlist[[1]][i]
    gene <- genevect[i] 
    #gene <- "ENSG00000204314.10"
    #gene <- "ENSG00000246089.3"
    #print(paste0("processing gene: ", i))
    print(paste0("genename: ", gene));
    gene_df <- coloc_df %>% filter(phenotype_id == gene)
    gene_df <- gene_df %>% drop_na
    snpid <- list(snpid=gene_df$hg19_snp, snp=paste0("SNP.", seq(1:nrow(gene_df))))
    
    betas <- gene_df %>% 
        select(slope_EA, slope_AA) %>%
        rename(EA = slope_EA, AA = slope_AA)

    ses <- gene_df %>% 
        select(slope_se_EA, slope_se_AA) %>%
        rename(EA = slope_se_EA, AA = slope_se_AA)

    maf <- gene_df %>% 
        select(maf_EA, maf_AA) %>%
        rename(EA = maf_EA, AA = maf_AA)

    #colocalization:
    res <- coloc.abf(dataset1=list(beta=betas$EA, varbeta=ses$EA^2, N=nrow(betas),sdY=1,type="quant"),
                        dataset2=list(beta=betas$AA, varbeta=ses$AA^2, N=nrow(betas),sdY=1,type="quant"),
                        MAF=maf$EA, p1=1e-04, p2=1e-04, p12=5e-05) # p12=1e-05

    ##sensitivity(res, rule="H4 > 0.5")

    # res <- coloc.abf(dataset1=list(beta=betas$EA, varbeta=ses$EA^2, N=nrow(betas),sdY=1,type="quant"),
    #                     dataset2=list(beta=betas$AA, varbeta=ses$AA^2, N=nrow(betas),sdY=1,type="quant"),
    #                     MAF=maf$EA, p1=1e-03,p2=1e-03, p12=1e-04) # p12=1e-05


    # res <- coloc.abf(dataset1=list(beta=betas$EA, varbeta=ses$EA^2, N=nrow(betas),sdY=1,type="quant"),
    #                     dataset2=list(beta=betas$AA, varbeta=ses$AA^2, N=nrow(betas),sdY=1,type="quant"),
    #                     MAF=maf$EA, p1=1e-03,p2=1e-03, p12=1e-05) # p12=1e-05
    
    # res <- coloc.abf(dataset1=list(beta=betas$EA, varbeta=ses$EA^2, N=nrow(betas),sdY=1,type="quant"),
    #                     dataset2=list(beta=betas$AA, varbeta=ses$AA^2, N=nrow(betas),sdY=1,type="quant"),
    #                     MAF=maf$EA, p1=1e-05,p2=1e-05, p12=1e-06) # p12=1e-05


    #colocalization probablility:
    gene_coloc[[gene]] <- res$summary[c(2:6)]

    # snp_coloc[[gene]] <- res$results %>% arrange(desc(SNP.PP.H4))
    # snp_ids[[gene]] <- data.frame(snpid)

    # #fine mapping (w/ Null SNP) i.e probability of no causal snp
    # fmap_ea <- finemap.abf(dataset=list(beta=as.numeric(betas$EA), varbeta=as.numeric(ses$EA)^2, N=nrow(betas),sdY=1,type="quant")) #p1 = 1e-04
    # fmap_aa <- finemap.abf(dataset=list(beta=betas$AA, varbeta=ses$AA^2, N=nrow(betas),sdY=1,type="quant")) #p1 = 1e-04
    
    # #0.95 credible set
    # fmap_ea_res <- fmap_ea %>% 
    #     arrange(desc(SNP.PP)) %>% 
    #     filter(!is.na(V.)) %>% 
    #     mutate(cumsum=cumsum(SNP.PP)) %>% 
    #     filter(cumsum < 0.95) %>% 
    #     select(snp, SNP.PP, cumsum)

    # fmap_ea_snps[[gene]] <- fmap_ea_res

    # #0.95 credible set
    # fmap_aa_res <- fmap_aa %>% 
    #     arrange(desc(SNP.PP)) %>% 
    #     filter(!is.na(V.)) %>% 
    #     mutate(cumsum=cumsum(SNP.PP)) %>%
    #     filter(cumsum < 0.95) %>%
    #     select(snp, SNP.PP, cumsum)

    # fmap_aa_snps[[gene]] <- fmap_aa_res
    # cat("\n")

}

genecoloc_df <- bind_rows(!!!gene_coloc, .id="phenotype_id") %>% as_tibble()
#snpid_df <- bind_rows(snp_ids, .id="phenotype_id") %>%  as_tibble()

#snpcoloc_df <- bind_rows(snp_coloc,, .id="phenotype_id") %>%  select(phenotype_id, snp, SNP.PP.H4) %>% as_tibble()
#merged_snpcoloc <- merge(snpid_df, snpcoloc_df, by=c("phenotype_id", "snp")) %>% as_tibble()

# fmap_ea_df <- bind_rows(fmap_ea_snps, .id="phenotype_id") %>%  as_tibble()
# merged_fmap_ea <- merge(snpid_df, fmap_ea_df, by=c("phenotype_id", "snp")) %>% as_tibble()

# fmap_aa_df <- bind_rows(fmap_aa_snps, .id="phenotype_id") %>%  as_tibble()
# merged_fmap_aa <- merge(snpid_df, fmap_aa_df, by=c("phenotype_id", "snp")) %>% as_tibble()

# genecoloc_df %>% head
# snpcoloc_df %>% head
# fmap_ea_df %>% head
# fmap_aa_df %>% head

#output1 <- "/Users/suryachhetri/datasets/prs_project/coloc_output/Cells_EBV-transformed_lymphocytes_genecoloc_betas_ses_eqtlgwas_1e4.txt"
output1 <- "/Users/suryachhetri/datasets/prs_project/coloc_output/Cells_EBV-transformed_lymphocytes_genecoloc_betas_ses_geuvisMontgo_eqtlgwas_1e4.txt"
output1 <- "/Users/suryachhetri/datasets/prs_project/coloc_output/Cells_EBV-transformed_lymphocytes_genecoloc_betas_ses_geuvisMontgo_ONLYshared_eQTLnogwas_1e4.txt"
# # output2 <- "/Users/suryachhetri/datasets/prs_project/coloc_output/Cells_EBV-transformed_lymphocytes_snpcoloc_betas_ses_eqtlgwas_1e4.txt"
# # output3 <- "/Users/suryachhetri/datasets/prs_project/coloc_output/Cells_EBV-transformed_lymphocytes_finemap_ea_betas_ses_eqtlgwas_1e4.txt"
# output3 <- "/Users/suryachhetri/datasets/prs_project/coloc_output/Cells_EBV-transformed_lymphocytes_finemap_ea_betas_ses_geuvisMontgo_eqtlgwas_1e4.txt"
# # output4 <- "/Users/suryachhetri/datasets/prs_project/coloc_output/Cells_EBV-transformed_lymphocytes_finemap_aa_betas_ses_eqtlgwas_1e4.txt"
# output4 <- "/Users/suryachhetri/datasets/prs_project/coloc_output/Cells_EBV-transformed_lymphocytes_finemap_ea_betas_ses_geuvisMontgo_eqtlgwas_1e4.txt"

write.table(genecoloc_df, output1, quote=F, col.names=T)
# write.table(snpcoloc_df, output2, quote=F, col.names=T)
# write.table(fmap_ea_df, output3, quote=F, col.names=T)
# write.table(fmap_aa_df, output4, quote=F, col.names=T)

# # Saving on object in RData format
# save(data1, file = "data.RData")
# # Save multiple objects
# save(data1, data2, file = "data.RData")
# # To load the data again
# load("data.RData")

# output5 <- "/Users/suryachhetri/datasets/prs_project/coloc_output/Cells_EBV-transformed_lymphocytes_coloc_snpids_1e4.Rdata"
# save(snp_ids, file=output5)
# load(output5)

# #Save workspace image
# output5 <- "/Users/suryachhetri/datasets/prs_project/coloc_output/my_LCL_coloc_workspace.Rdata"
# save.image(file = output5)

# #To restore your workspace
# load("my_work_space.RData")

#write.table(genecoloc_df, output1, quote=F, col.names=T)
#write.table(merged_snpcoloc, output2, quote=F, col.names=T)
#write.table(merged_fmap_ea, output3, quote=F, col.names=T)
#write.table(merged_fmap_aa, output4, quote=F, col.names=T)

end_time <- Sys.time()
print(paste0("total run time: ", end_time - start_time))




#######################################################
# Download results from marcc:

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
####
#### Coloc profile with each coloc category:
####
##########################################################################


## Coloc heatmap
input_file <- "/Users/suryachhetri/datasets/prs_project/coloc_output/Adipose_subcutaneous_genecoloc_eqtl_1e5_final.txt"
input_file <- "/Users/suryachhetri/datasets/prs_project/coloc_output/Cells_EBV-transformed_lymphocytes_genecoloc_betas_ses_eqtlgwas_1e4.txt"
input_file <- "/Users/suryachhetri/datasets/prs_project/coloc_output/Cells_EBV-transformed_lymphocytes_genecoloc_betas_ses_geuvisMontgo_eqtlgwas_1e4.txt"
input_file <- "/Users/suryachhetri/datasets/prs_project/coloc_output/Cells_EBV-transformed_lymphocytes_genecoloc_betas_ses_geuvisMontgo_ONLYshared_eQTLnogwas_1e4.txt"
#input_file <- output1
read_df <- fread(input_file, header=F) %>% as.data.frame()
data <- read_df[,3:ncol(read_df)]
mat_data <- as.matrix(data)

## Rowwise scaling - interested in meth variance across cis-reg regions for any samples
#mat_data <- t(apply(mat_data, 1, scale))
mat_data[is.na(mat_data)] <- 0

## Naming of cols and rows
#colnames(mat_data) <- colnames(data)
#rownames(mat_data) <- read_df$sample
colnames(mat_data) <- c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")
rownames(mat_data) <- read_df[,2]
rownames(mat_data)
colnames(mat_data)

output_dir <- "/Users/suryachhetri/datasets/prs_project/coloc_output"
output_file_name <- file.path(output_dir, "Cells_EBV-transformed_lymphocytes_genecoloc_eqtl_1e4.pdf")       
output_file_name <- file.path(output_dir, "Cells_EBV-transformed_lymphocytes_genecoloc_geuvisMontgo_eqtl_1e4.pdf")       
output_file_name <- file.path(output_dir, "Cells_EBV-transformed_lymphocytes_genecoloc_geuvisMontgo_eqtl_1e4.pdf")       
pdf(output_file_name)
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

#########################
## Miscellaneous correlation plot:
input_file <- "/Users/suryachhetri/datasets/prs_project/coloc_output/ea_aa_cis_effects_correlation.txt"
read_df <- fread(input_file, header=T) %>% as.data.frame()
data <- read_df[,3:ncol(read_df)]
mat_data <- as.matrix(data)

## Rowwise scaling - interested in meth variance across cis-reg regions for any samples
#mat_data <- t(apply(mat_data, 1, scale))
mat_data[is.na(mat_data)] <- 0

## Naming of cols and rows
#colnames(mat_data) <- colnames(data)
#rownames(mat_data) <- read_df$sample
colnames(mat_data) <- colnames(data)
rownames(mat_data) <- read_df[,2]
rownames(mat_data)
colnames(mat_data)

output_dir <- "/Users/suryachhetri/datasets/prs_project/coloc_output"
output_file_name <- file.path(output_dir, "ea_aa_cis_effects_correlation.txt.pdf")       
pdf(output_file_name)
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



# Saving row names of cluster one
output_dir <- "~/Dropbox/DNAme_project/DMR_methprofiles/plots"
total_clusters <- length(row_order(ht7))
hmap_object <- ht7

extract_clusters <- function(heatmap_object, read_df, totalcluster){
    cluster <- list()
    for(i in 1:total_clusters){
        print(i)
        read_cluster <- read_df[row_order(ht7)[[i]],]
        cluster_df <- data.frame(do.call(rbind, strsplit(read_cluster[,c("DMR_idx")], "_")))
        names(cluster_df) <- c("DMR_idx", "Chrom", "start", "end")
        cluster_df <- cluster_df[,c("Chrom", "start", "end", "DMR_idx")]
        cluster_df$cluster <- i
        cluster[[i]] <- cluster_df
        write.table(cluster_df,file.path(output_dir, paste0("kmeans",total_clusters, "_cluster_" , i, "_df_rev.txt")), sep="\t", quote=F, row.names=F,col.names=F)
        #print(head(cluster_df))
    }
    return(cluster)
}

cluster_combined_outfile <- extract_clusters(hmap_object, read_df, total_clusters)
cluster_combined_df <- do.call(rbind, cluster_combined_outfile)
write.table(cluster_combined_df,file.path(output_dir, paste0("kmeans",total_clusters, "_cluster_all_df.txt")), sep="\t", quote=F, row.names=F,col.names=F)


# > genecoloc_df %>% data.frame %>% filter(PP.H0.abf > 0.5) %>% nrow
# [1] 1804
# > genecoloc_df %>% data.frame %>% filter(PP.H1.abf > 0.5) %>% nrow
# [1] 296
# > genecoloc_df %>% data.frame %>% filter(PP.H2.abf > 0.5) %>% nrow
# [1] 975
# > genecoloc_df %>% data.frame %>% filter(PP.H3.abf > 0.5) %>% nrow
# [1] 1198
# > genecoloc_df %>% data.frame %>% filter(PP.H4.abf > 0.5) %>% nrow
# [1] 4144

# > read_df %>% filter(V3>0.5) %>% nrow
# [1] 1577
# > read_df %>% filter(V4>0.5) %>% nrow
# [1] 98
# > read_df %>% filter(V5>0.5) %>% nrow
# [1] 1582
# > read_df %>% filter(V6>0.5) %>% nrow
# [1] 521
# > read_df %>% filter(V7>0.5) %>% nrow
# [1] 2723

#####
#without GWAS intersect:
# > genecoloc_df %>% data.frame %>% filter(PP.H0.abf > 0.5) %>% nrow
# [1] 417
# > genecoloc_df %>% data.frame %>% filter(PP.H1.abf > 0.5) %>% nrow
# [1] 204
# > genecoloc_df %>% data.frame %>% filter(PP.H2.abf > 0.5) %>% nrow
# [1] 750
# > genecoloc_df %>% data.frame %>% filter(PP.H3.abf > 0.5) %>% nrow
# [1] 1417
# > genecoloc_df %>% data.frame %>% filter(PP.H4.abf > 0.5) %>% nrow
# [1] 4093

