#!/usr/bin/env R

#load libraries:
library(data.table)
library(tidyverse)
library(coloc)

#timeit
start_time <- Sys.time()

library("optparse")
 
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# input file:
file_path = "/Users/suryachhetri/datasets/prs_project/final_hg19"
tissue = "Adipose_Subcutaneous"
eqtlgwas_colocfile <- file.path(file_path, paste0(tissue, "_coloc_betas_ses_eqtlgwas_1e5.txt"))
eqtlgwas_colocfile <- file.path(file_path, paste0(tissue, "_coloc_betas_ses_eqtlgwas_1e4.txt"))

#read file:
coloc_df = fread(eqtlgwas_colocfile, sep="\t")

#create unique genelist:
genevect <- coloc_df %>% distinct(phenotype_id) %>% unlist(use.name=F) %>% head

#create empty list for rbind:
gene_coloc <- list()
snp_coloc <- list()
snp_ids <- list()

fmap_ea_snps <- list()
fmap_aa_snps <- list()

idx_vect = c(1:length(genevect))
#masterlist <- list(genevect, idx_vect)

for(i in idx_vect) {

    #gene <- masterlist[[1]][i]
    gene <- genevect[i] 
    print(paste0("processing gene: ", i))
    print(paste0("genename: ", gene));
    gene_df <- coloc_df %>% filter(phenotype_id == gene)
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
                        MAF=maf$EA)

    #colocalization probablility:
    gene_coloc[[gene]] <- res$summary[c(2:6)]
    snp_coloc[[gene]] <- res$results %>% arrange(desc(SNP.PP.H4))
    snp_ids[[gene]] <- data.frame(snpid)

    #fine mapping (w/ Null SNP) i.e probability of no causal snp
    fmap_ea <- finemap.abf(dataset=list(beta=betas$EA, varbeta=ses$EA^2, N=nrow(betas),sdY=1,type="quant"))
    fmap_aa <- finemap.abf(dataset=list(beta=betas$AA, varbeta=ses$AA^2, N=nrow(betas),sdY=1,type="quant"))
    
    #0.95 credible set
    fmap_ea_res <- fmap_ea %>% 
        arrange(desc(SNP.PP)) %>% 
        filter(!is.na(V.)) %>% 
        mutate(cumsum=cumsum(SNP.PP)) %>% 
        filter(cumsum < 0.95) %>% 
        select(snp, SNP.PP, cumsum)

    fmap_ea_snps[[gene]] <- fmap_ea_res

    #0.95 credible set
    fmap_aa_res <- fmap_aa %>% 
        arrange(desc(SNP.PP)) %>% 
        filter(!is.na(V.)) %>% 
        mutate(cumsum=cumsum(SNP.PP)) %>%
        filter(cumsum < 0.95) %>%
        select(snp, SNP.PP, cumsum)

    fmap_aa_snps[[gene]] <- fmap_aa_res
    cat("\n")

 }


genecoloc_df <- bind_rows(gene_coloc, .id="phenotype_id") %>%  as_tibble()
snpid_df <- bind_rows(snp_ids, .id="phenotype_id") %>%  as_tibble()

snpcoloc_df <- bind_rows(snp_coloc,, .id="phenotype_id") %>%  select(phenotype_id, snp, SNP.PP.H4) %>% as_tibble()
merged_snpcoloc <- merge(snpid_df, snpcoloc_df, by=c("phenotype_id", "snp")) %>% as_tibble()

fmap_ea_df <- bind_rows(fmap_ea_snps, .id="phenotype_id") %>%  as_tibble()
merged_fmap_ea <- merge(snpid_df, fmap_ea_df, by=c("phenotype_id", "snp")) %>% as_tibble()

fmap_aa_df <- bind_rows(fmap_aa_snps, .id="phenotype_id") %>%  as_tibble()
merged_fmap_aa <- merge(snpid_df, fmap_aa_df, by=c("phenotype_id", "snp")) %>% as_tibble()

write.table(genecoloc_df, file.path(file_path, paste0(tissue, "_genecoloc_1e4.txt")), quote=F, col.names=T)
write.table(merged_snpcoloc, file.path(file_path, paste0(tissue, "_snpcoloc_1e4.txt")), quote=F, col.names=T)
write.table(merged_fmap_ea, file.path(file_path, paste0(tissue, "_finemap_ea_1e4.txt")), quote=F, col.names=T)
write.table(merged_fmap_aa, file.path(file_path, paste0(tissue, "_finemap_ea_1e4.txt")), quote=F, col.names=T)

end_time <- Sys.time()
end_time - start_time

