#!/usr/bin/env Rscript

suppressMessages({

    #load libraries:
    library(data.table)
    library(tidyverse)
    library(coloc)
    library(optparse, quietly=TRUE)

})

#Usage: Rscript step2_run_coloc.R -i ~/datasets/prs_project/final_hg19/Adipose_Subcutaneous_coloc_betas_ses_eqtlgwas_1e4.txt -d ~/datasets/prs_project/final_hg19/coloc_test -o eqtl_1e4_final -t Adipose_Subcutaneous

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
#opt <- parse_args(OptionParser(option_list=option_list))
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
if (!dir.exists(output_dir)){
dir.create(output_dir)
} else {
    print("Target dir exists...")
}

output1 <- file.path(output_dir, paste0(tissue, "_genecoloc_", filebase))
output2 <- file.path(output_dir, paste0(tissue, "_snpcoloc_", filebase))
output3 <- file.path(output_dir, paste0(tissue, "_finemap_ea_", filebase))
output4 <- file.path(output_dir, paste0(tissue, "_finemap_aa_", filebase))


#read input file:
#snpcount <- 10000
#colocfile <- "/Users/suryachhetri/datasets/prs_project/final_hg19/Adipose_subcutaneous_coloc_betas_ses_eqtlgwas_1e5.txt"
coloc_df <- fread(colocfile, sep="\t")

#filter for genes with > 100 snps
coloc_filt_df <- coloc_df %>%
    count(phenotype_id, sort=T, name="count") %>%
    filter(count < snpcount) %>%
    as_tibble()

print("Dimension of Coloc dataframe...")
print(coloc_filt_df); cat("\n")

#create unique genelist:
#genevect <- coloc_df %>% distinct(phenotype_id) %>% unlist(use.name=F) %>% head(20)
genevect <- coloc_df %>% distinct(phenotype_id) %>% unlist(use.name=F)

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
    #gene <- "ENSG00000204314.10"
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
                        MAF=maf$EA) # p12=1e-05

    #colocalization probablility:
    gene_coloc[[gene]] <- res$summary[c(2:6)]
    snp_coloc[[gene]] <- res$results %>% arrange(desc(SNP.PP.H4))
    snp_ids[[gene]] <- data.frame(snpid)

    #fine mapping (w/ Null SNP) i.e probability of no causal snp
    fmap_ea <- finemap.abf(dataset=list(beta=as.numeric(betas$EA), varbeta=as.numeric(ses$EA)^2, N=nrow(betas),sdY=1,type="quant")) #p1 = 1e-04
    fmap_aa <- finemap.abf(dataset=list(beta=betas$AA, varbeta=ses$AA^2, N=nrow(betas),sdY=1,type="quant")) #p1 = 1e-04
    
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

print("test complete")
genecoloc_df <- bind_rows(!!!gene_coloc, .id="phenotype_id") %>%  as_tibble()
genecoloc_df %>% head
write.table(genecoloc_df, output1, quote=F, col.names=T)
snpid_df <- bind_rows(snp_ids, .id="phenotype_id") %>%  as_tibble()

snpcoloc_df <- bind_rows(snp_coloc,, .id="phenotype_id") %>%  select(phenotype_id, snp, SNP.PP.H4) %>% as_tibble()
merged_snpcoloc <- merge(snpid_df, snpcoloc_df, by=c("phenotype_id", "snp")) %>% as_tibble()

fmap_ea_df <- bind_rows(fmap_ea_snps, .id="phenotype_id") %>%  as_tibble()
merged_fmap_ea <- merge(snpid_df, fmap_ea_df, by=c("phenotype_id", "snp")) %>% as_tibble()

fmap_aa_df <- bind_rows(fmap_aa_snps, .id="phenotype_id") %>%  as_tibble()
merged_fmap_aa <- merge(snpid_df, fmap_aa_df, by=c("phenotype_id", "snp")) %>% as_tibble()

genecoloc_df %>% head
snpcoloc_df %>% head
fmap_ea_df %>% head
fmap_aa_df %>% head

#write.table(genecoloc_df, output1, quote=F, col.names=T)
write.table(merged_snpcoloc, output2, quote=F, col.names=T)
write.table(merged_fmap_ea, output3, quote=F, col.names=T)
write.table(merged_fmap_aa, output4, quote=F, col.names=T)

end_time <- Sys.time()
print(paste0("total runtime: ", end_time - start_time))

