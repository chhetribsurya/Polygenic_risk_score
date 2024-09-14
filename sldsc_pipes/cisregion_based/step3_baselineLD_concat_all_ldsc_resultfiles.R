#!/usr/bin/env R

#####################
# R Concat dataframe
#####################

library(tidyverse)
library(data.table)
#library(stringr)

args <- commandArgs(trailingOnly=TRUE)

#argsvariable assign
prefix <- args[1]
print(paste("Prefix name:", prefix))

annotpath <- args[2]
print(paste("Annotpath name:", annotpath))

output_dir <- dirname(annotpath)
print(paste("Output_dir name:", output_dir))

#list files for processing
#results_dirlist <- list.files(file.path(annotpath, prefix), "^highly_heritableTrait_results*", include.dirs = TRUE, full.names=TRUE)
results_dirlist <- list.files(file.path(annotpath, prefix), "*results*", include.dirs = TRUE, full.names=TRUE)
print(paste("Results_dirlist:", results_dirlist))

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
print(head(result_concat_df))

#file names
#outfile1 <- file.path(output_dir, paste0("combined_highly_heritableTrait_results_", prefix, ".txt"))
outfile2 <- file.path(output_dir, paste0("combined_ALL_results_", prefix, ".txt"))

#write final df
#fwrite(result_concat_df, outfile1, sep="\t", row.names=F, quote=F, col.names=T)
fwrite(result_concat_df, outfile2, sep="\t", row.names=F, quote=F, col.names=T)

