#!/usr/bin/env R

import pandas as pd

df_meta = pd.read_csv("full_variant_qc_metrics.txt", sep="\t")
df_data = pd.read_csv("continuous-SBP-both_sexes-auto_medadj_irnt.tsv", sep="\t")

select_metacol = ['chrom', 'pos', 'ref', 'alt', 'rsid', 'varid', 'nearest_genes']
df_meta_filt = df_meta.loc[:,select_metacol]
df_meta_filt.rename(columns={"chrom": "chr"}, inplace=True)

select_datacol = ["chr", "pos","alt", "af_AFR", "beta_AFR", "se_AFR", "pval_AFR", "af_EUR","beta_EUR","se_AFR", "pval_EUR"]
df_data_filt = df_data.loc[:,select_datacol]
df_data_filtdrop = df_data_filt.dropna(subset=["pval_AFR", "pval_EUR"], how="any")
df_data_filtdrop.shape

#df_data_filt.shape
#df_data_filt.head()
#df_meta_filt.shape
#df_meta_filt.head()
#df_data_head = df_data_filt.head(50)
#df_data_head.dropna(subset=["pval_AFR", "pval_EUR"], how="any")

df_merged = pd.merge(df_data_filtdrop, df_meta_filt, on=["chr", "pos","alt"], how="inner")
df_merged.to_csv("continuous_SBP_gwas_merged_final.txt", sep="\t", index=False)


#######

library(tidyverse)
library(data.table)
library("CMplot")

df <- fread("continuous_SBP_gwas_merged_final.txt", sep="\t")
df_sbp_pval <- df %>% data.frame %>% dplyr::select(nearest_genes,chr,pos,pval_AFR, pval_EUR)
df_sbp_pvalAFR <- df_sbp_pval %>% dplyr::select(nearest_genes,chr,pos,pval_AFR)
df_sbp_pvalEUR <- df_sbp_pval %>% dplyr::select(nearest_genes,chr,pos,pval_EUR)

df_sbp_beta <- df %>% data.frame %>% dplyr::select(nearest_genes,chr,pos,beta_AFR, beta_EUR)
df_sbp_betaAFR <- df_sbp_beta %>% dplyr::select(nearest_genes,chr,pos,beta_AFR)
df_sbp_betaEUR <- df_sbp_beta %>% dplyr::select(nearest_genes,chr,pos,beta_EUR)

#df_sbp_pval <- df_sbp_pval %>% replace_na(list(pval_AFR=0, pval_EUR=0, nearest_genes="unknown"))

# df_test <- df_sbp_pval %>% filter(pval_EUR < 1e-6)
# CMplot(df_test, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=c(1e-8,1e-10),threshold.lty=c(1,2),
#         threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
#         chr.den.col=c("darkgreen", "yellow", "red"),signal.col=c("red","green"),cex=c(0.4), signal.cex=c(0.5,0.5),
#         signal.pch=c(19,19),file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
#         width=14,height=6)

SNPs <- list(
	df_sbp_pval$nearest_genes[df_sbp_pval$pval_AFR<1e-6]
	#df_sbp_pval$nearest_genes[df_sbp_pval$pval_EUR<1e-70]
)

#Single plot with AFR thresholds and SNP highlights
CMplot(df_sbp_pvalAFR, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=c(1e-6,1e-4),threshold.lty=c(1,2),
        threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
        chr.den.col=c("darkgreen", "yellow", "red"),cex=c(0.3), signal.col=c("red","green"),signal.cex=c(0.5,0.5),
        signal.pch=c(19,19),file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
        highlight=SNPs, highlight.text=SNPs, highlight.text.cex=0.7, main="AFR SBP GWAS", width=14,height=6)

#FinalSingle plot with AFR thresholds and without SNP highlights
CMplot(df_sbp_pvalAFR, plot.type="m", LOG10=TRUE, ylim=NULL, threshold=c(1e-6,1e-4),threshold.lty=c(1,2),
        threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
        chr.den.col=c("darkgreen", "yellow", "red"),cex=c(0.4), signal.col=c("red","green"),signal.cex=c(1,1),
        signal.pch=c(19,19),file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
        main="AFR SBP GWAS", width=14,height=6)

#Single plot without (EUR threshold, SNP highlights, Y limits)
CMplot(df_sbp_pvalEUR, plot.type="m", LOG10=TRUE, ylim=NULL, 
         threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
        chr.den.col=c("darkgreen", "yellow", "red"),cex=c(0.2), signal.col=c("red","green"),signal.cex=c(0.3,0.3),
        signal.pch=c(19,19), threshold.lwd=c(1,1), threshold=c(1e-6,1e-4),threshold.lty=c(1,2),, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
        main="EUR SBP GWAS", width=14,height=6)

#FinalSingle plot without EUR thresholds and without SNP highlights and with ylimits
CMplot(df_sbp_pvalEUR, plot.type="m", LOG10=TRUE, ylim=c(0,20), 
         threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
        chr.den.col=c("darkgreen", "yellow", "red"),cex=c(0.3),threshold.lwd=c(1.5,1.5), threshold=c(1e-6,1e-4),threshold.lty=c(1,2),, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
        main="EUR SBP GWAS", width=14,height=6)

#Multiplot track with AFR EUR SNP highlights only
SNPs <- list(
	df_sbp_pval$nearest_genes[df_sbp_pval$pval_AFR<1e-6],
	df_sbp_pval$nearest_genes[df_sbp_pval$pval_EUR<1e-70]
)

CMplot(df_sbp_pval, plot.type="m",multracks=TRUE,threshold=c(1e-6,1e-4),threshold.lty=c(1,2), 
        threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
        chr.den.col=c("darkgreen", "yellow", "red"), cex=c(0.3), signal.col=c("red","green","blue"),
        signal.cex=0.4, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
        highlight=SNPs, highlight.text=SNPs, highlight.text.cex=0.7)

#FinalMultiplot track without SNP highlights
CMplot(df_sbp_pval, plot.type="m",multracks=TRUE,threshold=c(1e-6,1e-4),threshold.lty=c(1,2), 
        threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
        chr.den.col=c("darkgreen", "yellow", "red"), cex=c(0.3), signal.col=c("red","green","blue"),
        signal.cex=0.4, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)

#Note: if you are not supposed to change the color of signal, 
#          please set signal.col=NULL and highlight.col=NULL.

#FinalDensity plot
CMplot(df_sbp_pval,type="p",plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,
    main="SBP UKBB GWAS",file.output=TRUE,verbose=TRUE,width=9,height=6)

CMplot(df_sbp_betaAFR, type="h",plot.type="m", band=0.5, LOG10=FALSE, ylab="SNP effect",ylim=c(-0.02,0.02),
        threshold.lty=2, threshold.lwd=1, threshold.col="red", amplify=FALSE,cex=0.5,
        chr.den.col=NULL, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


CMplot(df_sbp_betaEUR, type="h",plot.type="m", band=0.5, LOG10=FALSE, ylab="SNP effect",ylim=c(-0.02,0.02),
        threshold.lty=2, threshold.lwd=1, threshold.col="red", amplify=FALSE,cex=0.5,
        chr.den.col=NULL, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE)


# users can personally set the windowsize and the min/max of legend by:
# bin.size=1e6
# bin.range=c(min, max)
# memo: add a character to the output file name
# chr.labels: change the chromosome names
# main: change the title of the plots, for manhattan plot, if there are more than one trait, main can be
#       assigned as a character vector containing the desired title for each trait