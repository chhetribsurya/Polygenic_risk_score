library(data.table)
library(argparser)
library(dplyr)

p<-arg_parser("Combine phenotypes and create predictor")
p<-add_argument(p, "--pheno", help="phenotype to load",  default="pheno.SAIGE.txt")
p<-add_argument(p, "--geno", help="genotype matrix to load",  default="mendrand.raw")
p<-add_argument(p, "--out", help="phenotype file to generate", default="pheno_mendrand_UKBB.txt")
p<-add_argument(p, "--prs", help="polygenic risk score file", default=tempfile())
argv<-parse_args(p)

data<-fread(argv$pheno)
geno<-fread(argv$geno)

geno<-geno[,c(2,eval(parse(text="7:ncol(geno)")))]
data<-merge(data, geno, by="IID")

#data1<-fread(paste("/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/", argv$prs, sep=""))
data1<-fread(argv$prs, header=F)
data1$pos<-paste(data1$V1, data1$V2, sep="_")
data1<-subset(data1, data1$pos %in% colnames(geno))
x<-paste(paste("(",data1$V3, ")", sep=""), paste("'", data1$pos, "'", sep=""), sep="*data$", collapse="+")

data$predictor<-eval(parse(text=x))

data$quartile<-ntile(data$predictor, 5)

fwrite(data, argv$out, quote=F, row.names=F, sep="\t", na="NA")

