
get_interaction<-function(data, snpset){
 d<-data
 colnames(d)[colnames(d) %in% snpset]<-c("one", "two")
 if(length(which(!is.na(d$one)))>0&length(which(!is.na(d$two)))>0){
 if(sd(d$one, na.rm=T)>0&sd(d$two, na.rm=T)>0){
 newd<-d[!is.na(one)&(!is.na(two)),]
 if(length(unique(newd$gender))>1&length(unique(newd$study))>1){
 mod<-lm(DBP~one*two+factor(gender)+factor(study)+as.matrix(d[,..x]), data=d)
 }else if(length(unique(newd$gender))==1&length(unique(newd$study))>1){
 mod<-lm(DBP~one*two+factor(study)+as.matrix(d[,..x]), data=d)
 }else if(length(unique(newd$gender))>1&length(unique(newd$study))==1){
 mod<-lm(DBP~one*two+factor(gender)+as.matrix(d[,..x]), data=d)
 }else{
 mod<-lm(DBP~one*two+as.matrix(d[,..x]), data=d)
 }
 res<-summary(mod)$coefficients
 res_snp1_p<-res[2,4]
 res_snp1_b<-res[2,1]
 res_snp2_p<-res[3,4]
 res_snp2_b<-res[3,1]
 res_inter_b<-res[dim(res)[1],1]
 res_inter_p<-res[dim(res)[1],4]
 return(c(snpset, res_snp1_b, res_snp1_p, res_snp2_b, res_snp2_p, res_inter_b, res_inter_p))
}}}

library(data.table)

load("Revision_multivariable_forest_data.rdata")
data<-multivar_data
x<-colnames(data)[grepl("PC", colnames(data))]
snps<-colnames(data)[grepl(":", colnames(data))]
k<-data.frame(a1=snps, a2="age")
a<-as.matrix(t(k))

res<-apply(a, 2, function(x) get_interaction(data, x))

save(data, res, x, file="DBP_META_prs_age_interaction_testing.rdata")

print("Done beauty!!")
