#####################################################################################
# test codes:
gene_list %>% slice(1:3) %>% as.list()
test_func <- function(x){paste0(x, "test")}
gene_list %>% slice(1:3) %>% map(., function(x){paste0(x, "test")})
gene_list %>% slice(1:3) %>% map(., test_func)

# dataframe is a list of vectors:
list(A = c("a1","a2"), B = c("b1", "b2")) %>% data.frame()
$A
[1] "a1" "a2"
$B
[1] "b1" "b2"

#converts list to vectors
unlist(list(A = c("a1","a2", "2"), B = c("b1", "b2", "1")))
  A1   A2   A3   B1   B2   B3
"a1" "a2"  "2" "b1" "b2"  "1"

list("a", "b", "c") %>% data.frame
#converts list to vectors
unlist(list("a", "b", "c")) 
[1] "a" "b" "c"

top3 <- gene_list %>% slice(1:3) %>% as.list()

#gene_list %>% as.list() %>% map_chr(1)
#gene_list %>% map_chr(1)

# subsetting list and vectors
gene_list %>% head() %>% as.list() %>% unlist(use.names = FALSE) %>% .[c(1:2)]
gene_list %>% head() %>% as.list() %>% .[[1]] %>% .[c(1:2)]

# subsetting the list of genes of interest:
gene_list[gene_list$phenotype_id %in% (top3 %>% as.list() %>% unlist())]
gene_list[which(gene_list$phenotype_id == (top3 %>% as.list() %>% unlist()))]
target <- top3 %>% as.list() %>% unlist(use.name=F)
gene_list %>% unlist(use.name=F) %>% .[c(1:3)]


filter(gene_list, gene_list$phenotype_id %in% target)

gene_list %>% 
    filter(gene_list$phenotype_id %in% target)

####################################################################################
gene <- gene_list %>% slice(1) %>% unlist(use.name=F)

filter(coloc_df, phenotype_id == gene)

pheno_count <- coloc_df %>%
    group_by(phenotype_id) %>%
    summarise(count=n())

pergene_count <- coloc_df %>%
    count(phenotype_id, sort=T) %>%
    mutate(quant = quantile(n))

quantile(pergene_count$n, probs=c(seq(0,1, 0.05)))
quantile(pergene_count$n, 0.05) #5th percentile

toppergene_count <- coloc_df %>%
    count(phenotype_id, sort=T, name="count") %>%
    filter(!(count > quantile(count, 0.005)))
######################################################################

pergene_count <- coloc_df %>%
    count(phenotype_id, sort=T, name="count") %>%
    filter(count > 100)

for(gene in gene_list){

gene="ENSG00000204314.10"
gene_coloc <- coloc_df %>%
    filter(phenotype_id == gene)

traits <- c("EA", "AA", "GWAS")
snpid <-  gene_coloc$hg19_snp

betas <- gene_coloc %>% 
    select(slope_EA, slope_AA, beta) %>% 
    as.matrix()

#set index with column:
rownames(betas) <- snpid
colnames(betas) <- traits

ses <- gene_coloc %>% 
    select(slope_se_EA, slope_se_AA, se) %>% 
    as.matrix()

#set index with column:
rownames(ses) <- snpid
colnames(ses) <- traits

res <- hyprcoloc(betas, ses, trait.names=traits, snp.id=snpid, snpscores=TRUE);
res <- hyprcoloc(betas, ses, trait.names=traits, snp.id=snpid, snpscores=TRUE, uniform.priors = TRUE, trait.subset = c("EA","AA"));

}

write.table(betas, "betas.txt", sep="\t", quote=FALSE)
write.table(ses, "ses.txt", sep="\t", quote=FALSE)
write.table(snpid, "snpid.txt", sep="\t", quote=FALSE)

#colocalization:
my.res <- coloc.abf(dataset1=list(beta=beta_df$EA, varbeta=ses_df$EA^2, N=length(beta_df$EA),sdY=1,type="quant"),
                    dataset2=list(beta=beta_df$AA, varbeta=ses_df$AA^2, N=length(beta_df$EA),sdY=1,type="quant"),
                    MAF=runif(length(beta_df$EA)), 0, 0.5)

#colocalization probablility:
my.res$summary
my.res$summary[6] 

#snps coloclization probability:
new_res <- my.res$results
new_res %>% arrange(desc(SNP.PP.H4)) %>% head

#####################
> my.res$summary
       nsnps    PP.H0.abf    PP.H1.abf    PP.H2.abf    PP.H3.abf    PP.H4.abf
8.860000e+02 2.156012e-03 0.000000e+00 9.978251e-01 0.000000e+00 1.887330e-05

> new_res %>% arrange(desc(SNP.PP.H4)) %>% head
      snp       V.df1    z.df1     r.df1 lABF.df1      V.df2       z.df2
1 SNP.857 0.006191828 3.707016 0.7841954 4.621502 0.04017783 -0.07739199
2 SNP.864 0.013397042 2.685644 0.6267926 1.767618 0.07011801  1.58824968
3 SNP.866 0.013397042 2.685644 0.6267926 1.767618 0.07011801  1.58824968

      r.df2    lABF.df2 internal.sum.lABF   SNP.PP.H4
1 0.3589786 -0.22127115          4.400231 0.093067820
2 0.2429333  0.16725217          1.934870 0.007908734
3 0.2429333  0.16725217          1.934870 0.007908734
###############################################

#fine mapping (last snp Null SNP - probability that this list doesn't have a causal snp)
my.res <- finemap.abf(dataset=list(beta=beta_df$EA, varbeta=ses_df$EA^2, N=length(beta_df$EA),sdY=1,type="quant"))

#order the fine mapped loci: (filter(is.na))
head(my.res[order(my.res$SNP.PP,decreasing = T),])
my.res %>% arrange(desc(SNP.PP)) %>% head %>% filter(!is.na(V.))

####################
> my.res %>% arrange(desc(SNP.PP)) %>% head %>% filter(is.na(V.))
  V. z. r. lABF.  snp  prior    SNP.PP
1 NA NA NA     0 null 0.9114 0.9128988
####################


#Miscellanous with hyprcoloc:
#cbind :column bind; rbind for row bind:
betas.1 = cbind(betas[,1],betas[,1]); ses.1 = cbind(ses[,1],ses[,1])
cor((betas[,1]/ses[,1]), (betas[,2]/ses[,2]))

res <- hyprcoloc(betas, ses, trait.names=c("EA", "AA"), snp.id=snpid, snpscores=TRUE, bb.alg=F);
res <- hyprcoloc(betas.1, ses.1, trait.names=c("EA", "AA"), snp.id=snpid, snpscores=TRUE, uniform.priors = TRUE, bb.alg=F);



prior.options = c(1e-4, 1e-10, 1e-20, 1e-25, 1e-100);
for(i in prior.options){
  res <- hyprcoloc(betas.1, ses.1, trait.names=traits, snp.id=snpid, snpscores=TRUE, bb.alg=F,
                   uniform.priors = TRUE, prior.1 = i, reg.steps = 2);
  print(paste0("prior.1 = ",i));
  print(res);
  }


###############################################
# using dplyr functions:
library(data.table)
library(tidyverse)
library(hyprcoloc)
library(coloc)

# input file:
output_path = "/Users/suryachhetri/datasets/prs_project/final_hg19"
tissue = "Adipose_Subcutaneous"
eqtlgwas_colocfile <- file.path(output_path, paste0(tissue, "_coloc_betas_ses_eqtlgwas.txt"))

# read file:
coloc_df = fread(eqtlgwas_colocfile, sep="\t")
# coloc_df %>% head
# coloc_df %>%
#     distinct(phenotype_id) %>% 
#     dim()

# #create unique genelist:
# gene_list <- coloc_df %>%
#     distinct(phenotype_id)

# #colocalization:
# my.res <- coloc.abf(dataset1=list(beta=beta_df$EA, varbeta=ses_df$EA^2, N=length(beta_df$EA),sdY=1,type="quant"),
#                     dataset2=list(beta=beta_df$AA, varbeta=ses_df$AA^2, N=length(beta_df$EA),sdY=1,type="quant"),
#                     MAF=runif(length(beta_df$EA)), 0, 0.5)

# #colocalization probablility:
# my.res$summary
# my.res$summary[6] 

# #snps coloclization probability:
# new_res <- my.res$results
# new_res %>% arrange(desc(SNP.PP.H4)) %>% head


# #fine mapping (last snp Null SNP - probability that this list doesn't have a causal snp)
# my.res <- finemap.abf(dataset=list(beta=beta_df$EA, varbeta=ses_df$EA^2, N=length(beta_df$EA),sdY=1,type="quant"))
# my.res %>% arrange(desc(SNP.PP)) %>% head %>% filter(!is.na(V.))


# #iteration with genes
# gene_subset <- gene_list %>% slice(1:5) %>% unlist(use.name=F)
# cols <- colnames(coloc_data) %>% unlist()
# gene="ENSG00000227232.5"
# gene="ENSG00000268903.1"
# gene="ENSG00000269981.1"

#Colocalization:
#create unique genelist:
start_time <- Sys.time()

genelist <- coloc_df %>% distinct(phenotype_id) %>% unlist(use.name=F)
gene_coloc <- list()
snp_coloc <- list()
snp_ids <- list()

index <- c(1:length(genelist))
for(gene in genelist){

    gene_df <- coloc_df %>% filter(phenotype_id == gene)
    snpid <- list(snpid=gene_df$hg19_snp, snp=paste0("SNP.", seq(1:nrow(gene_df))))
    

    betas <- gene_df %>% 
        select(slope_EA, slope_AA) %>%
        rename(EA = slope_EA, AA = slope_AA)

    ses <- gene_df %>% 
        select(slope_se_EA, slope_se_AA) %>%
        rename(EA = slope_se_EA, AA = slope_se_AA)

    #colocalization:
    res <- coloc.abf(dataset1=list(beta=betas$EA, varbeta=ses$EA^2, N=nrow(betas),sdY=1,type="quant"),
                        dataset2=list(beta=betas$AA, varbeta=ses$AA^2, N=nrow(betas),sdY=1,type="quant"),
                        MAF=runif(length(betas$EA)), 0, 0.5)

    #colocalization probablility:
    gene_coloc[[gene]] <- res$summary[c(2:6)]
    res_snp <- res$results
    snp_coloc[[gene]] <- res_snp %>% arrange(desc(SNP.PP.H4)) %>% head
    snp_ids[[gene]] <- data.frame(snpid)

    # coloc_output[[gene]][["trait_prob"]] <- res$summary[c(2:6)]
    # res_snp <- res$results
    # coloc_output[[gene]][["snp_prob"]] <- res_snp %>% arrange(desc(SNP.PP.H4)) %>% head

    #snps coloclization probability:
    res_snp <- res$results
    res_snpcoloc <- res_snp %>% arrange(desc(SNP.PP.H4)) %>% head

    #fine mapping (last snp Null SNP - probability that this list doesn't have a causal snp)
    fmap_ea <- finemap.abf(dataset=list(beta=betas$EA, varbeta=ses$EA^2, N=nrow(betas),sdY=1,type="quant"))
    fmap_aa <- finemap.abf(dataset=list(beta=betas$AA, varbeta=ses$AA^2, N=nrow(betas),sdY=1,type="quant"))
    fmap_ea_res <- fmap_ea %>% arrange(desc(SNP.PP)) %>% head %>% filter(!is.na(V.))
    fmap_aa_res <- fmap_aa %>% arrange(desc(SNP.PP)) %>% head %>% filter(!is.na(V.))
    cat("\n")
    #cat(paste0("\n\n", gene, "; coloc prob = ", res$summary[6]))

  }

genecoloc_df <- do.call(rbind, gene_coloc)
end_time <- Sys.time()
end_time - start_time

snpid_df <- do.call(rbind, snp_ids)
snpcoloc_df <- do.call(rbind, snp_coloc) %>% select(snp, SNP.PP.H4)
end_time <- Sys.time()
end_time - start_time

#zip/enumerate equivalent 
#method 1 
d <- list(num=c(1,2,3), char=c("foo", "bar", "baz"))
for(i in c(1:length(d[[1]]))) {
    elem1 = d[[1]][i]
    elem2 = d[[2]][i]
  print(c(elem1, elem2))
}


#method 2
genelist=c("foo", "bar", "baz")
index=c(1:length(genelist))
d = list(index=index, genelist=genelist)
for(i in index) {
    elem1 = d[[1]][i]
    elem2 = d[[2]][i]
  print(c(elem1, elem2))
}


#method 3
x = seq(1,10)
j = seq(11,20)

for (i in 1:length(x)){

    print (c(x[i],j[i]))
}


index <- c(1, 2, 3)
chars <- c('foo', 'bar', 'baz')
vect <- c('a', 'b', 'c')

for(i in index){
    print(c(index[i], chars[i], vect[i]))
}

index <- c(1, 2, 3)
chars <- c('foo', 'bar', 'baz')
vect <- c('a', 'b', 'c')

for(i in index){
    print(c(index[i], chars[i], vect[i]))
}

vect1 <- c(1, 2, 3)
vect1 <- c('foo', 'bar', 'baz')
vect2 <- c('a', 'b', 'c')


idx_list <- list(vect1, vect2)
idx_vect <- c(1:length(idx_list[[1]]))

for(i in idx_vect){
    x <- idx_list[[1]][i]
    j <- idx_list[[2]][i]
    print(c(i, x, j))
}


    vect1 <- c('A', 'B', 'C')
    vect2 <- c('a', 'b', 'c')
    
    # eqiv to zip values:
    idx_list <- list(vect1, vect2)
    idx_vect <- c(1:length(idx_list[[1]]))
    
    for(i in idx_vect){
        x <- idx_list[[1]][i]
        j <- idx_list[[2]][i]
        print(c(i, x, j))
    }

