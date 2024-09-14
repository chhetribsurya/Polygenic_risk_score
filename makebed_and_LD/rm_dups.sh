#sometimes reference vcf contains duplicate variant IDs
#remove all those including Varinat "." IDs

#remove duplicated IDs from European sample
cut -f 2 1000G_bial_nochild_europeans.bim | sort | uniq -d > /home-4/schhetr1@jhu.edu/surya/datasets/1KG/dups_eur.txt

#remove duplicated IDs from African sample
cut -f 2 1000G_bial_nochild_africans.bim | sort | uniq -d > /home-4/schhetr1@jhu.edu/surya/datasets/1KG/dups_afr.txt

#NOTE: dups_afr.txt and dups_eur.txt are infact same set of variants are exact copy as the source is same. Thus one run is suffice to remove duplicates 

#Regenerate bed file with filtered duplicates:
#make bedfile for africans with filtered out dups
plink --bfile dipt_hg38/1KG/1000G_bial_nochild_africans --exclude dups_afr.txt --make-bed --out ./1000G_bial_africans
 
#make bedfile for europeans with filtered out dups
plink --bfile dipt_hg38/1KG/1000G_bial_nochild_europeans --exclude dups_afr.txt --make-bed --out ./1000G_bial_europeans

