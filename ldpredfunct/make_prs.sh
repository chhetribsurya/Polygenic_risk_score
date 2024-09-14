input=$1
output=$2

cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/JHS/merged_genotypes/server/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/CHS/merged_genotypes/CNV370v1/preimpute/race1/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/CHS/merged_genotypes/HumanOmni1_Quadv1/preimpute/imputed/race2/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/CHS/merged_genotypes/Hiseq2000/perlfiles/race1/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/MESA/merged_genotypes/race1/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/MESA/merged_genotypes/race2/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/MESA/merged_genotypes/race3/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/MESA/merged_genotypes/race4/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/ARIC/merged_genotypes/Hiseq2000/preimpute/race2/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/ARIC/merged_genotypes/Hiseq2000/preimpute/race1/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/ARIC/merged_genotypes/server/race1/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/ARIC/merged_genotypes/server/race2/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/Framingham/merged_genotypes/SHARE/server_imputed/affy/unrelated/
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/Framingham/merged_genotypes/SHARE/server_imputed/omni/unrelated/
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/GARNET/race5/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/SHARE/race3/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/PAGE/race3/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/PAGE/race4/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}
cd /work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/WHIMS/race5/unrelated/ 
/work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/plink2 --pfile unrelated_pfile --extract ${input} --ref-allele force ${input} --export A --out ${output} --rm-dup exclude-all
Rscript /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/change_pheno.r --pheno pheno_mendrand.txt --geno ${output}.raw --out ${output}.pheno.txt --prs ${input}

cd /work-zfs/abattle4/marios/BP_MR

