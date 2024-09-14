prefix="EUR"
prefix="AA"

plinkfile="/work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/WHIMS/race5/unrelated"
plinkfile="/work-zfs/abattle4/marios/HF_GWAS/cohorts/WHI/merged_genotypes/server_download/PAGE/race3/unrelated"

plinkpath="/work-zfs/abattle4/surya/tools/plink2.9"

for chr in {1..22}; do \
${plinkpath}/plink2 --pgen $plinkfile/unrelated_pfile.pgen --pvar $plinkfile/unrelated_pfile.pvar --psam $plinkfile/unrelated_pfile.psam --chr $chr --make-bed --out ./whimsUnrelated_${prefix}_${chr}; \
done
