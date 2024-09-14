input=$1
output=$2
snp_field=$3
p_field=$4

/work-zfs/abattle4/marios/workfolder/marios/plink/plink --bfile /work-zfs/abattle4/marios/annotations/1kG_plink/1000G_hg19_chrpos_EUR --clump ${input} --clump-field ${p_field} --clump-p1 0.00000005 --clump-kb 1000 --clump-r2 0.001 --clump-snp-field ${snp_field} --exclude /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/duplicated_snps.snplist --out ${output}_clumped

awk '{if (NR>1) print $3}' ${output}_clumped.clumped | head -n -2 > ${output}_variants.clumped

head -n 1 ${input} > ${output}_pass1.txt
grep -w -F -f ${output}_variants.clumped ${input} | awk '{OFS="\t"; print}' >> ${output}_pass1.txt

rm ${output}_clumped* ${output}_variants.clumped

input=${output}_pass1.txt

/work-zfs/abattle4/marios/workfolder/marios/plink/plink --bfile /work-zfs/abattle4/marios/annotations/1kG_plink/1000G_hg19_chrpos_EUR --clump ${input} --clump-field ${p_field} --clump-p1 0.00000005 --clump-kb 1000000000 --clump-r2 0.1 --clump-snp-field ${snp_field} --exclude /work-zfs/abattle4/marios/BP_MR/UKBB_DBP65_GWAS/duplicated_snps.snplist --out ${output}_clumped

awk '{if (NR>1) print $3}' ${output}_clumped.clumped | head -n -2 > ${output}_variants.clumped

grep -w -F -f ${output}_variants.clumped ${input} | awk '{OFS="\t"; print $1, $2, $4}' > ${output}_prs.txt

rm ${output}_clumped* ${output}_variants.clumped ${output}_pass1.txt

