#!/usr/bin/bash

#usage: ./ld_value_eur.sh <snplistfile>
#example: ./ld_values_eur.sh PP.H.7/PP.H1.7.coloc.abf.topsnpList_forLDSC.txt

snplistfile=$1
snpfilename=$(basename $snplistfile)
plinkpath="/home-4/schhetr1@jhu.edu/surya/tools/plink1.9"
bfileset="/home-4/schhetr1@jhu.edu/surya/datasets/1KG/1000G_bial_europeans"
outfile="/home-4/schhetr1@jhu.edu/surya/datasets/1KG/coloc_topsnps/PP.H.75/${snpfilename}"

#for SNPlist:
#Obtaining LD values from a group of SNPs with other SNPs i.e with in 1MB of each SNP (--ld-window-kb 1MB: default parameter)
$plinkpath/plink --bfile $bfileset \
            --r2 \
            --ld-snp-list $snplistfile \
            --ld-window-kb 1000 \
            --ld-window-r2 0.8 \
            --threads 10 \
            --out $outfile

#sumtotal snpslist - oroginal snps + snp_a + snp_b
cat  $snplistfile <(tail -n+2 ${outfile}.ld| awk '{print $3}') <(tail -n+2 ${outfile}.ld| awk '{print $6}') | sort| uniq > ${outfile}.ldSnps.txt

#            --list-duplicate-vars suppress-first \
#for single snp
#Obtaining LD values for a specific SNP versus all others
#To obtain all LD values for a set of SNPs versus one specific SNP, use the --ld-snp command in conjunction with --r2. For example, to get a list of all values for every SNP within 1Mb of rs12345, use the command
#plink --file mydata \
#          --r2 \
#          --ld-snp rs12345 \
#          --ld-window-kb 1000 \ 
#          --ld-window-r2 0.8

#Obtaining a matrix of LD values
#Alternatively, it is possible to add the --matrix option, which creates a matrix of LD values rather than a list:
