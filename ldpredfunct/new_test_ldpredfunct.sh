# Input files:
outputdir=debug

mkdir -p $outputdir

plinkfile="test/TSI_[1:22]"
phenotype="test/TSI_simulated_trait.txt"
phenotype="debug/TSI_simulated_trait.txt"
statsfile="test/summary_statistics_traininig.txt"
functfile=test/test_functfile.txt

outCoord="$outputdir/Coord_Final"; 
outLdpredfunct="$outputdir/ldpredfunct_posterior_means"; 
outValidate="$outputdir/ldpredfunct_prs"; 

# Output files
#outCoord="test/Coord_Final"
#outLdpredfunct="test/ldpredfunct_posterior_means"
#outValidate="test/ldpredfunct_prs"

N=100
h2=0.0000000001
K=5

#N=373000;
#h2=0.1512;
#K=100

#source activate ldsc-new
eval "$(conda shell.bash hook)"
conda activate ldpred_funct_new

rm $outCoord
#cd /work-zfs/abattle4/surya/tools/LDpred-funct
python /work-zfs/abattle4/surya/tools/LDpred-funct/ldpredfunct.py \
    --gf=${plinkfile} \
    --pf=${phenotype} \
    --ssf=${statsfile} --FUNCT_FILE=${functfile} --N=${N} \
    --H2=${h2} \
    --K=${K} \
    --coord=${outCoord} \
    --posterior_means=${outLdpredfunct} \
    --out=${outValidate} > ${outValidate}.log


# Note: Make sure that the file ${outCoord}  does not exists already.
# rm ${outCoord}
#python /work-zfs/abattle4/surya/tools/LDpred-funct/ldpredfunct.py --gf=${plinkfile} --pf=${phenotype} --FUNCT_FILE=${functfile}  --coord=${outCoord} --ssf=${statsfile} --N=${N} --posterior_means=${outLdpredfunct}  --H2=${h2} --out=${outValidate} > ${outValidate}.log
