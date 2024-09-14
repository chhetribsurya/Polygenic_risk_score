import sys
import time
import os
import pdb
import numpy as np
import pandas as pd
from collections import Counter

def readin_peak_samples(PEAK_dir):
    peak_dat = pd.read_csv('%s/peak_by_sample_matrix_RPKM.txt' % (PEAK_dir), sep=' ', nrows=10)
    samples = [x for x in peak_dat.columns if x.startswith('HG')]
    return samples


def read_in_peaks_Zenodo(PEAK_dir, chromosome):
    peak_dat = pd.read_csv('%s/chromosome%d_corrected_fpkm.txt' % (PEAK_dir, chromosome), sep='\t')
    peak_dat.columns = ['PEAK'] + list(peak_dat.columns[1:])

    peak_loc = pd.read_csv('%s/chromosome%d_loc.bed' % (PEAK_dir, chromosome), sep='\t')
    peak_loc.columns = ['PEAK', 'CHR', 'START', 'END']

    peak_dat = peak_dat.merge(peak_loc, on = 'PEAK')

    samples = [x for x in peak_dat.columns if x.startswith('HG')]
    return [peak_dat, samples]



def convert_gt_to_number(arri):
    '''
    convert genotype type A/A to 0, 1, 2 ...
    '''
    arri_idx = np.where(arri!='0')[0]
    nts = np.unique([a for b in [x.split('/') for x in arri[arri_idx]] for a in b])

    # remove the sites with more than two alleles, keep biallelic variants
    if len(nts) > 2:
        return np.ones(len(arri)) * (-1)
    nts_dict = {}
    i = 0
    for ni in nts:
        nts_dict[ni] = i
        i += 1
    numeric_gt = [np.sum([nts_dict[a] for a in x.split('/')]) for x in arri[arri_idx]]
    numeric_gt_arr = np.ones(len(arri)) * (-1)
    numeric_gt_arr[arri_idx] = numeric_gt
    return numeric_gt_arr


def obtain_numerical_gt(gt_dat, samples):
    'deprecated'
    gt_numerical_dat = []
    print('In total %d variants are converted' % len(gt_dat))
    start = time.time()
    gt_numerical_dat = list(map(convert_gt_to_number, np.array(gt_dat[samples])))
    end = time.time()
    print('    Used %f s' % (end - start))
    
    gt_numerical_dat = pd.DataFrame(gt_numerical_dat)

    # remove rows with more than two alleles
    validQTLsnps = np.where([len(set(x))>1 for x in np.array(gt_numerical_dat)])[0]
    gt_numerical_dat = gt_numerical_dat.iloc[validQTLsnps].reset_index(drop=True)
    gt_dat = gt_dat.iloc[validQTLsnps].reset_index(drop=True)
 
    gt_numerical_dat.columns = samples
    gt_numerical_dat['CHR_POS'] = gt_dat['CHR_POS']
    gt_numerical_dat['CHR'] = gt_dat['CHR']
    gt_numerical_dat['POS'] = gt_dat['POS']
    gt_numerical_dat = gt_numerical_dat[['CHR_POS', 'CHR', 'POS']  + list(samples)]
   
    gt_numerical_dat = gt_numerical_dat.iloc[np.where(~gt_numerical_dat['CHR_POS'].duplicated())[0]]
 
    return gt_numerical_dat
    


def readin_genotype(Genotype_dir, chromosome, samples, suffix = ''):
    ## Read in Genotpye
    gt_dat = pd.read_csv('%s/dosage_by_sample_matrix_chr%d%s.txt' % (Genotype_dir, chromosome, suffix), sep=' ', low_memory=False)
    [gt_numerical_dat, numbers] = qc_genotype_dat(gt_dat, samples)
    return [gt_numerical_dat, numbers]



def read_in_1k_variants(chromosome):
    fn = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR//Genotype/ALL.chr%d.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.maf005.variants.txt' % chromosome
    wgs_1k = pd.read_csv(fn, usecols = [0], sep=' ')
    wgs_1k_snps = np.array(wgs_1k['CHR_POS'])
    return wgs_1k_snps



def read_in_WGS_GT(WGS_dir, chromosome, samples_peaks = None):
    WGS_fn = '%s/ALL.chr%d.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.GBRsamples.maf005.recode.GT.FORMAT' % (WGS_dir, chromosome)
    WGS_result = pd.read_csv(WGS_fn, comment = '$', sep='\t', low_memory=False)
    WGS_result = WGS_result.drop_duplicates()

    WGS_result.index = WGS_result[['CHROM', 'POS']].apply(lambda x: '_'.join((str(x[0]), str(x[1]))), axis=1)
    WGS_result = WGS_result.replace('./.', '-1/-1')
    WGS_result = WGS_result.replace(0, '-1/-1')

    WGS_result.columns = ['CHR', 'POS'] + list(WGS_result.columns[2:])
    WGS_result.index.name = 'CHR_POS'
    WGS_result = WGS_result.reset_index()

    samples = [x for x in WGS_result.columns if x.startswith('HG')]
    if samples_peaks is not None:
    	samples = list(np.intersect1d(samples_peaks, samples))
    for s in samples:
        WGS_result[s] = [np.sum(list(map(int, x.split('|')))) for x in np.array(WGS_result[s])]

    [WGS_dat, numbers] = qc_genotype_dat(WGS_result, samples)

    WGS_dat = WGS_dat[~WGS_dat['CHR_POS'].duplicated()].reset_index(drop=True) 

    return [WGS_dat, samples, numbers]



def qc_genotype_dat(df, samples, MAC = 3, snps = None):

    df = df.copy()
    df['POS'] = list(map(str, np.array(df['POS'])))
    df = df.iloc[np.where([x.isdigit() for x in np.array(df['POS'])])[0]]
    df = df.copy()
    df['POS'] = list(map(int, np.array(df['POS'])))

    statistic = [len(df)]
    df = df[list(df.columns[:3]) + samples]

    if snps is not None:
        df = df.set_index('CHR_POS').loc[np.intersect1d(snps, df['CHR_POS'])]
        df['CHR_POS'] = df.index
        df = df[['CHR_POS', 'CHR', 'POS']  + list(samples)]
        df = df.reset_index(drop = True)
        print('Remove variants with not exist in 1000 Genome Project with MAF >= 0.05: %d --> %d' % (statistic[0], len(df)))
        statistic.append(len(df))

    # remove variants with one genotype
    validQTLsnps = np.where([len(set(np.round(x[x!=-1]))) > 1 for x in np.array(df[samples])])[0]
    print('Remove variants with one genotype: %d --> %d' % (len(df), len(validQTLsnps)))
    statistic.append(len(validQTLsnps))
    df = df.iloc[validQTLsnps].reset_index(drop=True)

    # Filter on number of each allele > 3
    valid_snps = np.where([np.sum(x[x!=-1]) >= MAC for x in np.array(df[samples])])[0]
    print('Remove variants with minor allele count < %d: %d --> %d' % (MAC, len(df), len(valid_snps)))
    statistic.append(len(valid_snps))
    df = df.iloc[valid_snps].reset_index(drop=True)

    # Filter on number of each allele < 3
    valid_snps = np.where([np.sum(2-x[x!=-1]) >= MAC for x in np.array(df[samples])])[0]
    print('Remove variants with major allele count < %d: %d --> %d' % (MAC, len(df), len(valid_snps)))
    statistic.append(len(valid_snps))
    df = df.iloc[valid_snps].reset_index(drop=True)


    return [df, statistic] 



def achieve_ll(rowi):
    '''
    Map PL scores back to likelihood
    '''
    pl_scores = [list(map(int, str(x).split(','))) for x in rowi]
    pl_scores = [np.array(x) for x in pl_scores]
    pl_scores = [x[x<10] for x in pl_scores]

    post_pp = [list(map(lambda x: 1/10 ** x, pl_s)) for pl_s in pl_scores]    
    max_post_pp = [np.max(np.array(x)/float(np.sum(x))) for x in post_pp]
    
    return max_post_pp


def derive_ll(info_dat, samples):
    info_arr = np.array(info_dat[samples])
    idx_score = np.where(info_arr!='0')
    post_pp_arr = achieve_ll(info_arr[idx_score])

    post_pp = np.zeros(info_arr.shape)
    post_pp[idx_score] = post_pp_arr
    
    post_pp_dat = pd.DataFrame(post_pp)
    post_pp_dat.columns = samples
    post_pp_dat['CHR_POS'] = info_dat['CHR_POS']
    post_pp_dat['CHR'] = info_dat['CHR']
    post_pp_dat['POS'] = info_dat['POS']
    post_pp_dat = post_pp_dat[['CHR_POS', 'CHR', 'POS']  + samples]
    
    return post_pp_dat




def readin_genotype_info(gt_dat, VCF_dir, chromosome, samples):
    ## Read in Genotpye INFO
    info_dat = pd.read_csv('%s/gt_info_by_sample_matrix_chr%d_atac.txt' % (VCF_dir, chromosome), sep=' ', low_memory=False)
    info_dat = info_dat[[x for x in info_dat.columns if 'Unnamed' not in x]]
    info_dat = info_dat.set_index('CHR_POS').loc[gt_dat['CHR_POS']].reset_index()
    info_dat = info_dat.replace(0, '0')
   
    # match sample order in genotype matrix with the peak matrix
    info_dat = info_dat[['CHR_POS', 'CHR', 'POS'] + samples]
 
    print('Derive posterior probability for the genotypes')
    start = time.time()
    post_pp_dat = derive_ll(info_dat=info_dat, samples=samples)
    end = time.time()
    print('    Used %f seconds' % (end - start))
        
    return post_pp_dat



if __name__ == "__main__":
    CHROMOSOME = int(sys.argv[1])
    GT_DIR = 'minDP2'
    SUFFIX = '_minDP2'

    PEAK_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie/Peaks'
    WGS_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/Genotype'

    # root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_subsample_0.5'
    root_dir = '/work-zfs/abattle4/heyuan/Variant_calling/datasets/GBR/ATAC_seq/alignment_bowtie'
    Genotype_dir = '%s/Called_GT/%s' % (root_dir, GT_DIR)
    VCF_dir = '%s/VCF_files' % root_dir
    QTL_dir = '%s/QTLs' % root_dir
    if not os.path.exists(QTL_dir):
        os.makedirs(QTL_dir)

    ## read in data
    SAMPLES = readin_peak_samples(PEAK_dir)
    GT_DAT = readin_genotype(Genotype_dir = Genotype_dir, chromosome=CHROMOSOME, samples=SAMPLES)
    WEIGHT_DAT = readin_genotype_info(gt_dat=GT_DAT, VCF_dir = VCF_dir, chromosome=CHROMOSOME, samples=SAMPLES)


    [WGS_DAT, SAMPLES_WGS] = read_in_WGS_GT('1k_genome_chr%d' % CHROMOSOME, WGS_dir, SAMPLES, snps = np.array(GT_DAT['CHR_POS']))
    SAMPLES = list(np.intersect1d(SAMPLES, SAMPLES_WGS))

    [WGS_dat, _] = read_in_WGS_GT('1k_genome_chr%d' % CHROMOSOME, WGS_dir, SAMPLES)
    WGS_DAT_all = obtain_numerical_gt(WGS_dat_all, SAMPLES)


    save_matrix = False
    if save_matrix:
        save_dir = os.path.join(Genotype_dir, 'Save_Matrix')
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        GT_DAT[['CHR_POS'] + SAMPLES].to_csv('%s/called_genotypes_chromosome%d%s.txt' % (save_dir, CHROMOSOME, SUFFIX), sep='\t', index = False)
        GT_DAT[['CHR_POS','CHR', 'POS']].to_csv('%s/called_genotypes_chromosome%d%s_loc.bed' % (save_dir, CHROMOSOME, SUFFIX), sep='\t', index=False)

        WEIGHT_DAT[SAMPLES].to_csv('%s/called_genotypes_chromosome%d%s_weights.txt' % (save_dir, CHROMOSOME, SUFFIX), sep='\t')

        WGS_DAT[['CHR_POS'] + SAMPLES].to_csv('%s/called_genotypes_chromosome%d%s_realGT.txt' % (save_dir, CHROMOSOME, SUFFIX), sep='\t', index = False)
        WGS_DAT[['CHR_POS','CHR', 'POS']].to_csv('%s/called_genotypes_chromosome%d%s_loc_realGT.bed' % (save_dir, CHROMOSOME, SUFFIX), sep='\t', index=False)

        WGS_DAT_all[['CHR_POS'] + SAMPLES].to_csv('%s/called_genotypes_chromosome%d%s_realGT_all.txt' % (save_dir, CHROMOSOME, SUFFIX), sep='\t', index = False)
        WGS_DAT_all[['CHR_POS','CHR', 'POS']].to_csv('%s/called_genotypes_chromosome%d%s_loc_realGT_all.bed' % (save_dir, CHROMOSOME, SUFFIX), sep='\t', index=False)

    

