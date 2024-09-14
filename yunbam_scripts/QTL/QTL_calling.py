import sys
import time
import os
import pdb
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from statsmodels.regression import linear_model as sm
from statsmodels.stats.multitest import multipletests
from scipy import stats

#### Compute QTL

# Genotype data: gt_numerical_dat
# Peak data: RPKM_dat
# Genotype data weights: post_pp_dat


def compute_QTL_gti_peaki(datapoint):

    [peak_sample, gt_sample, weight_sample] = datapoint
    valid_samples = np.where(gt_sample!= -1)[0]
    
    y = np.array(peak_sample[valid_samples])
    y = y.astype(float)
    x = np.array(gt_sample[valid_samples])
    x_weights = np.array(weight_sample[valid_samples])
   
    x = sm.add_constant(x) 
    wls_model = sm.WLS(y, x, weights = x_weights)
    results = wls_model.fit()

    return results.pvalues[1]



def compute_QTL_peaki(WINDOW, peaki, gt_numerical_dat, post_pp_dat, samples):
    QTL_result_peaki = []
    [start, end] = [peaki['START'], peaki['END']]
    
    SNPs_close = np.where((np.array(gt_numerical_dat['POS']) < int(end + WINDOW)) * (np.array(gt_numerical_dat['POS']) > int(start - WINDOW)))[0]
   
    bb = np.array(gt_numerical_dat.iloc[SNPs_close][samples])
    cc = np.array(post_pp_dat.iloc[SNPs_close][samples])
    aa = np.repeat(np.array(peaki[samples])[np.newaxis], len(bb), axis=0)

    datapoints = zip(aa,bb,cc)
    try:
    	pvalues = list(map(compute_QTL_gti_peaki, datapoints)) 
    except:
        pdb.set_trace()
    QTL_result_peaki = pd.DataFrame({"PeakID": peaki['PEAK'], "Start": peaki['START'], "End": peaki['END'], "CHR_POS": gt_numerical_dat.iloc[SNPs_close]['CHR_POS'], "POS": gt_numerical_dat.iloc[SNPs_close]['POS'], "P-value": pvalues})

    return QTL_result_peaki



def compute_QTLs(chromosome, WINDOW, peak_df, genotype_df, save_dir, saveSuffix = '', weight_df = None, use_GT = False):

    samples1 = [x for x in peak_df.columns if x.startswith('HG')]
    samples2 = [x for x in genotype_df.columns if x.startswith('HG')]
    samples = np.intersect1d(samples1, samples2)
    print('Run QTL for %d samples' % len(samples))

    if weight_df is None:
        weight_df = genotype_df.copy()
        weight_df[samples] = 1

    if use_GT:
        gt = np.array([map(round, x) for x in np.array(genotype_df[samples])])    
        genotype_df[samples] = gt

    QTL_results = pd.DataFrame()
    print('Compute QTLs for %d peaks ...' % len(peak_df))
    start = time.time()
    for p in range(len(peak_df)):
        token = compute_QTL_peaki(WINDOW, peak_df.iloc[p], genotype_df, weight_df, samples)
        QTL_results = QTL_results.append(token)
        if (p+1) % 500 == 0:
            print('    %d peaks finished' % (p+1))
    end = time.time()
    print('Tested %d associations for chromosome%d' % (len(QTL_results), chromosome))
    print('    Used %f miniutes' % ((end - start)/60))
    print("") 
    QTL_results = QTL_results.sort_values('P-value')
    QTL_results.to_csv('%s/CHR%d_caQTLs_WINDOW_%skb%s.txt' % (save_dir, chromosome, str(WINDOW/1000.0), saveSuffix), sep='\t', index = False)

    return QTL_results



def compute_QTL_peakLevel(WINDOW, save_dir, saveSuffix = ''):

    try:
        fn2 = '%s/caQTLs_WINDOW_%skb%s_peakLevel.txt' % (save_dir, str(WINDOW/1000.0), saveSuffix)
        df_peakLevel = pd.read_csv(fn2, sep='\t')
        df_peakLevel = df_peakLevel[df_peakLevel['BH_pvalue'] < 0.05]
        return df_peakLevel[['PeakID', 'CHR_POS', 'Peak_level_Pvalue', 'BH_pvalue']]

    except:
        print('Read in QTL data by chromosome')

    df_peakLevel = pd.DataFrame()
    for chromosome in range(1,23):
        fn1 = '%s/CHR%d_caQTLs_WINDOW_%skb%s.txt' % (save_dir, chromosome, str(WINDOW/1000.0), saveSuffix)
        try:
            dfi_peakLevel = pd.read_csv('%s_peakLevel' % fn1, sep='\t')
        except:
            dfi = pd.read_csv(fn1, sep='\t')
            dfi = dfi[~dfi['P-value'].isnull()].reset_index(drop = True)
            dfi_peakLevel = dfi.loc[dfi.groupby('PeakID')["P-value"].idxmin()]
            peak_test_number = pd.DataFrame(dfi.groupby('PeakID').size())
            peak_test_number.columns = ['Test']
            dfi_peakLevel = peak_test_number.merge(dfi_peakLevel, on = ['PeakID'])
            dfi_peakLevel['Peak_level_Pvalue'] = [np.min([1,x]) for x in np.array(dfi_peakLevel['P-value']) * np.array(dfi_peakLevel['Test'])]
            dfi_peakLevel.to_csv('%s_peakLevel' % fn1, sep='\t', index = False)
        df_peakLevel = df_peakLevel.append(dfi_peakLevel)

    df_peakLevel = df_peakLevel[['PeakID','CHR_POS', 'Peak_level_Pvalue', 'Start','End']]
    df_peakLevel['BH_pvalue'] = multipletests(df_peakLevel['Peak_level_Pvalue'], method='fdr_bh')[1]
    df_peakLevel = df_peakLevel.sort_values('BH_pvalue')

    fn2 = '%s/caQTLs_WINDOW_%skb%s_peakLevel.txt' % (save_dir, str(WINDOW/1000.0), saveSuffix)
    df_peakLevel.to_csv(fn2,  sep='\t', index = False)

    return df_peakLevel[['PeakID', 'CHR_POS', 'Peak_level_Pvalue', 'BH_pvalue']]



def compute_QTL_Qvector(WINDOW, save_dir, saveSuffix = ''):

    try:
        fn2 = pd.read_csv('%s/caQTLs_WINDOW_%skb%s_peakLevel_qvector.txt' % (save_dir, str(WINDOW/1000.0), saveSuffix), nrows=10)
        return

    except:
        print('Read in QTL data by chromosome')

    df_peakLevel = pd.DataFrame()
    for chromosome in range(1,23):
        fn1 = '%s/CHR%d_caQTLs_WINDOW_%skb%s.txt' % (save_dir, chromosome, str(WINDOW/1000.0), saveSuffix)
        try:
            dfi_peakLevel = pd.read_csv('%s_peakLevel' % fn1, sep='\t')
        except:
            dfi = pd.read_csv(fn1, sep='\t')
            dfi = dfi[~dfi['P-value'].isnull()].reset_index(drop = True)
            pdb.set_trace()
            dfi_peakLevel = dfi.loc[dfi.groupby('PeakID')["P-value"].idxmin()]
            peak_test_number = pd.DataFrame(dfi.groupby('PeakID').size())
            peak_test_number.columns = ['Test']
            dfi_peakLevel = peak_test_number.merge(dfi_peakLevel, on = ['PeakID'])
            dfi_peakLevel['Peak_level_Pvalue'] = [np.min([1,x]) for x in np.array(dfi_peakLevel['P-value']) * np.array(dfi_peakLevel['Test'])]
            dfi_peakLevel.to_csv('%s_peakLevel' % fn1, sep='\t', index = False)
        df_peakLevel = df_peakLevel.append(dfi_peakLevel)

    df_peakLevel = df_peakLevel[['PeakID','CHR_POS', 'Peak_level_Pvalue', 'Start','End']]
    df_peakLevel['BH_pvalue'] = multipletests(df_peakLevel['Peak_level_Pvalue'], method='fdr_bh')[1]
    df_peakLevel = df_peakLevel.sort_values('BH_pvalue')

    fn2 = '%s/caQTLs_WINDOW_%skb%s_peakLevel.txt' % (save_dir, str(WINDOW/1000.0), saveSuffix)
    df_peakLevel.to_csv(fn2,  sep='\t', index = False)

    return df_peakLevel[['PeakID', 'CHR_POS', 'Peak_level_Pvalue', 'BH_pvalue']]

