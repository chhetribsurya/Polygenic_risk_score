#!/usr/bin/env python

""" enrichmentAnalysis using scipy.stats.Fishers_test """

# cat kinases.txt | enrichmentAnalysis -10 4445 151 > /home/sfischer/junk/goenrich

import sys, traceback 
import scipy as sp # Open Source Library of Scientific Tools, including numpy
from scipy import stats as sciStats # stats module of scipy/shortcut for ease of access contains fisher's exact test
import statsmodels.sandbox.stats.multicomp as sm # statistical models package for python; contains multiple hypothesis testing
import pandas as pd # Python Data Analysis Library; provides high performance data structures; data frames and sorting

#=======================================================
# usage: prints usage message
#=======================================================

def usage():
    sys.stderr.write("enrichmentAnalysis:\n")
    sys.stderr.write("Compute an enrichment score for each of a set of terms provided\n")
    sys.stderr.write("on standard input.  Uses the Fisher's Exact Test.\n")
    sys.stderr.write("\nWrites results to standard output.\n")
    sys.stderr.write("\nUsage:  enrichmentAnalysis  pValueCutoff  backgroundTotal resultTotal\n")
    sys.stderr.write("\nWhere:\n")
    sys.stderr.write("\tpValueCutoff: terms with a p-value > cutoff not returned; expects value (0,1]\n")
    sys.stderr.write("\tbackgroundTotal: the total number of annotated entities in the background\n")
    sys.stderr.write("\tresultTotal: the total number of annotated entities in the result\n")
    sys.stderr.write("\nProvide input data on standard in, tab delimited with three columns:\n")
    sys.stderr.write("\ttermID, backgroundCount, resultCount\n")
    sys.stderr.write("(Additional trailing columns provided on input are passed through as trailing columns on output)\n")
    sys.stderr.write("\nWrites output to standard out.\n")
    sys.stderr.write("\tExcludes terms with p-value > cutoff.\n")
    sys.stderr.write("\tSorted by p-value in ascending order.\n")
    sys.stderr.write("\tTab delimited.\n")
    sys.stderr.write("\tIncludes the following (in order):\n")
    sys.stderr.write("\t\tFold Enrichment: ratio of two proportions:\n")
    sys.stderr.write("\t\t\tfraction of genes annotated by the term in result set\n")
    sys.stderr.write("\t\t\tto fraction of annotated genes in the background.\n")
    sys.stderr.write("\t\tOdds Ratio: statistic from the fisher's exact test\n")
    sys.stderr.write("\t\tpercent: percentage of genes annotated by this term in the result set,\n")
    sys.stderr.write("\t\tout of all annotated genes in the result set\n") 
    sys.stderr.write("\tp-value: pvalue from fisher's exact test\n")
    sys.stderr.write("\tBenjamini: Benjamini-Hochberg FDR\n")
    sys.stderr.write("\tBonferroni: Bonferroni adjusted p-value\n")
    sys.stderr.write("\tcolumns from input: termID, backgroundCount, resultCount, <trailing columns>\n\n")
    sys.exit(0)

#=======================================================
# toNum: wrapper for converting strings to numeric values
#=======================================================

def toNum(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

#=======================================================
# processArgs: processes command line arguments
#=======================================================

def processArgs(argv):
    if (len(argv) < 4):
        sys.stderr.write("ERROR: Provided " + str(len(argv) -1) + " arguments.  Expected 3.\n")
        usage()
        
    cutoff = toNum(argv[1])
    if (cutoff > 1 or cutoff <= 0):
        sys.stderr.write("ERROR: pValueCutoff must be >0 and <=1 \n")
        usage()

    return(cutoff, toNum(argv[2]), toNum(argv[3]))

#=======================================================
# fishersExactTest: uses scipy.stats fisher_exact to perform test
# returns pandas.DataFrame indexed on termID containing pvalue, oddsRatio,
# and percent entities annotated by term in the result set
#=======================================================

def fishersExactTest(data, rTotal, bTotal, threshold):
    fisherResults = {}
    for idx, row in data.iterrows():
        termID = row[0]
        baCount = row[1] # background annotated
        raCount = row[2] # result annotated
        
        bNotAnnotated = bTotal - baCount
        rNotAnnotated = rTotal - raCount
        percentAnnotated = (float(raCount) / float(rTotal)) * 100.0
        foldEnrichment = (float(raCount) / float(rTotal)) / (float(baCount) / float(bTotal))
        oddsRatio, pValue = sciStats.fisher_exact([[raCount, baCount],[rNotAnnotated, bNotAnnotated]], alternative = 'greater') # looking for over-enrichment
        if (pValue <= threshold) :
            fisherResults[idx] = {'pValue' : pValue, 'percent' : "{0:.2f}".format(percentAnnotated), 'oddsRatio' : "{0:.2f}".format(oddsRatio), 'foldEnrichment' : "{0:.2f}".format(foldEnrichment)}

    return (pd.DataFrame.from_dict(fisherResults, orient = 'index', dtype=float))
  
#=======================================================
# main
#=======================================================

def main(argv = None):

    if argv is None:
        argv = sys.argv
 
    pValueCutoff, backgroundTotalCount, resultTotalCount = processArgs(argv)
    
    try:
        data = pd.read_table(sys.stdin, header=None ) # read data to a pandas data frame 
      
        # Fisher Exact Test
        fisherResults = fishersExactTest(data, resultTotalCount, backgroundTotalCount, pValueCutoff)
        
        # Multiple Hypothesis Testing
        bonferroni = sm.multipletests(fisherResults['pValue'], method = 'bonferroni')
        benjamini = sm.multipletests(fisherResults['pValue'], method = 'fdr_bh')

        adjustedStats = pd.DataFrame({'bonferroni' : bonferroni[1], # second element of the tuple is the pvalues_corrected array
                                      'benjamini' : benjamini[1]}, index = fisherResults.index)
     
        # combine DataFrames based on indexes
        adjustedResults = pd.merge(fisherResults, adjustedStats, how="left", left_index = True, right_index = True) 
        results = pd.merge(adjustedResults, data, how = "left", left_index = True, right_index = True)

        # sort on p-value asc
        results.sort(columns = 'pValue', ascending = True, inplace = True)

        # print results.shape

        results.to_csv(sys.stdout, sep="\t", header = False, index = False)

    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        errorMessage = traceback.format_exception(exc_type, exc_value, exc_traceback)
        sys.stderr.write('Error: '.join(line for line in errorMessage)) 
 
if __name__ == "__main__":
    main()

