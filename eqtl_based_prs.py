
import numpy as np, pandas as pd
import pybedtools
import os, re
from glob import glob
from os.path import basename, join, splitext
from collections import defaultdict
import time
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy import stats



def perform_ldprune():
    """either use ld prune or LD clumping
        to eliminate LD effects
    """
    pass

def perform_ldclump():
    """either use ld prune or LD clumping
        to eliminate LD effects
    """
    pass

def eqtl_shared_variants():
    filter within 1MB
    pass

def weighted_beta_both():
    """apply func for generating shared Beta across
       population for prediction of population agnostic
       prs prediction
    """
    filter cohorts based on pval
    pass

def pop_agnostic_prs():
    """ uses weighted_beta_both func() inorder to calc
        population agnostic prs
    """
    pass

def ancestry_specific_prs():
    """ uses ancestry specific beta to calc population 
        specific prs
    """
    pass


aa_filepath = "/Users/suryachhetri/my_data/gtex/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_AFR_eQTL_all_associations"
ea_filepath = "/Users/suryachhetri/my_data/gtex/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations"
tissue = "Heart_Atrial"
tissue = "Adipose_Subcutaneous"
ldgenome_path = ""
chromrange = 1
chrom_X = False
pval = 10e-5 
output_path = "/Users/suryachhetri/Desktop/prs_project/datasets"


def pop_shared_variants(aa_filepath, ea_filepath, tissue, chromrange):
    """ Finds shared variants plus specfic variants between 
        the population cohorts without similar effects/eQTL info
    """
    ## Enter the chromosome no. range that you would like to analyse the data on. By default, it would take the autosomal 
    ## gene model from (chr1:chr22, excluding chrX & chrY), however, set chrom_XY=True for normal gene model analysis.
    chrom_list = []
    for chrom_no in range(chromrange):
        chrom_list.append("chr" + str(chrom_no + 1))
    if chrom_X:
        chrom_list.append("chrX")
        #chrom_list.append("chrY")

    # read chromwise for AFR and EUR:
    aa_filelist = glob(join(aa_filepath, str(tissue) + "*.parquet"))
    aa_dict = {basename(file).split(".")[4] :file for file in aa_filelist}
    ea_filelist = glob(join(ea_filepath, str(tissue) + "*.parquet"))
    ea_dict = {basename(file).split(".")[4] :file for file in ea_filelist}

    shared_variants_df = []
    aa_specific_df = []
    ea_specific_df = []
    for chrom in chrom_list:
        print("processing : {}".format(chrom))

        """ parse afr cohorts"""
        read_aa = pd.read_parquet(aa_dict.get(chrom), engine="fastparquet")
        read_aa.reset_index(inplace=True, drop=True)
        read_aa.drop_duplicates(["variant_id"], inplace=True)
        read_aa["new_col"] = read_aa["variant_id"].map(lambda x:x.split("_"))
        aa_df = pd.DataFrame(read_aa['new_col'].values.tolist(), columns=['chrom','start', 'ref', 'alt', 'assemb'])

        """ parse eur cohorts"""
        read_ea = pd.read_parquet(ea_dict.get(chrom), engine="fastparquet")
        read_ea.reset_index(inplace=True, drop=True)
        read_ea.drop_duplicates(["variant_id"], inplace=True)
        read_ea["new_col"] = read_ea["variant_id"].map(lambda x:x.split("_"))
        ea_df = pd.DataFrame(read_ea['new_col'].values.tolist(), columns=['chrom','start', 'ref', 'alt', 'assemb'])

        """ merge afr eur cohorts """
        merged_df = pd.merge(aa_df, ea_df, on=aa_df.columns.to_list(), how = "outer", indicator=True)
        shared_variants = merged_df[merged_df["_merge"] == "both"]
        aa_specific = merged_df[merged_df["_merge"] == "left_only"]
        ea_specific = merged_df[merged_df["_merge"] == "right_only"]

    print("combining shared variants for chroms...")
    shared_variants_all = pd.concat(shared_variants_df, ignore_index=True)
    print("combining AA specific variants for chroms...")
    aa_specific_all = pd.concat(aa_specific_df, ignore_index=True)
    print("combining AA specific variants for chroms...")
    ea_specific_all = pd.concat(ea_specific_df, ignore_index=True)
    return(shared_variants, a_specific, ea_specific)

# Variant segregation:
shared, afr_specific, eur_specific =  pop_shared_variants(aa_filepath, ea_filepath, tissue, 1)


def pop_shared_variants_betabased(self):
    """ 
    Finds shared variants and specfic variants between 
    population cohorts without similar effects/eQTL info

    """
    ## Enter the chromosome no. range that you would like to analyse the data on. By default, it would take the autosomal 
    ## gene model from (chr1:chr22, excluding chrX & chrY), however, set chrom_XY=True for normal gene model analysis.
    
    chrom_list = []
    for chrom_no in range(chromrange):
        chrom_list.append("chr" + str(chrom_no + 1))
    if chrom_X:
        chrom_list.append("chrX")
        #chrom_list.append("chrY")

    # read chromwise for AFR and EUR:
    aa_filelist = glob(join(aa_filepath, str(tissue) + "*.parquet"))
    aa_dict = {basename(file).split(".")[4] :file for file in aa_filelist}
    ea_filelist = glob(join(ea_filepath, str(tissue) + "*.parquet"))
    ea_dict = {basename(file).split(".")[4] :file for file in ea_filelist}

    shared_variants_set1 = []
    shared_variants_set2 = []
    aa_specific_df = []
    ea_specific_df = []
    for chrom in chrom_list:
        print("processing : {}".format(chrom))

        """ parse afr cohorts """
        read_aa = pd.read_parquet(aa_dict.get(chrom), engine="fastparquet")
        #read_aa.to_csv(join(output_path, "Adipose_Subcutaneous.v8.AFR.eqtl_allpairs.chr1.txt"), sep="\t", header=True)
        read_aa.reset_index(inplace=True, drop=True)
        sorted_aa = read_aa.sort_values(["variant_id", "pval_nominal"], ascending=True).reset_index(drop=True)
        sorted_aa.drop_duplicates(["variant_id"], inplace=True, keep="first") #drop based on most sig pval

        """ parse eur cohorts"""
        read_ea = pd.read_parquet(ea_dict.get(chrom), engine="fastparquet")
        #read_ea.to_csv(join(output_path, "Adipose_Subcutaneous.v8.EUR.eqtl_allpairs.chr1.txt"), sep="\t", header=True)
        read_ea.reset_index(inplace=True, drop=True)
        sorted_ea = read_ea.sort_values(["variant_id", "pval_nominal"], ascending=True).reset_index(drop=True)
        sorted_ea.drop_duplicates(["variant_id"], inplace=True, keep="first") #drop based on most sig pval

        """ merge afr eur cohorts """
        merged_df = pd.merge(sorted_aa, sorted_ea, on=["variant_id", "phenotype_id"], how = "outer", suffixes=["_AA", "_EA"], indicator=True)
        shared_variants = merged_df[merged_df["_merge"] == "both"] #filter based on sig pval 0.0001 of at least in one cohort
        print("total shared variants across population without thresholds : {}".format(shared_variants.shape[0]))
        shared_vars_thres = shared_variants.loc[(shared_variants['pval_nominal_AA'] <= pval) | (shared_variants['pval_nominal_EA'] <= pval)]
        aa_specific = merged_df[merged_df["_merge"] == "left_only"]
        aa_specific_thresh = aa_specific.loc[aa_specific['pval_nominal_AA'] <= pval]
        ea_specific = merged_df[merged_df["_merge"] == "right_only"]
        ea_specific_thresh = ea_specific.loc[ea_specific['pval_nominal_EA'] <= pval]
        
        """Basic summary stats"""
        print('''Total unique variants detected across cohorts ::
        European cohorts : {0}
        African cohorts : {1}'''.format(sorted_aa.shape[0], sorted_ea.shape[0]))
        print('''\n\nVariants dist without thresholds ::
        Shared (at least one pop thresh) : {0}
        AFR specific : {1}
        EUR specific : {2}'''.format(shared_variants.shape[0],aa_specific.shape[0],ea_specific.shape[0]))
        print('''\nVariants with pval {0} thresh ::
        Shared : {1}
        AFR specific : {2}
        EUR specific : {3}'''.format(pval, shared_vars_thres.shape[0],aa_specific_thresh.shape[0],ea_specific_thresh.shape[0]))

        """calculate correlation"""
        shared_vars_thres["logpval_nominal_AA"] = -(np.log10(shared_vars_thres["pval_nominal_AA"]))
        shared_vars_thres["logpval_nominal_EA"] = -(np.log10(shared_vars_thres["pval_nominal_EA"]))
        
        """calculate variant-gene pairwise slope corr(r) and coeff of var(rsquare) for each gene"""
        def scipystats_linregress(df): return stats.linregress(df["slope_AA"], df["slope_EA"])[2] #slope,intercept,r_value,p_value,std_err
        stats_df = shared_vars_thres.groupby(["phenotype_id"]).apply(scipystats_linregress)
        stats_df = stats_df.reset_index(name="beta_correlation")
        stats_filt_df = stats_df[~(stats_df["beta_correlation"] == 0)].reset_index(drop=True)
        stats_filt_df["rsquare"] = stats_filt_df["beta_correlation"]**2
        corr_median = stats_filt_df["beta_correlation"].median()
        rsquare_median = stats_filt_df["rsquare"].median()
        print("total genes, filtered genes :  {}, {}".format(stats_df.shape[0], stats_filt_df.shape[0]))

        def pandas_linregress_corr(df):
            corr = df.groupby("phenotype_id")[['slope_AA','slope_EA']].corr()
            df_corr = corr.iloc[0::2,-1]
            corr_df = df_corr.reset_index()
            nans = lambda df1: df1.loc[df1.isnull().any(axis=1)]
            nancorr_rows = nans(corr_df)
            corr_df.dropna(inplace=True)
            corr_df.columns = ["phenotype_id", "slope_level", "beta_correlation"]
            corr_df = corr_df.loc[:, ["phenotype_id", "beta_correlation"]]
            return corr_df

        """ calculate correlation using pandas"""
        corr_df = pandas_linregress_corr(shared_vars_thres)

        """ plot histogram for correlation distribution """
        fig, ax = plt.subplots(1,1) # set up figure and axes; concise than : fig = plt.figure(); ax = fig.add_subplot(111)
        ax.hist(stats_filt_df["beta_correlation"], color="red", bins=20)
        plt.grid(True)
        plt.axvline(corr_median, color='b', linestyle='dashed', linewidth=2)
        plt.xlabel('corr(beta_AA vs beta_EA)  |  SharedVars(n = {})'.format(str(shared_vars_thres.shape[0])))
        plt.title('Beta Correlation for ' + str(chrom) + "(n = {} genes)".format(str(stats_filt_df.shape[0])))
        median = "Median : {}".format(str(round(corr_median, 2)))
        median_text = AnchoredText(median, loc=2)
        ax.add_artist(median_text)
        plt.savefig('/Users/suryachhetri/Desktop/prs_project/{}_betaStats_corr.pdf'.format(chrom)) #saves current figure
        plt.close()

        fig, ax = plt.subplots(1,1)
        ax.hist(stats_filt_df["rsquare"], color="red", bins=20)
        plt.grid(True)
        plt.axvline(rsquare_median, color='b', linestyle='dashed', linewidth=2)
        plt.xlabel('coeff(beta_AA vs beta_EA)  |  SharedVars(n = {})'.format(str(shared_vars_thres.shape[0])))
        plt.title('Beta Rsquare for ' + str(chrom) + "(n = {} genes)".format(str(stats_filt_df.shape[0])))
        median = "Median : {}".format(str(round(rsquare_median, 2)))
        median_text = AnchoredText(median, loc=2)
        ax.add_artist(median_text)
        plt.savefig('/Users/suryachhetri/Desktop/prs_project/{}_betaStats_coeff.pdf'.format(chrom)) #saves current figure
        plt.close()
        
        """ merge stats filt df for retention of variants and tss info"""
        shared_df = shared_vars_thres.reset_index(drop=True)
        merged_df1 = pd.merge(stats_filt_df, shared_df, on=["phenotype_id"])
        select_cols = ["phenotype_id", "variant_id", "tss_distance_AA", "tss_distance_EA", "slope_AA", "slope_EA",
                        "logpval_nominal_AA", "logpval_nominal_EA", "beta_correlation", 'rsquare']
        shared_var_set1 = merged_df1.loc[:,select_cols]
        shared_var_set1[["tss_distance_AA", "tss_distance_EA"]] = shared_var_set1[["tss_distance_AA", "tss_distance_EA"]].astype(int)

        #filtered merged datasets: # coeff of var of at least 0.5
        stats_05 = stats_filt_df[stats_filt_df["rsquare"] >= 0.5]
        merged_df2 = pd.merge(stats_05, shared_df, on=["phenotype_id"])
        shared_var_set2 = merged_df2.loc[:,select_cols]
        shared_var_set2[["tss_distance_AA", "tss_distance_EA"]] = shared_var_set1[["tss_distance_AA", "tss_distance_EA"]].astype(int)
        

        """ find variant to highlight """
        pheno_idx = shared_var_set1.groupby("phenotype_id").size().sort_values().idxmax()
        #pheno_idx = "ENSG00000001460.17"
        variant_df = shared_var_set1[shared_var_set1["phenotype_id"] == pheno_idx]
        
        # generate scatter plot:
        cm = plt.cm.get_cmap('RdYlBu_r') #["autumn", "viridis"]
        plt.subplot(211)
        plt.scatter(x=variant_df["tss_distance_AA"], y=variant_df["logpval_nominal_AA"], 
                    c=variant_df["slope_AA"], cmap=cm, s=4, alpha=0.8)
        plt.colorbar().set_label('Beta value')
        text = AnchoredText("AFR Cohort", loc=1)
        plt.gca().add_artist(text)
        plt.ylabel("-Log10(pval)")
        #plt.legend(loc=1, prop={'size': 8}) #requires label="" in plt.scatter()

        plt.subplot(212)
        plt.scatter(x=variant_df["tss_distance_EA"], y=variant_df["logpval_nominal_EA"], 
                    c=variant_df["slope_EA"], cmap=cm, s=4, alpha=0.8)
        plt.colorbar().set_label('Beta value')
        text = AnchoredText("EUR Cohort", loc=1)
        plt.gca().add_artist(text)
        plt.xlabel("TSS Dist(bp)")
        plt.ylabel("-Log10(pval)")
        #plt.legend(loc=1, prop={'size': 8}) #requires label="" in plt.scatter()
        
        # super title
        plt.suptitle('''Shared Variants Dist (n = {0})
        {1}'''.format(variant_df.shape[0], str(pheno_idx)), fontsize=8, y=0.98, ha='center')
        plt.savefig('/Users/suryachhetri/Desktop/prs_project/{}_variant_scatter.pdf'.format(pheno_idx)) #saves current figure
        plt.close()

    return(shared_var_set1, shared_var_set2)



###################
""" Custom spine bounds"""

# https://matplotlib.org/3.1.1/gallery/ticks_and_spines/spines_bounds.html
# Fixing random state for reproducibility
np.random.seed(19680801)

x = np.linspace(0, 2*np.pi, 50)
y = np.sin(x)
y2 = y + 0.1 * np.random.normal(size=x.shape)

fig, ax = plt.subplots()
ax.plot(x, y)
ax.plot(x, y2)

# set ticks and tick labels
ax.set_xlim((0, 2*np.pi))
ax.set_xticks([0, np.pi, 2*np.pi])
ax.set_xticklabels(['0', r'$\pi$', r'2$\pi$'])

ax.set_ylim((-1.5, 1.5))
ax.set_yticks([-1, 0, 1]) # add yticks ax.set_yticks([-1, 0, 1, 1.5])

# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Only draw spine between the y-ticks
ax.spines['left'].set_bounds(-1, 1.5)

# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

plt.close()
#plt.show()


###################

# For repeatable "random" data
np.random.seed(0)
# Specify the mean and standard deviation for each mock data group
data_specs = [(2, 2), (7, 1), (4, 2.5), (10, 0.5), (5.5, 0.1)]

# Generate data and place into a pandas DataFrame
data = [np.random.normal(mu, sigma, 10) for mu, sigma in data_specs]
data = pd.DataFrame(data).T
data.columns = ['Group_%s' % n for n in range(1,6)]

fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(8,6))

y = data.mean()
y_all = data.values
x = np.arange(len(means))
error = data.std()

xlims = (-1, 5)
ylims = (-5, 15)
bar_ylims = (0, 15)

custom_lineplot(ax[0][0], x, y, error, xlims, ylims)
custom_scatterplot(ax[0][1], x, y, error, xlims, ylims)
custom_barchart(ax[1][0], x, y, error, xlims, bar_ylims, error_kw)
custom_boxplot(ax[1][1], x, y_all, error, xlims, ylims)

titles = ['Line Plot', 'Scatter Plot', 'Bar Chart', 'Box Plot']
xlabel = 'Group'
ylabel = 'Value ($units^2$)'
xticks = x
xticklabels = range(1,6)

for i, axes in enumerate(ax.flat):
    # Customize y ticks on a per-axes basis
    yticks = np.linspace(axes.get_ylim()[0], axes.get_ylim()[1], 5)
    yticklabels = yticks
    stylize_axes(axes, titles[i], xlabel, ylabel, xticks, yticks, xticklabels, yticklabels)
    
fig.tight_layout()

################################

# Fixing random state for reproducibility
np.random.seed(19680801)

x = np.linspace(0, 2*np.pi, 50)
y = np.sin(x)
y2 = y + 0.1 * np.random.normal(size=x.shape)

fig, ax = plt.subplots()
ax.plot(x, y)
ax.plot(x, y2)

# set ticks and tick labels
ax.set_xlim((0, 2*np.pi))
ax.set_xticks([0, np.pi, 2*np.pi])
ax.set_xticklabels(['0', r'$\pi$', r'2$\pi$'])
ax.set_ylim((-1.5, 1.5))
ax.set_yticks([-1, 0, 1])

# Only draw spine between the y-ticks
ax.spines['left'].set_bounds(-1, 1)
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

plt.show()



################################
# Listed colormap with categorical coloring on colorbar matplotlib 
# https://stackoverflow.com/questions/15908371/matplotlib-colorbars-and-its-text-labels

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

#discrete color scheme
cMap = ListedColormap(['white', 'green', 'blue','red'])

#data
np.random.seed(42)
data = np.random.rand(4, 4)
fig, ax = plt.subplots()
heatmap = ax.pcolor(data, cmap=cMap)

#legend
cbar = plt.colorbar(heatmap)
cbar.ax.set_yticklabels(['0','1','2','>3'])
cbar.set_label('# of contacts', rotation=270)

# put the major ticks at the middle of each cell
ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor=False)
ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor=False)
ax.invert_yaxis()

#labels
column_labels = list('ABCD')
row_labels = list('WXYZ')
ax.set_xticklabels(column_labels, minor=False)
ax.set_yticklabels(row_labels, minor=False)

plt.show()

#################################

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Create dummy dataframe, or load your own with pd.read_csv()

columns = ["sex", "age", "BMI", "smoke", "type"]
data = pd.DataFrame(np.array([[1,0,0,1,0], [23,16,94,18,24], [32, 26, 28, 23, 19], [0,1,1,1,0], [1,2,2,2,1]]).T, columns=columns)


x_col = "sex"
y_columns = ["age", "BMI", "smoke"]


for y_col in y_columns:

    figure = plt.figure
    ax = plt.gca()
    ax.scatter(data[x_col], data[y_col])
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title("{} vs {}".format(x_col, y_col))

    plt.legend()
    plt.show()

#################################

#sample of scatter plot:
import matplotlib.pyplot as plt

xyc = range(20)

plt.subplot(121)
plt.scatter(xyc[:13], xyc[:13], c=xyc[:13], s=35, vmin=0, vmax=20)
plt.colorbar()
plt.xlim(0, 20)
plt.ylim(0, 20)

plt.subplot(122)
plt.scatter(xyc[8:20], xyc[8:20], c=xyc[8:20], s=35, vmin=0, vmax=20)   
plt.colorbar()
plt.xlim(0, 20)
plt.ylim(0, 20)

plt.show()


##################################
Datavmax = max(max(CurrentsArray))
Datavmin = min(min(CurrentsArray))

plt.subplot(121)
plt.scatter(growthTarray, CuSearray, PFarray, CurrentsArray, vmin=Datavmin, vmax=Datavmax, alpha=0.75)
plt.colorbar()
plt.xlim(600,760)
plt.ylim(0,2.5)

plt.subplot(122)
plt.scatter(rf85growthTarray, rf85CuSearray, rf85PFarray, rf85CurrentsArray, vmin=Datavmin, vmax=Datavmax, alpha=0.75)
plt.colorbar()
plt.xlim(600,760)
plt.ylim(0,2.5)

plt.show()

##############################################

# https://stackoverflow.com/questions/42973223/how-share-x-axis-of-two-subplots-after-they-are-created/42974975
# The usual way to share axes is to create the shared properties at creation. Either

# fig=plt.figure()
# ax1 = plt.subplot(211)
# ax2 = plt.subplot(212, sharex = ax1)
# or

# fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
# Sharing the axes after they have been created should therefore not be necessary.

# However if for any reason, you need to share axes after they have been created (actually, using a different library which creates some subplots, like here, or sharing an inset axes might be a reason), there would still be a solution:

# Using

# ax1.get_shared_x_axes().join(ax1, ax2)



# draw scatter plot based on -log(pval)
# for 'ENSG00000280670.2'
# shared_vars_thres.groupby("phenotype_id")
#                            .size()
#                            .sort_values()
#                            .idxmax()

#############################################


In [1939]: shared_vars_thres.groupby("phenotype_id").size().sort_values()
Out[1939]:
phenotype_id
ENSG00000127472.10      1
ENSG00000215859.8       1
ENSG00000121940.15      1
ENSG00000215915.9       1
ENSG00000118729.11      1
                     ...
ENSG00000143452.15    281
ENSG00000224468.3     352
ENSG00000204138.12    425
ENSG00000162782.15    493
ENSG00000280670.2     522
Length: 546, dtype: int64

In [1940]: shared_vars_thres.groupby("phenotype_id").size().sort_values().idxmax()
Out[1940]: 'ENSG00000280670.2'

In [1941]: idx=shared_vars_thres.groupby("phenotype_id").size().sort_values().idxmax()

In [1942]: corr_df[corr_df["phenotype_id"] == idx]
Out[1942]:
          phenotype_id  beta_correlation
542  ENSG00000280670.2          0.968443

In [1943]: shared_vars_thres[shared_vars_thres["phenotype_id"] == idx]






# set up figure and axes
f, ax = plt.subplots(1,1)

ax.scatter(xa,ya, marker='o', s=20, c="lightgreen", alpha=0.9)
ax.scatter(xb,yb, marker='o', s=20, c="dodgerblue", alpha=0.9)
ax.scatter(xc,yc marker='o', s=20, c="firebrick", alpha=1.0)
ax.scatter(xd,xd,xd, marker='o', s=20, c="goldenrod", alpha=0.9)
line1 = Line2D(range(10), range(10), marker='o', color="goldenrod")
line2 = Line2D(range(10), range(10), marker='o',color="firebrick")
line3 = Line2D(range(10), range(10), marker='o',color="lightgreen")
line4 = Line2D(range(10), range(10), marker='o',color="dodgerblue")
plt.legend

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

# make some data
x = np.arange(10)
y = x

# set up figure and axes
f, ax = plt.subplots(1,1)

# loc works the same as it does with figures (though best doesn't work)
# pad=5 will increase the size of padding between the border and text
# borderpad=5 will increase the distance between the border and the axes
# frameon=False will remove the box around the text

anchored_text = AnchoredText("Test", loc=2)
ax.plot(x,y)
ax.add_artist(anchored_text)

plt.show()



shared_variants_corr = shared_vars_thres.groupby("phenotype_id")[['slope_AA','slope_EA']].corr()
#shared_variants_corr = shared_vars_thres.groupby("phenotype_id")[['slope_AA','slope_EA']]
#shared_vars_thres.groupby("phenotype_id").size().sort_values()

variants_corr = shared_variants_corr.iloc[0::2,-1]
corr_df = variants_corr.reset_index()
nans = lambda df: df.loc[df.isnull().any(axis=1)]
nancorr_rows = nans(corr_df)
corr_df.dropna(inplace=True)
corr_df.columns = ["phenotype_id", "slope_level", "beta_correlation"]
corr_df = corr_df.loc[:, ["phenotype_id", "beta_correlation"]]
#corr_df[(corr_df["beta_correlation"] < 1) & (corr_df["beta_correlation"] > 0.8)]
#corr_df[corr_df["beta_correlation"] == 1]
corr_median = corr_df["beta_correlation"].median()
print("filtered genes, dropped genes :  {}, {}".format(corr_df.shape[0], nancorr_rows.shape[0]))



        # #shared_variants_corr = shared_vars_thres.groupby("phenotype_id")[['logpval_nominal_AA','logpval_nominal_EA']].corr()
        # #shared_vars_thres[shared_vars_thres.isnull().any(axis=1)]

        # variants_corr = shared_variants_corr.iloc[0::2,-1]
        # corr_df = variants_corr.reset_index(drop=True)

        # """ select variants with same geneid """
        # merged_df1 = pd.merge(sorted_aa, sorted_ea, on=["variant_id"], how = "outer", suffixes=["_AA", "_EA"], indicator=True)

        # """ filter out variants with disparate geneid"""
        # shared_variants = merged_df1[merged_df1["_merge"] == "both"]
        # shared_variants_set2 = shared_variants2[shared_variants2["phenotype_id_AA"] == shared_variants1["phenotype_id_EA"]]
        # shared_variants_set2 = shared_variants1[~(shared_variants1["phenotype_id_AA"] == shared_variants1["phenotype_id_EA"])]
        # #shared_variants[~(shared_variants["phenotype_id_AA"] == shared_variants["phenotype_id_EA"])].iloc[:,[0,1,2,9,10,11]]

        """ calculate variant-gene pairwise pval_correlation for each gene"""
        shared_variants_corr = shared_variants_set1.groupby("phenotype_id_AA")[['logpval_nominal_AA','logpval_nominal_EA']].corr()


        # test_df = shared_variants_set1[shared_variants_set1["phenotype_id_AA"] == "ENSG00000281571.2"]
        # test_df.plot.scatter(y="logpval_nominal_AA", x="tss_distance_AA", color="blue")
        # plt.title('ENSG00000281571.2')
        # plt.savefig('/Users/suryachhetri/Desktop/scatter_pval_AA.pdf') #saves current figure
        
        # test_df.plot.scatter(y="logpval_nominal_EA", x="tss_distance_EA", color="red")
        # plt.title('ENSG00000281571.2')
        # plt.savefig('/Users/suryachhetri/Desktop/scatter_pval_EA.pdf') #saves current figure
        # plot_df.plot('scatter', y="logpval_nominal_EA", x="tss_distance_EA")

        # Get mean as standard deviation
        negative_corr = corr_df[corr_df["logpval_nominal_EA"] < 0]
        negative_corr["logpval_nominal_EA"].mean()
        neg_mean = negative_corr["logpval_nominal_EA"].mean()
        neg_std = negative_corr["logpval_nominal_EA"].std()
        negative_corr
        

        positive_corr = corr_df[corr_df["logpval_nominal_EA"] > 0]
        positive_corr["logpval_nominal_EA"].mean()
        pos_mean = positive_corr["logpval_nominal_EA"].mean()
        pos_std = positive_corr["logpval_nominal_EA"].std()
        positive_corr
        
        """ calculate variant-gene pairwise slopeSE_correlation for each gene"""
        shared_variants_corr2 = shared_variants_set1.groupby("phenotype_id_AA")[['slope_se_AA','slope_se_EA']].corr()
        variants_corr2 = shared_variants_corr2.iloc[0::2,-1]
        corr_df2 = variants_corr2.reset_index()

        negative_corr2 = corr_df2[corr_df2["slope_se_EA"] < 0]
        negative_corr2["slope_se_EA"].mean()
        neg_mean2 = negative_corr2["slope_se_EA"].mean()
        neg_std2 = negative_corr2["slope_se_EA"].std()
        negative_corr2

        positive_corr2 = corr_df2[corr_df2["slope_se_EA"] > 0]
        positive_corr2["slope_se_EA"].mean()
        pos_mean2 = positive_corr2["slope_se_EA"].mean()
        pos_std2 = positive_corr2["slope_se_EA"].std()
        positive_corr2

        """ calculate variant-gene pairwise slopeSE_correlation for each gene"""
        shared_variants_corr3 = shared_variants_set1.groupby("phenotype_id_AA")[['slope_AA','slope_EA']].corr()
        variants_corr3 = shared_variants_corr3.iloc[0::2,-1]
        corr_df3 = variants_corr3.reset_index()

        negative_corr3 = corr_df3[corr_df3["slope_EA"] < 0]
        negative_corr3["slope_EA"].mean()
        neg_mean3 = negative_corr3["slope_EA"].mean()
        neg_std3 = negative_corr3["slope_EA"].std()
        negative_corr3

        positive_corr3 = corr_df3[corr_df3["slope_EA"] > 0]
        positive_corr3["slope_EA"].mean()
        pos_mean3 = positive_corr3["slope_EA"].mean()
        pos_std3 = positive_corr3["slope_EA"].std()
        positive_corr3

        #histogram
        corr_df.hist()
        plt.axvline(pos_mean, color='b', linestyle='dashed', linewidth=2)
        plt.axvline(neg_mean, color='b', linestyle='dashed', linewidth=2)
        plt.xlabel('corr(pval_AA vs pvalEA)')
        plt.title('Pval Correlation')
        plt.savefig('/Users/suryachhetri/Desktop/genebased_pval_corr.pdf') #saves current figure

        corr_df2.hist(color="red")
        plt.axvline(pos_mean2, color='b', linestyle='dashed', linewidth=2)
        plt.axvline(neg_mean2, color='b', linestyle='dashed', linewidth=2)
        plt.title('slopeSE Correlation')
        plt.xlabel('corr(slopeSE_AA vs slopeSE_EA)')
        plt.savefig('/Users/suryachhetri/Desktop/genebased_slopeSE_corr.pdf') #saves current figure

        corr_df3.hist(color="grey")
        plt.axvline(pos_mean3, color='b', linestyle='dashed', linewidth=2)
        plt.axvline(neg_mean3, color='b', linestyle='dashed', linewidth=2)
        plt.title('slope Correlation')
        plt.xlabel('corr(slope_AA vs slope_EA)')
        plt.savefig('/Users/suryachhetri/Desktop//genebased_slope_corr.pdf') #saves current figure


    #print("combining shared variants set1 for chroms...")
    shared_variants_set1 = pd.concat(shared_variants_set1, ignore_index=True)
    #print("combining shared variants set2 for chroms...")
    shared_variants_set1 = pd.concat(shared_variants_set2, ignore_index=True)
    #print("combining AA specific variants for chroms...")
    aa_specific_all = pd.concat(aa_specific_df, ignore_index=True)
    #print("combining AA specific variants for chroms...")
    ea_specific_all = pd.concat(ea_specific_df, ignore_index=True)
    return(shared_variants, a_specific, ea_specific)

shared_variants[~(shared_variants["phenotype_id_AA"] == shared_variants["phenotype_id_EA"])].iloc[:,[0,1,2,9,10,11]]

def pop_shared_variants_pvalbased(self):
    """ Finds shared variants and specfic variants between 
        population cohorts without similar effects/eQTL info
    """
    ## Enter the chromosome no. range that you would like to analyse the data on. By default, it would take the autosomal 
    ## gene model from (chr1:chr22, excluding chrX & chrY), however, set chrom_XY=True for normal gene model analysis.
    chrom_list = []
    for chrom_no in range(chromrange):
        chrom_list.append("chr" + str(chrom_no + 1))
    if chrom_X:
        chrom_list.append("chrX")
        #chrom_list.append("chrY")

    # read chromwise for AFR and EUR:
    aa_filelist = glob(join(aa_filepath, str(tissue) + "*.parquet"))
    aa_dict = {basename(file).split(".")[4] :file for file in aa_filelist}
    ea_filelist = glob(join(ea_filepath, str(tissue) + "*.parquet"))
    ea_dict = {basename(file).split(".")[4] :file for file in ea_filelist}

    shared_variants_df = []
    aa_specific_df = []
    ea_specific_df = []
    for chrom in chrom_list:
        print("processing : {}".format(chrom))

        """ parse afr cohorts"""
        read_aa = pd.read_parquet(aa_dict.get(chrom), engine="fastparquet")
        read_aa.reset_index(inplace=True, drop=True)
        read_aa.drop_duplicates(["variant_id"], inplace=True)
        read_aa["new_col"] = read_aa["variant_id"].map(lambda x:x.split("_"))
        aa_df = pd.DataFrame(read_aa['new_col'].values.tolist(), columns=['chrom','start', 'ref', 'alt', 'assemb'])

        """ parse eur cohorts"""
        read_ea = pd.read_parquet(ea_dict.get(chrom), engine="fastparquet")
        read_ea.reset_index(inplace=True, drop=True)
        read_ea.drop_duplicates(["variant_id"], inplace=True)
        read_ea["new_col"] = read_ea["variant_id"].map(lambda x:x.split("_"))
        ea_df = pd.DataFrame(read_ea['new_col'].values.tolist(), columns=['chrom','start', 'ref', 'alt', 'assemb'])

        """ merge afr eur cohorts """
        merged_df = pd.merge(aa_df, ea_df, on=aa_df.columns.to_list(), how = "outer", indicator=True)
        shared_variants = merged_df[merged_df["_merge"] == "both"]
        aa_specific = merged_df[merged_df["_merge"] == "left_only"]
        ea_specific = merged_df[merged_df["_merge"] == "right_only"]

    print("combining shared variants for chroms...")
    shared_variants_all = pd.concat(shared_variants_df, ignore_index=True)
    print("combining AA specific variants for chroms...")
    aa_specific_all = pd.concat(aa_specific_df, ignore_index=True)
    print("combining AA specific variants for chroms...")
    ea_specific_all = pd.concat(ea_specific_df, ignore_index=True)
    return(shared_variants, a_specific, ea_specific)


def pop_shared_variants_LDpruned(self):
    """ Finds shared variants and specfic variants between 
        population cohorts without similar effects/eQTL info
    """
    ## Enter the chromosome no. range that you would like to analyse the data on. By default, it would take the autosomal 
    ## gene model from (chr1:chr22, excluding chrX & chrY), however, set chrom_XY=True for normal gene model analysis.
    chrom_list = []
    for chrom_no in range(chromrange):
        chrom_list.append("chr" + str(chrom_no + 1))
    if chrom_X:
        chrom_list.append("chrX")
        #chrom_list.append("chrY")

    # read chromwise for AFR and EUR:
    aa_filelist = glob(join(aa_filepath, str(tissue) + "*.parquet"))
    aa_dict = {basename(file).split(".")[4] :file for file in aa_filelist}
    ea_filelist = glob(join(ea_filepath, str(tissue) + "*.parquet"))
    ea_dict = {basename(file).split(".")[4] :file for file in ea_filelist}

    shared_variants_df = []
    aa_specific_df = []
    ea_specific_df = []
    for chrom in chrom_list:
        print("processing : {}".format(chrom))

        """ parse afr cohorts"""
        read_aa = pd.read_parquet(aa_dict.get(chrom), engine="fastparquet")
        read_aa.reset_index(inplace=True, drop=True)
        read_aa.drop_duplicates(["variant_id"], inplace=True)
        read_aa["new_col"] = read_aa["variant_id"].map(lambda x:x.split("_"))
        aa_df = pd.DataFrame(read_aa['new_col'].values.tolist(), columns=['chrom','start', 'ref', 'alt', 'assemb'])

        """ parse eur cohorts"""
        read_ea = pd.read_parquet(ea_dict.get(chrom), engine="fastparquet")
        read_ea.reset_index(inplace=True, drop=True)
        read_ea.drop_duplicates(["variant_id"], inplace=True)
        read_ea["new_col"] = read_ea["variant_id"].map(lambda x:x.split("_"))
        ea_df = pd.DataFrame(read_ea['new_col'].values.tolist(), columns=['chrom','start', 'ref', 'alt', 'assemb'])

        """ merge afr eur cohorts """
        merged_df = pd.merge(aa_df, ea_df, on=aa_df.columns.to_list(), how = "outer", indicator=True)
        shared_variants = merged_df[merged_df["_merge"] == "both"]
        aa_specific = merged_df[merged_df["_merge"] == "left_only"]
        ea_specific = merged_df[merged_df["_merge"] == "right_only"]

    print("combining shared variants for chroms...")
    shared_variants_all = pd.concat(shared_variants_df, ignore_index=True)
    print("combining AA specific variants for chroms...")
    aa_specific_all = pd.concat(aa_specific_df, ignore_index=True)
    print("combining AA specific variants for chroms...")
    ea_specific_all = pd.concat(ea_specific_df, ignore_index=True)
    return(shared_variants, a_specific, ea_specific)


test_df = shared_variants[~(shared_variants["phenotype_id_AA"] == shared_variants["phenotype_id_EA"])].iloc[:,[0,1
     ...: ,2,9,10,11]]


    for chrom in chrom_list:
        print("processing : {}".format(chrom))

        """ parse afr cohorts """
        read_aa = pd.read_parquet(aa_dict.get(chrom), engine="fastparquet")
        read_aa.reset_index(inplace=True, drop=True)
        read_aa.drop_duplicates(["variant_id"], inplace=True, keep="first") #drop based on most sig pval

        """ parse eur cohorts"""
        read_ea = pd.read_parquet(ea_dict.get(chrom), engine="fastparquet")
        read_ea.reset_index(inplace=True, drop=True)
        read_ea.drop_duplicates(["variant_id"], inplace=True, keep="first") #drop based on most sig pval

        """ merge afr eur cohorts """
        merged_df1 = pd.merge(read_aa, read_ea, on=["variant_id"], how = "outer", suffixes=["_AA", "_EA"], indicator=True)
        shared_variants = merged_df1[merged_df1["_merge"] == "both"] #filter based on sig pval 0.0001 of at least in one cohort
        shared_variants1 = shared_variants.loc[(shared_variants['pval_nominal_AA'] <= pval) | (shared_variants['pval_nominal_EA'] <= pval)]
        aa_specific = merged_df[merged_df["_merge"] == "left_only"]
        ea_specific = merged_df[merged_df["_merge"] == "right_only"]

        """ select variants with same geneid """
        shared_variants_set1 = shared_variants[shared_variants["phenotype_id_AA"] == shared_variants["phenotype_id_EA"]]
        shared_variants_set1["logpval_nominal_AA"] = -(np.log10(shared_variants_set1["pval_nominal_AA"]))
        shared_variants_set1["logpval_nominal_EA"] = -(np.log10(shared_variants_set1["pval_nominal_EA"]))

        """ filter out variants with disparate geneid"""
        shared_variants_set2 = shared_variants[~(shared_variants["phenotype_id_AA"] == shared_variants["phenotype_id_EA"])]
        #shared_variants[~(shared_variants["phenotype_id_AA"] == shared_variants["phenotype_id_EA"])].iloc[:,[0,1,2,9,10,11]]

        """ calculate variant-gene pairwise pval_correlation for each gene"""
        shared_variants_corr = shared_variants_set1.groupby("phenotype_id_AA")[['logpval_nominal_AA','logpval_nominal_EA']].corr()
        variants_corr = shared_variants_corr.iloc[0::2,-1]
        corr_df = variants_corr.reset_index()





In [1939]: shared_vars_thres.groupby("phenotype_id").size().sort_values()
Out[1939]:
phenotype_id
ENSG00000127472.10      1
ENSG00000215859.8       1
ENSG00000121940.15      1
ENSG00000215915.9       1
ENSG00000118729.11      1
                     ...
ENSG00000143452.15    281
ENSG00000224468.3     352
ENSG00000204138.12    425
ENSG00000162782.15    493
ENSG00000280670.2     522
Length: 546, dtype: int64

In [1940]: shared_vars_thres.groupby("phenotype_id").size().sort_values().idxmax()
Out[1940]: 'ENSG00000280670.2'

In [1941]: idx=shared_vars_thres.groupby("phenotype_id").size().sort_values().idxmax()

In [1942]: corr_df[corr_df["phenotype_id"] == idx]
Out[1942]:
          phenotype_id  beta_correlation
542  ENSG00000280670.2          0.968443

In [1943]: shared_vars_thres[shared_vars_thres["phenotype_id"] == idx]
Out[1943]:
             phenotype_id               variant_id  ...  logpval_nominal_AA  logpval_nominal_EA
422542  ENSG00000280670.2    chr1_45309100_A_G_b38  ...            2.214371           10.143605
422589  ENSG00000280670.2    chr1_45325482_G_A_b38  ...            2.467902           18.581221
422612  ENSG00000280670.2  chr1_45351837_TTA_T_b38  ...            5.285209           37.050771
422769  ENSG00000280670.2    chr1_45421219_A_G_b38  ...            3.072773           34.737830
422776  ENSG00000280670.2    chr1_45424483_T_C_b38  ...            2.199721           15.473680
...                   ...                      ...  ...                 ...                 ...
424490  ENSG00000280670.2    chr1_46140056_A_C_b38  ...            2.725505           26.037682
424491  ENSG00000280670.2   chr1_46140643_T_TA_b38  ...            1.216548           19.420809
424516  ENSG00000280670.2    chr1_46155214_C_A_b38  ...            1.769337           16.489946
424535  ENSG00000280670.2    chr1_46163981_A_C_b38  ...            1.769337           16.052593
424546  ENSG00000280670.2    chr1_46168202_T_C_b38  ...            1.913662           15.502484

[522 rows x 19 columns]

In [1944]: corr_df.iloc[np.where(corr_df["beta_correlation"] >= 0.1)]
Out[1944]:
           phenotype_id  beta_correlation
0    ENSG00000000971.15          0.488487
1    ENSG00000001460.17          0.959270
2    ENSG00000007341.18          0.985578
7    ENSG00000011021.21          0.238936
8    ENSG00000025800.13          0.538663
..                  ...               ...
537   ENSG00000278811.4          0.663445
538   ENSG00000279778.1          0.890519
540   ENSG00000280378.1          0.124307
542   ENSG00000280670.2          0.968443
545   ENSG00000284237.1          0.807309

[301 rows x 2 columns]




def my_plotter(ax, data1, data2, param_dict):
    """
    A helper function to make a graph

    Parameters
    ----------
    ax : Axes
        The axes to draw to

    data1 : array
       The x data

    data2 : array
       The y data

    param_dict : dict
       Dictionary of kwargs to pass to ax.plot

    Returns
    -------
    out : list
        list of artists added
    """
    out = ax.plot(data1, data2, **param_dict)
    return out

fig, ax = plt.subplots(1, 1)
my_plotter(ax, data1, data2, {'marker':'x'})

##################################################

fig, axes = plt.subplots(nrows=nrows)
fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
axes[0].set_title(cmap_category + ' colormaps', fontsize=14)
fig.save_fig()
fig.tight_layout()
##################################################

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 2, 100)

plt.plot(x, x, label='linear')
plt.plot(x, x**2, label='quadratic')
plt.plot(x, x**3, label='cubic')

plt.xlabel('x label')
plt.ylabel('y label')

plt.title("Simple Plot")

plt.legend()

plt.show()



