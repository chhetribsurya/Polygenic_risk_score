
# initial settings
N = 100
I = 1
S = N - I
R = 0
beta = 0.2
gamma = 1./10



                                     
class StochasticSIR:

    def __init__(self,beta,gamma,S,I,R):
        self.S = S
        self.I = I
        self.R = R
        self.beta = beta
        self.gamma = gamma
        self.t = 0.
        self.N = S + I + R
        self.trajectory = np.array([[self.S, self.I, self.R]])
        self.times = None

    def run(self, T=None, make_traj=True):
        """The Gillespie algorithm."""
        if T is None:
            T = sys.maxsize
        self.times = [0.]
        t0 = self.t
        transition = 1
        while self.t < t0 + T:
            transition, t = self.step()
            if not transition:
                return self.t
            if make_traj: self.trajectory = np.concatenate(
                (self.trajectory, [[self.S,self.I,self.R]]), axis=0)
            self.times.append(self.t)
return self.tB


# input files and initial settings:

aa_infile = ""
ea_infile = ""
pattern = ""
ldgenome_infile = ""
ea_pval = ""
aa_pval = ""


class PRSpredictAA:
class PRSpredictEA:
class PRSpredictALL:

def __init__(self, ea_filepath, aa_filepath, tissue, chromrange, chromX=False, ldgenome_infile, pval):
    self.ea_filepath = ea_filepath
    self.aa_filepath = aa_filepath
    self.tissue = tissue
    self.ldgenome_infile = ldgenome_infile
    self.pval = pval
    self.chromrange = chromrange
    self.chromX = chromX
    pass


def intersect_parquet_chromwise(self):

    """ Enter the chromosome no. range that you would like to analyse the data on. By default, it would take the autosomal 
        gene model from (chr1:chr22, excluding chrX & chrY), however, set chrom_XY=True for normal gene model analysis.
    """
    chrom_list = []
    for chrom_no in range(self.chromrange):
        chrom_list.append("chr" + str(chrom_no + 1))
    if chrom_X:
        chrom_list.append("chrX")
        #chrom_list.append("chrY")

    # read chromwise for AFR and EUR:
    ea_filelist = glob(join(self.ea_filepath, str(self.tissue) + "*.parquet"))
    ea_dict = {basename(file).split(".")[4] :file for file in ea_filelist}
    aa_filelist = glob(join(self.aa_filepath, str(self.tissue) + "*.parquet"))
    aa_dict = {basename(file).split(".")[4] :file for file in aa_filelist}

    for chrom in chrom_list:
        print("processing : {}".format(chrom))

        """ parse afr cohorts"""
        read_aa = pd.read_parquet(aa_dict.key(chrom), engine="fastparquet")
        read_aa.reset_index(inplace=True, drop=True)
        read_aa["new_col"] = read_aa["variant_id"].map(lambda x:x.split("_"))
        aa_newdf = pd.DataFrame(read_aa['new_col'].values.tolist(), columns=['chrom','start', 'ref', 'alt', 'assemb'])
        #aa_newdf["end"] = aa_newdf["start"].astype(int) + 1
        aa_df = pd.concat([aa_newdf.loc[:,["chrom", "start", "stop"]], read_aa.loc[:,["tss_distance", "pval_nominal", "slope", "slope_se"]]], axis=1)
        #aa_bed = pybedtools.BedTool.from_dataframe(aa_df)

        """ parse eur cohorts"""
        read_ea = pd.read_parquet(ea_dict.key(chrom), engine="fastparquet")
        read_ea.reset_index(inplace=True, drop=True)
        read_ea["new_col"] = read_ea["variant_id"].map(lambda x:x.split("_"))
        ea_newdf = pd.DataFrame(read_ea['new_col'].values.tolist(), columns=['chrom','start', 'ref', 'alt', 'assemb'])
        #ea_newdf["end"] = ea_newdf["start"].astype(int) + 1
        ea_df = pd.concat([ea_newdf.loc[:,["chrom", "start", "stop"]], read_ea.loc[:,["tss_distance", "pval_nominal", "slope"]]], axis=1)
        #ea_bed = pybedtools.BedTool.from_dataframe(ea_df)

        """ filter shared variants """
        cpg_bed_intersect = ea_bed.intersect(aa_bed wa = True, wb = True) 
    pass


        """ filter that didn't merge """
        #merged_df = pd.merge(aa_df, ea_df, on=aa_df.columns.to_list(), how="inner", suffixes=("_AA", "_EA"))
        #merged_df = pd.merge(merged_df, merged_df1, on=['chrom','start'], how = "outer", indicator=True)
        #merged_df[~(merged_df["alt_AA"] == merged_df["alt_EA"])].iloc[:, [0,1,3,6]]


new_df = read_aa[read_aa.groupby(["variant_id"])["pval_nominal"].idxmin()]

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


import numpy as np, pandas as pd
import pybedtools
import os, re
from glob import glob
from os.path import basename, join, splitext
from collections import defaultdict
import time

start_time = time.time()

## Set output dir:
output_dir = "/gpfs/gpfs1/home/schhetri/DNAme_paper"

## Can be ChromHMM file or IDEAS segmentation file:
ideas_file_path = "/gpfs/gpfs1/home/schhetri/DNAme_paper/ideas_seg_files"

## Can be ChromHMM file or IDEAS segmentation file:
cpg_file_path = "/gpfs/gpfs1/home/schhetri/DNAme_paper/CpG_files"

## IDEAS mnemonics file:
ideas_mnemonics_file = os.path.expanduser("~/for_encode/spp/analysis_scripts/ideas_table.txt")

## Cell-Tissue type to analyse:
samples = ["GM12878", "H1hESC", "HepG2", "K562"]

## Generate dict for ideas files:
ideas_file_list = glob(join(ideas_file_path, "*.bed.srt"))
ideas_filedict = {basename(file).split("_ideas")[0] :file for file in ideas_file_list}

## Generate dict for cpg files:
cpg_file_list = glob(join(cpg_file_path, "*.cov.srt"))
cpg_filedict = {basename(file).split("_CpG")[0] :file for file in cpg_file_list}

file_index_dict = defaultdict(list)
for key,value in cpg_filedict.iteritems():
    key_term =  key.split("_")[0]
    file_index_dict[key_term].append(key)

## Generate output dirs:
dname_ideas_output = join(output_dir, "DNAme_ideas_output") 
plot_output = join(dname_ideas_output, "plots")

if not os.path.exists(dname_ideas_output):
    os.makedirs(dname_ideas_output)

if not os.path.exists(plot_output):
    os.makedirs(plot_output)

## Select the state number for analysis:
## range() is exclusive:
select_state_num = range(1,10) + [15,16] + range(17,20+1) # regulatory regions
select_state_num = range(1,37) # all states

## Read the mnemonics file, create a dict of states with state ID
read_file = pd.read_csv(ideas_mnemonics_file, sep="\t")
state_list = read_file["Mnemonics"]
state_num_list = [ i for i in range(1,len(state_list)+1) ]
ideas_state_dict = dict(zip(state_num_list,state_list))
Target_state = [ideas_state_dict[each] for each in select_state_num] # Depends on select_state_num var
vals_to_replace_dict = {value:str(key)+"_"+value for key,value in ideas_state_dict.iteritems()}


def final_peaks_model(peak_input_file):
    input_file = peak_input_file
    peak_df = pd.read_csv(input_file, sep="\t", header=None) # skiprows=[0]
    peak_select_df = peak_df.iloc[:, [0,1,2,3]]
    peak_select_df.columns = ["chrom", "start", "end", "state"]
    print "\nCurrent dimension of the peak model:\n", peak_select_df.shape
    peak_select_df = peak_select_df.drop_duplicates()
    print "Dropping duplicates if any, current dimension of the peak model:\n", peak_select_df.shape
    return(peak_select_df)


def load_cpg_pybedtool_object(file_name_with_full_path):
    print " \nProcessing Cpg bed file\n "
    cpg_df =  pd.read_csv(file_name_with_full_path, sep="\t", header=None)
    cpg_df = cpg_df.iloc[:,0:6]
    Cpg_bed_file = pybedtools.BedTool.from_dataframe(cpg_df)
    return(Cpg_bed_file)


def generate_peaks_binned_coords(peaks_coordinates_info, upstream_range, downstream_range, bin_size):
    peaks_df =  peaks_coordinates_info.sort_values(["chrom", "start", "end"])
    upstream = upstream_range
    downstream = downstream_range
    bin_size = bin_size
    nrows =  peaks_df.shape[0]

    bins = range(-upstream, (downstream), bin_size)
    bin_len = len(bins)
    peaks_concat_df = pd.concat([peaks_df]*bin_len, ignore_index="TRUE")
    peaks_sorted_df = peaks_concat_df.sort_values(["chrom","start","end"])

    ### Copy the bin list that is deep copy:
    bin_start_list = bins[:]
    bin_end_list = []
    for each in bin_start_list:
        bin_end_list.append(each+bin_size)

    bin_df = pd.DataFrame({"bin_start":bin_start_list, "bin_end":bin_end_list})
    bin_df = bin_df.iloc[:,[1,0]] # Arrange bin_start as first col and bin_start as 2nd col.
    bin_concat_df = pd.concat([bin_df]*nrows, ignore_index="TRUE")

    ### Combine the peaks df and bin df by cbind or column wise:
    temp_peaks_df = pd.concat([peaks_sorted_df.reset_index(), bin_concat_df], axis = 1)
    temp_peaks_df["peaks_midpoint"] = (temp_peaks_df["start"] + temp_peaks_df["end"])/2
    temp_peaks_df["peaks_midpoint"] = temp_peaks_df["peaks_midpoint"].round().astype(int)
    final_peaks_df = temp_peaks_df.loc[:,["chrom", "peaks_midpoint", "bin_start", "bin_end", "start", "end", "state"]]

    """ 
    chrom_start = tss_midpt + (bin_start); 
    chrom_end = tss_midpt + (bin_end) """
    final_peaks_df["chrom_start"] = final_peaks_df["peaks_midpoint"] + final_peaks_df["bin_start"]
    final_peaks_df["chrom_end"] = final_peaks_df["peaks_midpoint"] + final_peaks_df["bin_end"]

    select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'bin_start', u'bin_end', u'peaks_midpoint', u'state']
    final_peaks_df = final_peaks_df.loc[:,select_cols]

    ### Resolve bins greater than actual loci position(eg: midpoint = 1000,leading to -ve chrom_start(crossing 0) coords 
    ### for 2kb upstream and 2kb downstream, particularly for 2kb downstream which exceeds the 0 position)
    final_peaks_df = final_peaks_df[~final_peaks_df["chrom_start"] < 0] 
    return(final_peaks_df)


def generate_peaks_binned_perc_meth(peaks_final_bedfile, meth_file_list, **kwargs):
    # file_name =  kwargs["files_basename"]
    print "kwargs: ", kwargs
    peaks_bedfile = peaks_final_bedfile

    # master_dict = {}
    for idx, each_file in enumerate(meth_file_list):
        cpg_bed_file = load_cpg_pybedtool_object(each_file)
        cpg_bed_intersect = cpg_bed_file.intersect(peaks_bedfile, wa = True, wb = True)     
        print(cpg_bed_intersect.head())  

        ### Working with the dataframes; reading the output file of pybedtool intersect:
        df_file = pd.read_csv(cpg_bed_intersect.fn, sep = "\t", header = None)
        df_ordered = df_file.iloc[:, 0:12]
        df_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "peak_chrom", "peak_start", "peak_end", "bin_start", "bin_end", "peaks_midpoint", "state"]
        df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
        df_grouped =  df_ordered.groupby(["bin_start", "bin_end", "state"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
        print "Dimension of currently intersected peak file is", df_grouped.reset_index().shape
        grouped_df = df_grouped.reset_index(name="meth_percent")
    return(grouped_df)


def generate_peaks_perc_meth(peaks_final_bedfile, meth_file_list, **kwargs):
    print "kwargs: ", kwargs
    peaks_bedfile = peaks_final_bedfile

    # master_dict = {}
    for idx, each_file in enumerate(meth_file_list):
        cpg_bed_file = load_cpg_pybedtool_object(each_file)
        cpg_bed_intersect = cpg_bed_file.intersect(peaks_bedfile, wa = True, wb = True)     
        print(cpg_bed_intersect.head())  

        ### Working with the dataframes; reading the output file of pybedtool intersect:
        df_file = pd.read_csv(cpg_bed_intersect.fn, sep = "\t", header = None)
        df_ordered = df_file.iloc[:, 0:9]
        df_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "peak_chrom", "peak_start", "peak_end", "state"]
        df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
        df_grouped =  df_ordered.groupby(["state"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
        print "Dimension of currently intersected peak file is", df_grouped.reset_index().shape
        grouped_df = df_grouped.reset_index(name="meth_percent")
    return(grouped_df)


## For binned methylation profile
def main():
    concat_list = []
    for key, value in file_index_dict.iteritems():
        print("Processing {}\n".format(key))
        peaks_coord_df = final_peaks_model(ideas_filedict.get(key))
        peaks_coord_df = peaks_coord_df[peaks_coord_df["state"].isin(Target_state)]
        peaks_coord_df = peaks_coord_df.replace({"state" : vals_to_replace_dict})
        peaks_midpoint_coord_df = generate_peaks_binned_coords(peaks_coord_df, 2000, 2000, 100)
        peaks_bed_file = pybedtools.BedTool.from_dataframe(peaks_midpoint_coord_df)
        sub_key = file_index_dict.get(key)
        for each_cpg_rep in sub_key:
            print("Processing {}\n".format(each_cpg_rep))
            list_of_cpg_files = [cpg_filedict.get(each_cpg_rep)]
            plot_data_set = generate_peaks_binned_perc_meth(peaks_bed_file, list_of_cpg_files)
            plot_data_set["sample"] = each_cpg_rep
            concat_list.append(plot_data_set)

    final_meth_df = pd.concat(concat_list)
    final_meth_df.to_csv(join(dname_ideas_output, "binned_ideas_methpercent.txt"), sep ="\t", header = True, index = True)
    final_meth_df["state"] =  final_meth_df["state"].str.replace(r"\d+_", "")
    final_methtab_df = final_meth_df.pivot(index="sample",columns="state",values="meth_percent").reset_index()
    final_methtab_df.to_csv(join(dname_ideas_output, "binned_ideas_methpercent_table.txt"), sep ="\t", header = True, index = False)


if __name__ == '__main__':
    main()


## For unbinned whole peak meth profile:
def main():
    concat_list = []
    for key, value in file_index_dict.iteritems():
        print("Processing {}\n".format(key))
        peaks_coord_df = final_peaks_model(ideas_filedict.get(key))
        peaks_coord_df = peaks_coord_df[peaks_coord_df["state"].isin(Target_state)]
        peaks_coord_df = peaks_coord_df.replace({"state" : vals_to_replace_dict})
        # peaks_midpoint_coord_df = generate_peaks_binned_coords(peaks_coord_df, 200, 200, 100)
        peaks_bed_file = pybedtools.BedTool.from_dataframe(peaks_coord_df)
        sub_key = file_index_dict.get(key)
        for each_cpg_rep in sub_key:
            print("Processing {}\n".format(each_cpg_rep))
            list_of_cpg_files = [cpg_filedict.get(each_cpg_rep)]
            plot_data_set = generate_peaks_perc_meth(peaks_bed_file, list_of_cpg_files)
            plot_data_set["sample"] = each_cpg_rep
            concat_list.append(plot_data_set)

    final_meth_df = pd.concat(concat_list)
    final_meth_df.to_csv(join(dname_ideas_output, "Ideas_methpercent_states.txt"), sep ="\t", header = True, index = True)
    final_meth_df["state"] =  final_meth_df["state"].str.replace(r"\d+_", "")
    final_methtab_df = final_meth_df.pivot(index="sample",columns="state",values="meth_percent").reset_index()
    final_methtab_df.to_csv(join(dname_ideas_output, "Ideas_methpercent_table_states.txt"), sep ="\t", header = True, index = False)


if __name__ == '__main__':
    main()

print("Time for analysis : {}".format(time.time()-start_time))
print("Task completed")


######################################################
######################################################



import re
def sorted_values(l):
    """ Sorts the given iterable in the way that is expected.
 
    Required arguments:
    l -- The iterable to be sorted.
 
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)

df_result = result_df[2]



# fig, ax = plt.subplots(1,1)
df_result = result_df[2]
df_result = df_result.groupby(["chrom"]).median().reset_index()
sort_list = df_result["chrom"]
df_new = pd.Series(sorted_values(sort_list), name="chrom")
df_final = pd.merge(df_result, df_new, on=["chrom"])

# set fig and axes:
fig, ax = plt.subplots()

# set width of bar
barWidth = 0.25
 
# set height of bar
bars1 = df_final[""]
bars2 = [28, 6, 16, 5, 10]
bars3 = [29, 3, 24, 25, 17]
 
# Set position of bar on X axis
r1 = np.arange(len(bars1))
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]


barWidth = 0.9
bar
base_x1 =  
ax.bar(df_final["chrom"], df_final["beta_correlation"], width = barWidth, color="orange")
ax.bar(df_final["chrom"], df_final["rsquare"], width = barWidth, color="red")
ax.xlabel("Value")
# ax.invert_yaxis() #plt.gca().invert_xaxis() # [::-1] reverses the list
fig.autofmt_xdate(rotation=45) #ax.tick_params(axis='x', labelrotation=45)


######################################################
######################################################
# set some xlim and ylim in Seaborn lmplot facetgrid
#You need to get hold of the axes themselves. Probably the cleanest way is to change your last row:
#https://stackoverflow.com/questions/25212986/how-to-set-some-xlim-and-ylim-in-seaborn-lmplot-facetgrid?noredirect=1&lq=1

lm = sns.lmplot('X','Y',df,col='Z',sharex=False,sharey=False)
#Then you can get hold of the axes objects (an array of axes):

axes = lm.axes
# After that you can tweak the axes properties

axes[0,0].set_ylim(0,)
axes[0,1].set_ylim(0,)

######################################################
######################################################

# Faceted logistic regression : seaborn plot
#https://seaborn.pydata.org/examples/logistic_regression.html

import seaborn as sns
sns.set(style="darkgrid")

# Load the example titanic dataset
df = sns.load_dataset("titanic")

# Make a custom palette with gendered colors
pal = dict(male="#6495ED", female="#F08080")

# Show the survival proability as a function of age and sex
g = sns.lmplot(x="age", y="survived", col="sex", hue="sex", data=df,
               palette=pal, y_jitter=.02, logistic=True)
g.set(xlim=(0, 80), ylim=(-.05, 1.05))

######################################################
######################################################

# Targeting a non-default axes with seaborn plots
# http://alanpryorjr.com/visualizations/seaborn/barplot/barplot/
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2)
sns.countplot(data=df,
                  y = 'Category',
                  hue = 'islong',
                  saturation=1,
                  ax=ax[1])

ax[1].tick_params(axis='x', labelrotation=45)
#fig.autofmt_xdate(rotation=45)

######################################################
######################################################

import numpy as np

np.random.seed(19680801)
data = np.random.randn(2, 100)
cm = plt.cm.get_cmap('RdYlBu')

value = np.linspace(0,11)
fig, axs = plt.subplots(2, 2, figsize=(5, 5))
axs[0, 0].hist(data[0])
im = axs[1, 0].scatter(data[0], data[1], c=data[0], cmap=cm, s=4, alpha=0.8)

axs[0, 1].plot(data[0], data[1])
axs[1, 1].hist2d(data[0], data[1])
sns.despine(offset=5, ax=axs[0, 1])
sns.despine(offset=5, ax=axs[1, 0])
axs[0,1].set_xlim(-2,2)
axs[0,1].set_ylabel("right")
axs[1, 0].set_xlabel("scatter test")
axs[1, 0].set_xlim(-2,2)

#plot colorbar to particular axes:
cbar = fig.colorbar(im, ax=axs[1, 0], 
                    orientation = "vertical", shrink=1)

#set the colorbar ticks and tick labels
# cbar.set_ticks(np.arange(0, 1.1, 0.5))

#For setting categorical values to colorbar:
max_num=max(data[0])
min_num=min(data[0])
medium = (max_num + min_num)/2

cbar.set_ticks([min_num, medium, max_num])
cbar.set_ticklabels(['low', 'medium', 'high'])

#####################################################
#####################################################






