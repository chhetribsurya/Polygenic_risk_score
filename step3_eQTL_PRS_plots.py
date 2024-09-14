#!/usr/bin/env python
import numpy as np, pandas as pd
import seaborn as sns
import pybedtools
import os, re, pickle

from glob import glob
from os.path import basename, join, splitext
from collections import defaultdict
from functools import reduce 

import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy import stats

"""supress warnings"""
import time, warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None

"""Input parameters"""
tissue = "Cells_EBV-transformed_lymphocytes"
output_path = "/Users/suryachhetri/datasets/prs_project/final_hg19/EBV_LCL" 

# chromrange = 22
# chrom_X = False
# pval = 10e-5
# zscorethresh = 0.1

# compressed formats: read/write with pickle
#cis_eqtl_unfilt_df = pop_shared_variants_betabased(aa_filepath, ea_filepath, chromrange, pval, tissue, chrom_X=False)
#cis_eqtl_unfilt_df.to_pickle(join(output_path, tissue + "_shared_allcis_eQTL_unfilt_stats.pkl"))
cis_eqtl_unfilt_df = pd.read_pickle(join(output_path, tissue + "_shared_allcis_eQTL_unfilt_stats.pkl"))


""" boxplot for genomewide beta corr and coeff variation """
df_result = cis_eqtl_unfilt_df
df_new = df_result.loc[:, ["chrom", "beta_correlation", "Zsquare_rsquare"]]
df_melt = pd.melt(df_new, "chrom", var_name="measurement")

""" using seaborn backend and placed in subplot"""
sns.set(style="ticks")
fig, ax = plt.subplots()

# Draw a nested boxplot:
# fig.set_size_inches(8, 6)
axs = sns.boxplot(x="chrom", y="value",
            hue="measurement", palette=["r", "b"],
            data=df_melt, \
            ax=ax, \
            #"""remove outlier"""
            flierprops = dict(markerfacecolor = '0.3', markersize = 1))

sns.despine(offset=5, trim=True)

# set axes and labels for sns subplot:
axs.set_xlabel("Chromosome", fontsize=10)
axs.set_ylabel("Correlation | Coeff of Variation", fontsize=10)
axs.set_title("Genome Wide Distribution (r|rsquare)", weight='bold', fontsize=12)
axs.tick_params(axis='x', labelrotation=45, labelcolor='black', labelsize=9)
axs.tick_params(axis='y', labelrotation=0, labelcolor='black', labelsize=9)

# Put the legend out of the figure
ax.legend(loc=4, prop={'size': 10}, edgecolor="black")
ax.axhline(0.8, color='tab:blue', linestyle='dashed', linewidth=2)

# Save figure:
fig.tight_layout()
fig.savefig(join(plot_path, '{}_genomewide_betaStats_coeff_boxplot.pdf'.format(tissue)))
plt.close()


"""PRS bar plots"""
input_dir = "/Users/suryachhetri/datasets/prs_project/final_hg19/EBV_LCL/prs_output"


"""eQTL-informed - pearson approach"""
id1 = "0.1_AA"; id2 = "AA"
aa_filt0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
aa_filt1 = pd.Series([282, 482, 1444, 3252, 9050, 31948, 82822, 113465], name="snp_count")
aa_filt = pd.concat([aa_filt0, pd.DataFrame(aa_filt1)], axis=1)
aa_filt["cohort"] = "AFR"

id1 = "0.1_EUR"; id2 = "EUR"
eur_filt0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
eur_filt1 = pd.Series([283, 483, 1450, 3269, 9114, 32107, 83218, 113957], name="snp_count")
eur_filt = pd.concat([eur_filt0, pd.DataFrame(eur_filt1)], axis=1)
eur_filt["cohort"] = "EUR"

"""eQTL-informed - bayesian approach""" #coloc-genes)
# id1 = "colocfilt0.5_AA"; id2 = "AA"
# aa_weighted0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
# aa_weighted1 = pd.Series([375, 657, 2092, 4925, 14521, 54374, 144821, 199633], name="snp_count")
# aa_weighted = pd.concat([aa_weighted0, pd.DataFrame(aa_weighted1)], axis=1)
# aa_weighted["cohort"] = "AFR"

# id1 = "colocfilt0.5_EUR"; id2 = "EUR"
# eur_weighted0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
# eur_weighted1 = pd.Series([376, 657, 2096, 4944, 14604, 54644, 145548, 200577], name="snp_count")
# eur_weighted = pd.concat([eur_weighted0, pd.DataFrame(eur_weighted1)], axis=1)
# eur_weighted["cohort"] = "EUR"

"""eQTL-informed - bayesian approach""" #coloc-genes)
id1 = "colocfilt0.1_AA"; id2 = "AA"
aa_weighted0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
aa_weighted1 = pd.Series([257, 458, 1468, 3405, 9814, 36272, 95674, 131639], name="snp_count")
aa_weighted = pd.concat([aa_weighted0, pd.DataFrame(aa_weighted1)], axis=1)
aa_weighted["cohort"] = "AFR"

id1 = "colocfilt0.1_EUR"; id2 = "EUR"
eur_weighted0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
eur_weighted1 = pd.Series([258, 458, 1471, 3412, 9871, 36398, 96022, 132110], name="snp_count")
eur_weighted = pd.concat([eur_weighted0, pd.DataFrame(eur_weighted1)], axis=1)
eur_weighted["cohort"] = "EUR"

"""cis-genes"""
id1 = "unfilt_AA"; id2 = "AA"
aa_unfilt0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
aa_unfilt1 = pd.Series([414, 724, 2389, 5666, 16884, 63837, 170515, 235009], name="snp_count")
aa_unfilt = pd.concat([aa_unfilt0, pd.DataFrame(aa_unfilt1)], axis=1)
aa_unfilt["cohort"] = "AFR"

id1 = "unfilt_EUR"; id2 = "EUR"
eur_unfilt0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
eur_unfilt1 = pd.Series([415, 725, 2396, 5691, 16971, 64146, 171316, 236020], name="snp_count")
eur_unfilt = pd.concat([eur_unfilt0, pd.DataFrame(eur_unfilt1)], axis=1)
eur_unfilt["cohort"] = "EUR"

"""gwas-baseline"""
id1 = "gwasbaseline_AA"; id2 = "AA"
aa_gwasbaseline0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
aa_gwasbaseline1 = pd.Series([484, 873, 2981, 7411, 23481, 93599, 259178, 363755], name="snp_count")
aa_gwasbaseline = pd.concat([aa_gwasbaseline0, pd.DataFrame(aa_gwasbaseline1)], axis=1)
aa_gwasbaseline["cohort"] = "AFR"

id1 = "gwasbaseline_EUR"; id2 = "EUR"
eur_gwasbaseline0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
eur_gwasbaseline1 = pd.Series([500, 913, 3305, 8717, 30363, 133546, 394937, 568503], name="snp_count")
eur_gwasbaseline = pd.concat([eur_gwasbaseline0, pd.DataFrame(eur_gwasbaseline1)], axis=1)
eur_gwasbaseline["cohort"] = "EUR"

# concat dataframes
frames = [aa_filt, eur_filt, aa_weighted, eur_weighted, aa_unfilt, eur_unfilt, aa_gwasbaseline, eur_gwasbaseline]
result = pd.concat(frames, keys=['aa_filt', 'eur_filt', 'aa_weighted', 'eur_weighted', 'aa_unfilt', 'eur_unfilt',  'aa_gwasbaseline', 'eur_gwasbaseline'])

final_df = result.reset_index().drop(columns="level_1", axis=0)
final_df = final_df.rename(columns={"level_0": "Category"})
final_df = final_df.loc[~final_df["Category"].str.contains("weighted")]


"""BMI"""
#bmi_df1 = bmi_df[~(bmi_df["Cohorts"] == "ALL")]

g = sns.barplot(x = 'pval_names', y = 'r2', hue = 'Category', data = final_df,
            palette = 'hls',
            #col="Cohorts",
            #order = ['aa_filt', 'eur_filt', 'aa_unfilt', 'eur_unfilt', 'aa_weighted', 'eur_weighted', 'aa_gwasbaseline', 'eur_gwasbaseline'],  
            capsize = 0.05,             
            saturation = 8,             
            errcolor = 'gray', errwidth = 2,  
            ci = 'sd'   
            )

box = g.get_position()
g.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # resize position

# Put a legend to the right side
g.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)

plot_path = input_dir
plt.savefig(join(plot_path, '1_plot_test.pdf'))
plt.close()


#########################

# LD predfunct PRS plots

#########################
import numpy as np, pandas as pd
import seaborn as sns
import pybedtools
import os, re, pickle

from glob import glob
from os.path import basename, join, splitext
from collections import defaultdict
from functools import reduce 

import matplotlib; 
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy import stats

# #final_select_df = final_df.loc[~final_df["Category"].str.contains("weighted")]
# final_select_df = final_df.copy()
# final_select_df["Group"] = ""
# final_select_df.loc[final_select_df["Category"].str.contains("_filt"), "Group"] = "eQTL-informed-genes(pearson)"
# final_select_df.loc[final_select_df["Category"].str.contains("_weighted"), "Group"] = "eQTL-informed-genes(bayesian)"
# final_select_df.loc[final_select_df["Category"].str.contains("_unfilt"), "Group"] = "cis-genes"
# final_select_df.loc[final_select_df["Category"].str.contains("gwasbaseline"), "Group"] = "gwas-baseline"
# final_select_df["R2perSNP"] = final_select_df["r2"]/final_df["snp_count"].astype(float)

input_dir = "/Users/suryachhetri/datasets/prs_project/final_hg19/EBV_LCL/prs_output/ldpredfunct_prs"
prefix = "WITHbaseOnlyPPH4_Final_ldpredfunct_PRS_computedR2_plus_persnpR2.txt"

infile = join(input_dir, prefix)
prs_df = pd.read_csv(infile, sep="\t")
finalprs_df = prs_df[prs_df["pvalcategory"] != 5e-8]


##########
#Unadjusted R2persnp --direct from LDpredfunct
melt_df1 = pd.melt(finalprs_df, id_vars = ['pvalcategory', "cohort"], 
                  value_vars= ['RawPRS_R2perSNP', 'WeightedPRS_R2perSNP', ])

aa_melt_df1  = melt_df1[melt_df1["cohort"] == "AA"]
eur_melt_df1  = melt_df1[melt_df1["cohort"] == "EUR"]


#################
#Adjusted R2perSNP by covariates 

melt_df2 = pd.melt(finalprs_df, id_vars = ['pvalcategory', "cohort"], 
                  value_vars= ['RawPRS_R2perSNP_adj', 'WeightedPRS_R2perSNP_adj', ])

aa_melt_df2  = melt_df2[melt_df2["cohort"] == "AA"]
eur_melt_df2  = melt_df2[melt_df2["cohort"] == "EUR"]

##########
#Unadjusted R2 --direct from LDpredfunct

melt_df3 = pd.melt(finalprs_df, id_vars = ['pvalcategory', "cohort"], 
                  value_vars= ['raw_prs', 'weighted_prs'])

aa_melt_df3  = melt_df3[melt_df3["cohort"] == "AA"]
eur_melt_df3  = melt_df3[melt_df3["cohort"] == "EUR"]


#########
#Adjusted R2 by covariates

melt_df4 = pd.melt(finalprs_df, id_vars = ['pvalcategory', "cohort"], 
                  value_vars= ['raw_prs_adj', 'weighted_prs_adj'])

aa_melt_df4  = melt_df4[melt_df4["cohort"] == "AA"]
eur_melt_df4  = melt_df4[melt_df4["cohort"] == "EUR"]

########
#Combined Adjusted and Unadjusted R2
melt_df5 = pd.melt(finalprs_df, id_vars = ['pvalcategory', "cohort"], 
                  value_vars= ['raw_prs', 'weighted_prs','raw_prs_adj', 'weighted_prs_adj'])

aa_melt_df5  = melt_df5[melt_df5["cohort"] == "AA"]
eur_melt_df5  = melt_df5[melt_df5["cohort"] == "EUR"]

##########
#color theme
color = [ sns.color_palette("Paired")[3] ] 
color.append((sns.color_palette("Paired")[7]))

new_color = sns.color_palette("Paired")[2:4] #index 4 exlcuded i.e. only 2 and 3 included
new_color.append((sns.color_palette("Paired")[6]))
new_color.append((sns.color_palette("Paired")[7]))

##########


##############
#concatenated categorical plot
sns.catplot(x = 'pvalcategory', y = 'predfunct_prs', hue = 'cohort', 
            data = finalprs_df,
            palette = color,
            col="cohort",
            kind="bar",
            height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
plt.legend(frameon=False)
plot_path = input_dir
plt.savefig(join(plot_path, '0.1_Plot_CombinedplotCohorts.pdf'))
plt.close()


#################
#Unadjusted R2 perSNP -- direct from LDpredfunct
g = sns.barplot(x = 'pvalcategory', y = 'value', hue = 'variable', 
            data = melt_df2, #melt_df1
            palette = color,
            ci=None
            #col= "cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8         
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )

plt.xlabel("GWAS P-value")
plt.ylabel("R2perSNP (SBP Trait)")

plt.legend(frameon=False)
plot_path = input_dir
plt.savefig(join(plot_path, 'Uanadj_1.1_Plot.pdf'))
plt.close()

g = sns.barplot(x = 'pvalcategory', y = 'value', hue = 'variable', 
            data = aa_melt_df2, #aa_melt_df1
            palette = color,
            ci=None
            #col= "cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8         
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   

            )

plt.xlabel("GWAS P-value")
plt.ylabel("R2perSNP (SBP Trait)")

plt.legend(frameon=False)
plot_path = input_dir
plt.savefig(join(plot_path, 'Unadj_aa_1.2_Plot.pdf'))
plt.close()

g = sns.barplot(x = 'pvalcategory', y = 'value', hue = 'variable', 
            data = eur_melt_df2, #eur_melt_df1
            palette = color,
            ci=None
            #col= "cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8         
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   

            )

plt.xlabel("GWAS P-value")
plt.ylabel("R2perSNP (SBP Trait)")

plt.legend(frameon=False)
plot_path = input_dir
plt.savefig(join(plot_path, 'Unadj_eur_1.3_Plot.pdf'))
plt.close()


#################
#Unadjusted R2 -- direct from LDpredfunct
g = sns.barplot(x = 'pvalcategory', y = 'value', hue = 'variable', 
            data = melt_df3,
            palette = color,
            ci=None
            #col= "cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8         
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )

plt.xlabel("GWAS P-value")
plt.ylabel("R2 (SBP Trait)")

plt.legend(frameon=False)
plot_path = input_dir
plt.savefig(join(plot_path, 'Uanadj_1.1_Plot.pdf'))
plt.close()

g = sns.barplot(x = 'pvalcategory', y = 'value', hue = 'variable', 
            data = aa_melt_df3,
            palette = color,
            ci=None
            #col= "cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8         
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   

            )

plt.legend(frameon=False)
plot_path = input_dir
plt.savefig(join(plot_path, 'Unadj_aa_1.2_Plot.pdf'))
plt.close()

g = sns.barplot(x = 'pvalcategory', y = 'value', hue = 'variable', 
            data = eur_melt_df3,
            palette = color,
            ci=None
            #col= "cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8         
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   

            )

plt.legend(frameon=False)
plot_path = input_dir
plt.savefig(join(plot_path, 'Unadj_eur_1.3_Plot.pdf'))
plt.close()


#################
#Adjusted R2 by covariates
g = sns.barplot(x = 'pvalcategory', y = 'value', hue = 'variable', 
            data = melt_df4,
            palette = color,
            ci=None
            #col= "cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8         
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )

plt.legend(frameon=False)
plot_path = input_dir
plt.savefig(join(plot_path, 'Adj_2.1_Plot.pdf'))
plt.close()

#aa_melt_df4  = melt_df4[melt_df4["cohort"] == "AA"]
g = sns.barplot(x = 'pvalcategory', y = 'value', hue = 'variable', 
            data = aa_melt_df4,
            palette = color,
            ci=None
            #col= "cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8         
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )

plt.legend(frameon=False)
plot_path = input_dir
plt.savefig(join(plot_path, 'Adj_aa_2.2_Plot.pdf'))
plt.close()

#eur_melt_df4  = melt_df4[melt_df4["cohort"] == "EUR"]
g = sns.barplot(x = 'pvalcategory', y = 'value', hue = 'variable', 
            data = eur_melt_df4,
            palette = color,
            ci=None
            #col= "cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8         
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )

plt.legend(frameon=False)
plot_path = input_dir
plt.savefig(join(plot_path, 'Adj_eur_2.3_Plot.pdf'))
plt.close()

#################
#Combined Adjusted and Unadjusted R2
#Adjusted by covariates

g = sns.barplot(x = 'pvalcategory', y = 'value', hue = 'variable', 
            data = melt_df5,
            palette = sns.color_palette("Paired"),
            ci=None
            #col= "cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8         
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )

plt.xlabel("GWAS P-value")
plt.ylabel("R2 (SBP Trait)")

plt.legend(frameon=False)
plot_path = input_dir
plt.savefig(join(plot_path, 'Aggregate_3.1_Plot.pdf'))
plt.close()

#aa_melt_df5  = melt_df5[melt_df5["cohort"] == "AA"]
g = sns.barplot(x = 'pvalcategory', y = 'value', hue = 'variable', 
            data = aa_melt_df5,
            palette = sns.color_palette("Paired"),
            ci=None
            #col= "cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8         
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )

plt.legend(frameon=False)
plot_path = input_dir
plt.savefig(join(plot_path, 'Aggregate_eur_3.2_Plot.pdf'))
plt.close()

#eur_melt_df5  = melt_df5[melt_df5["cohort"] == "EUR"]
g = sns.barplot(x = 'pvalcategory', y = 'value', hue = 'variable', 
            data = eur_melt_df5,
            palette = new_color,
            ci=None
            #col= "cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8         
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )


plt.xlabel("GWAS P-value")
plt.ylabel("R2 (SBP Trait)")

plt.legend(frameon=False)
plot_path = input_dir
plt.savefig(join(plot_path, 'Aggregate_eur_3.3_Plot.pdf'))
plt.close()





#################################################################
g = sns.barplot(x = 'pvalcategory', y = 'value', hue = 'variable', 
            data = melt_df5,
            palette = sns.color_palette("Paired"),
            ci=None
            #col= "cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8         
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )



plt.legend(title='Category', loc="upper right")
plt.legend(frameon=False)

plot_path = input_dir
plt.savefig(join(plot_path, '2.1_Plot.pdf'))
plt.close()


plot_df = melt_df[melt_df["cohort"] == "EUR"]
sns.barplot(x = 'pvalcategory', y = 'value', hue = 'variable', 
            data = plot_df,
            palette = color,
            #col="cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )

#plt.legend()
plot_path = input_dir
plt.savefig(join(plot_path, '2.1_Plot.pdf'))
plt.close()





sns.catplot(x = 'pval_names', y = 'R2perSNP', hue = 'Group', data = final_select_df,
            palette = 'hls',
            col="cohort",
            kind="bar",
            height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
#plt.legend()
plot_path = input_dir
plt.savefig(join(plot_path, '3_CombinedplotCohorts_R2perSNP.pdf'))
plt.close()


#final_select_df1 = final_select_df[~(final_select_df["cohort"] == "AA")]
sns.catplot(x = 'pval_names', y = 'snp_count', hue = 'Group', data = final_select_df,
            palette = 'hls',
            col="cohort",
            kind="bar",
            height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
#plt.legend()
plt.savefig(join(plot_path, '4_CombinedplotCohorts_EURsnpcount.pdf'))
plt.close()


##########################################
##########################################
## Related to plot1 (perSNP heritability):
"""BMI"""
g = sns.barplot(x = 'pval_names', y = 'R2perSNP', hue = 'Category', data = final_select_df,
            palette = 'hls',
            #col="Cohorts",
            #order = ['aa_filt', 'eur_filt', 'aa_unfilt', 'eur_unfilt', 'aa_weighted', 'eur_weighted', 'aa_gwasbaseline', 'eur_gwasbaseline'],  
            capsize = 0.05,             
            saturation = 8,             
            errcolor = 'gray', errwidth = 2,  
            ci = 'sd'   
            )

box = g.get_position()
g.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # resize position

# Put a legend to the right side
g.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)

plot_path=input_dir
plt.savefig(join(plot_path, '1.1_plot_test.pdf'))
plt.close()







###############################################
###############################################

# LDPRED-FUNCT PRS plots:

###############################################
###############################################


#final_select_df = final_df.loc[~final_df["Category"].str.contains("weighted")]
final_select_df = final_df.copy()
final_select_df["Group"] = ""
final_select_df.loc[final_select_df["Category"].str.contains("_filt"), "Group"] = "eQTL-informed-genes(pearson)"
final_select_df.loc[final_select_df["Category"].str.contains("_weighted"), "Group"] = "eQTL-informed-genes(bayesian)"
final_select_df.loc[final_select_df["Category"].str.contains("_unfilt"), "Group"] = "cis-genes"
final_select_df.loc[final_select_df["Category"].str.contains("gwasbaseline"), "Group"] = "gwas-baseline"
final_select_df["R2perSNP"] = final_select_df["r2"]/final_df["snp_count"].astype(float)

#final_df1 = final_df[final_df["Category"] == "eQTL Informed"]
sns.catplot(x = 'pval_names', y = 'r2', hue = 'Group', data = final_select_df,
            palette = 'hls',
            col="cohort",
            kind="bar",
            height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
#plt.legend()
plot_path = input_dir
plt.savefig(join(plot_path, '2.1_Plot_CombinedplotCohorts.pdf'))
plt.close()


sns.catplot(x = 'pval_names', y = 'R2perSNP', hue = 'Group', data = final_select_df,
            palette = 'hls',
            col="cohort",
            kind="bar",
            height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
#plt.legend()
plot_path = input_dir
plt.savefig(join(plot_path, '3_CombinedplotCohorts_R2perSNP.pdf'))
plt.close()


#final_select_df1 = final_select_df[~(final_select_df["cohort"] == "AA")]
sns.catplot(x = 'pval_names', y = 'snp_count', hue = 'Group', data = final_select_df,
            palette = 'hls',
            col="cohort",
            kind="bar",
            height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
#plt.legend()
plt.savefig(join(plot_path, '4_CombinedplotCohorts_EURsnpcount.pdf'))
plt.close()


##########################################
##########################################
## Related to plot1 (perSNP heritability):
"""BMI"""
g = sns.barplot(x = 'pval_names', y = 'R2perSNP', hue = 'Category', data = final_select_df,
            palette = 'hls',
            #col="Cohorts",
            #order = ['aa_filt', 'eur_filt', 'aa_unfilt', 'eur_unfilt', 'aa_weighted', 'eur_weighted', 'aa_gwasbaseline', 'eur_gwasbaseline'],  
            capsize = 0.05,             
            saturation = 8,             
            errcolor = 'gray', errwidth = 2,  
            ci = 'sd'   
            )

box = g.get_position()
g.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # resize position

# Put a legend to the right side
g.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)

plot_path=input_dir
plt.savefig(join(plot_path, '1.1_plot_test.pdf'))
plt.close()





# """PRS bar plots"""
# input_dir = "/Users/suryachhetri/datasets/prs_project/final_hg19/EBV_LCL/ldpredfunct_prs_output"


# """eQTL-informed - pearson approach"""
# id1 = "0.1_AA"; id2 = "AA"
# aa_filt0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")

# aa_filt1 = pd.Series([282, 482, 1444, 3252, 9050, 31948, 82822, 113465], name="snp_count")
# aa_filt = pd.concat([aa_filt0, pd.DataFrame(aa_filt1)], axis=1)
# aa_filt["cohort"] = "AFR"

# id1 = "0.1_EUR"; id2 = "EUR"
# eur_filt0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
# eur_filt1 = pd.Series([283, 483, 1450, 3269, 9114, 32107, 83218, 113957], name="snp_count")
# eur_filt = pd.concat([eur_filt0, pd.DataFrame(eur_filt1)], axis=1)
# eur_filt["cohort"] = "EUR"

# """eQTL-informed - bayesian approach""" #coloc-genes)
# # id1 = "colocfilt0.5_AA"; id2 = "AA"
# # aa_weighted0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
# # aa_weighted1 = pd.Series([375, 657, 2092, 4925, 14521, 54374, 144821, 199633], name="snp_count")
# # aa_weighted = pd.concat([aa_weighted0, pd.DataFrame(aa_weighted1)], axis=1)
# # aa_weighted["cohort"] = "AFR"

# # id1 = "colocfilt0.5_EUR"; id2 = "EUR"
# # eur_weighted0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
# # eur_weighted1 = pd.Series([376, 657, 2096, 4944, 14604, 54644, 145548, 200577], name="snp_count")
# # eur_weighted = pd.concat([eur_weighted0, pd.DataFrame(eur_weighted1)], axis=1)
# # eur_weighted["cohort"] = "EUR"

# """eQTL-informed - bayesian approach""" #coloc-genes)
# id1 = "colocfilt0.1_AA"; id2 = "AA"
# aa_weighted0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
# aa_weighted1 = pd.Series([257, 458, 1468, 3405, 9814, 36272, 95674, 131639], name="snp_count")
# aa_weighted = pd.concat([aa_weighted0, pd.DataFrame(aa_weighted1)], axis=1)
# aa_weighted["cohort"] = "AFR"

# id1 = "colocfilt0.1_EUR"; id2 = "EUR"
# eur_weighted0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
# eur_weighted1 = pd.Series([258, 458, 1471, 3412, 9871, 36398, 96022, 132110], name="snp_count")
# eur_weighted = pd.concat([eur_weighted0, pd.DataFrame(eur_weighted1)], axis=1)
# eur_weighted["cohort"] = "EUR"

# """cis-genes"""
# id1 = "unfilt_AA"; id2 = "AA"
# aa_unfilt0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
# aa_unfilt1 = pd.Series([414, 724, 2389, 5666, 16884, 63837, 170515, 235009], name="snp_count")
# aa_unfilt = pd.concat([aa_unfilt0, pd.DataFrame(aa_unfilt1)], axis=1)
# aa_unfilt["cohort"] = "AFR"

# id1 = "unfilt_EUR"; id2 = "EUR"
# eur_unfilt0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
# eur_unfilt1 = pd.Series([415, 725, 2396, 5691, 16971, 64146, 171316, 236020], name="snp_count")
# eur_unfilt = pd.concat([eur_unfilt0, pd.DataFrame(eur_unfilt1)], axis=1)
# eur_unfilt["cohort"] = "EUR"

# """gwas-baseline"""
# id1 = "gwasbaseline_AA"; id2 = "AA"
# aa_gwasbaseline0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
# aa_gwasbaseline1 = pd.Series([484, 873, 2981, 7411, 23481, 93599, 259178, 363755], name="snp_count")
# aa_gwasbaseline = pd.concat([aa_gwasbaseline0, pd.DataFrame(aa_gwasbaseline1)], axis=1)
# aa_gwasbaseline["cohort"] = "AFR"

# id1 = "gwasbaseline_EUR"; id2 = "EUR"
# eur_gwasbaseline0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
# eur_gwasbaseline1 = pd.Series([500, 913, 3305, 8717, 30363, 133546, 394937, 568503], name="snp_count")
# eur_gwasbaseline = pd.concat([eur_gwasbaseline0, pd.DataFrame(eur_gwasbaseline1)], axis=1)
# eur_gwasbaseline["cohort"] = "EUR"

# # concat dataframes
# frames = [aa_filt, eur_filt, aa_weighted, eur_weighted, aa_unfilt, eur_unfilt, aa_gwasbaseline, eur_gwasbaseline]
# result = pd.concat(frames, keys=['aa_filt', 'eur_filt', 'aa_weighted', 'eur_weighted', 'aa_unfilt', 'eur_unfilt',  'aa_gwasbaseline', 'eur_gwasbaseline'])

# final_df = result.reset_index().drop(columns="level_1", axis=0)
# final_df = final_df.rename(columns={"level_0": "Category"})
# final_df = final_df.loc[~final_df["Category"].str.contains("weighted")]


# """BMI"""
# #bmi_df1 = bmi_df[~(bmi_df["Cohorts"] == "ALL")]

# g = sns.barplot(x = 'pval_names', y = 'r2', hue = 'Category', data = final_df,
#             palette = 'hls',
#             #col="Cohorts",
#             #order = ['aa_filt', 'eur_filt', 'aa_unfilt', 'eur_unfilt', 'aa_weighted', 'eur_weighted', 'aa_gwasbaseline', 'eur_gwasbaseline'],  
#             capsize = 0.05,             
#             saturation = 8,             
#             errcolor = 'gray', errwidth = 2,  
#             ci = 'sd'   
#             )

# box = g.get_position()
# g.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # resize position

# # Put a legend to the right side
# g.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)

# plot_path = input_dir
# plt.savefig(join(plot_path, '1_plot_test.pdf'))
# plt.close()


#########################
#########################
# CLeaned Priority Plots 


"""PRS bar plots"""
input_dir = "/Users/suryachhetri/datasets/prs_project/final_hg19/EBV_LCL/ldpredfunct_prs_output"


"""gwas-baseline"""
#id1 = ""; id2 = "AA"
#aa_gwasbaseline0 = pd.read_csv(join(input_dir, tissue + "_prs_format_" + id1 + ".full_prs.scores_" + id2 + "_SBP_age_r2counts.tsv"), sep="\t")
aa_gwasbaseline0 = pd.Series([0.0008, 0.0010, 0.0010, 0.0010,  0.0002], name="raweffect_prs")
aa_gwasbaseline1 = pd.Series([0.0020, 0.0022, 0.0020, 0.0017, 0.0003], name="functprior_prs")
aa_gwasbaseline2 = pd.Series([182088, 147771, 128724, 113071, 47293], name="snp_count")
aa_gwasbaseline3 = pd.Series([">0%", ">25%", ">50%", ">75%", ">99%"], name="coloc_evidence")
aa_gwasbaseline = pd.concat([pd.DataFrame(aa_gwasbaseline0), pd.DataFrame(aa_gwasbaseline1),
                              pd.DataFrame(aa_gwasbaseline2), pd.DataFrame(aa_gwasbaseline3)], axis=1)
aa_gwasbaseline["cohort"] = "AFR"
aa_gwasbaseline["raweffect-R2perSNP"] = aa_gwasbaseline["raweffect_prs"]/aa_gwasbaseline["snp_count"].astype(float)
aa_gwasbaseline["functprior-R2perSNP"] = aa_gwasbaseline["functprior_prs"]/aa_gwasbaseline["snp_count"].astype(float)

#id1 = "gwasbaseline_EUR"; id2 = "EUR"
eur_gwasbaseline0 = pd.Series([0.0045, 0.0039, 0.0033, 0.0031,  0.0026], name="raweffect_prs")
eur_gwasbaseline1 = pd.Series([0.0117, 0.0112, 0.0089,  0.0077, 0.0047], name="functprior_prs")
eur_gwasbaseline2 = pd.Series([184177, 149648, 130401, 114614, 47667], name="snp_count")
eur_gwasbaseline3 = pd.Series([">0%", ">25%", ">50%", ">75%", ">99%"], name="coloc_evidence")
eur_gwasbaseline = pd.concat([pd.DataFrame(eur_gwasbaseline0), pd.DataFrame(eur_gwasbaseline1),
                              pd.DataFrame(eur_gwasbaseline2), pd.DataFrame(eur_gwasbaseline3)], axis=1)
eur_gwasbaseline["cohort"] = "EUR"
eur_gwasbaseline["raweffect-R2perSNP"] = eur_gwasbaseline["raweffect_prs"]/eur_gwasbaseline["snp_count"].astype(float)
eur_gwasbaseline["functprior-R2perSNP"] = eur_gwasbaseline["functprior_prs"]/eur_gwasbaseline["snp_count"].astype(float)

concat_df = pd.concat([aa_gwasbaseline, eur_gwasbaseline])

#df_select = concat_df.loc[:,["coloc_evidence", "cohort", "raweffect-R2perSNP", "functprior-R2perSNP"]]
pd.melt(concat_df, id_vars=["coloc_evidence", "cohort"], value_vars=["raweffect-R2perSNP", "functprior-R2perSNP"])

# #final_select_df = final_df.loc[~final_df["Category"].str.contains("weighted")]
# final_select_df = final_df.copy()
# final_select_df["Group"] = ""
# final_select_df.loc[final_select_df["Category"].str.contains("_filt"), "Group"] = "eQTL-informed-genes(pearson)"
# final_select_df.loc[final_select_df["Category"].str.contains("_weighted"), "Group"] = "eQTL-informed-genes(bayesian)"
# final_select_df.loc[final_select_df["Category"].str.contains("_unfilt"), "Group"] = "cis-genes"
# final_select_df.loc[final_select_df["Category"].str.contains("gwasbaseline"), "Group"] = "gwas-baseline"

# final_select_df["R2perSNP"] = final_select_df["r2"]/final_df["snp_count"].astype(float)

#final_df1 = final_df[final_df["Category"] == "eQTL Informed"]
sns.catplot(x = 'coloc_evidence', y = 'raweffect-R2perSNP', hue = 'cohort', data = concat_df,
            palette = 'hls',
            col="cohort", #remove this for adjacent barplot
            kind="bar",
            height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )

#plt.legend()
plot_path = input_dir
plt.savefig(join(plot_path, '2.1_Plot_CombinedplotCohorts.pdf'))
plt.close()

#grouped sided together bar plot
sns.barplot(x = 'coloc_evidence', y = 'raweffect-R2perSNP', hue = 'cohort', data = concat_df,
            palette = 'hls',
            #col="cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
plot_path = input_dir
plt.savefig(join(plot_path, '2.barPlot_CombinedplotCohorts.pdf'))
plt.close()


#grouped sided together bar plot
sns.barplot(x = 'coloc_evidence', y = 'functprior-R2perSNP', hue = 'cohort', data = concat_df,
            palette = 'hls',
            #col="cohort",
            #kind="bar",
            #height=7, aspect=1.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
plot_path = input_dir
plt.savefig(join(plot_path, '2.functprior_barPlot_CombinedplotCohorts.pdf'))
plt.close()


#df_select = concat_df.loc[:,["coloc_evidence", "cohort", "raweffect-R2perSNP", "functprior-R2perSNP"]]
df_final = pd.melt(concat_df, id_vars=["coloc_evidence", "cohort"], value_vars=["raweffect-R2perSNP", "functprior-R2perSNP"])


#final_select_df1 = final_select_df[~(final_select_df["cohort"] == "AA")]
sns.catplot(x = 'coloc_evidence', y = 'value', hue = 'variable', data = df_final,
            #palette = 'hls',
            col="cohort",
            kind="bar",
            height=9, aspect=2.5,
            #order = ['EA', 'AA', 'ALL'],  
            #capsize = 0.05,             
            #saturation = 8,             
            #errcolor = 'gray', errwidth = 2,  
            #ci = 'sd'   
            )
#plt.legend()
plt.savefig(join(plot_path, '4_raw-functprior-CombinedplotCohorts_EURsnpcount.pdf'))
plt.close()


# ##########################################
# ##########################################
# ## Related to plot1 (perSNP heritability):
# """BMI"""
# g = sns.barplot(x = 'pval_names', y = 'R2perSNP', hue = 'Category', data = final_select_df,
#             palette = 'hls',
#             #col="Cohorts",
#             #order = ['aa_filt', 'eur_filt', 'aa_unfilt', 'eur_unfilt', 'aa_weighted', 'eur_weighted', 'aa_gwasbaseline', 'eur_gwasbaseline'],  
#             capsize = 0.05,             
#             saturation = 8,             
#             errcolor = 'gray', errwidth = 2,  
#             ci = 'sd'   
#             )

# box = g.get_position()
# g.set_position([box.x0, box.y0, box.width * 0.95, box.height]) # resize position

# # Put a legend to the right side
# g.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)

# plot_path=input_dir
# plt.savefig(join(plot_path, '1.1_plot_test.pdf'))
# plt.close()

