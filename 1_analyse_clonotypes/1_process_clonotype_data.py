'''Clonotype analyses
'''

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob
import re
import math
import collections
import seaborn as sns
import matplotlib_venn
import scipy
import statsmodels.api as sm
import itertools
import functools

# Import (or reload) helper functions
import stackedBarGraph
import imp
import process_clonotype_data_helpers
imp.reload(process_clonotype_data_helpers)
from process_clonotype_data_helpers import *

# Set interactive plotting mode
plt.ion()
#
# Define clonotype field names
# - V gene, J gene, CDR3 length
clonotype_filter = ["V-GENE", "J-GENE", "CDR3_len"]
def clonotype_format(x): return "{}.{}.CDR3_len{}".format(*x)
#
# Get list of sample clonotypes
summary_files = list(glob.iglob(
    "/nfs/users/nfs_b/bb9/workspace/rotation2/team115_lustre/data/existing_bcr_data/MalariaSamplesBen/**/*_Summary.txt", 
    recursive=True, 
)) + ["/nfs/users/nfs_b/bb9/workspace/rotation2/team115_lustre/data/2017-02-06_IMGT_LEA_S20/Day63/LEA_S20_Summary.txt"]
#
# Extract sample day and number from filename
p = re.compile(r'Day(\d+)/LEA_S(\d+)_Summary.txt')
days = [p.search(f).group(1) for f in summary_files]
sample_nums = [p.search(f).group(2) for f in summary_files]
#
# Read in all summary files
summary_df = pd.concat(
    (pd.read_table(f, header=0, usecols=list(range(0, 30))) for f in summary_files),
    keys=zip(days, sample_nums),
    names=["day", "sample_num"]
)
summary_df.index.names = pd.core.base.FrozenList(['day', 'sample_num', 'seq_i'])
#
# Count unique BCRs per sample
#  print(summary_df.groupby(level="sample_num").apply(lambda x: x.shape))
# List possible Functionalities
#  print(summary_df["Functionality"].unique())
# Filter out unproductive and "No results"
summary_df = summary_df[summary_df["Functionality"].isin(['productive', 'productive (see comment)'])]
# Count unique BCRs per sample after filter
# print(summary_df.groupby(level="sample_num").apply(lambda x: x.shape))
#
# Convert row indices into columns
summary_df = summary_df.reset_index(level=summary_df.index.names)
del summary_df["seq_i"]
# 
summary_df["CDR3_len"] = summary_df["CDR3-IMGT length"].astype(int) + 2 # Add two to include conserved C/W residues
summary_df["CDR3_aa"] = summary_df["AA JUNCTION"].apply(lambda x: x.split()[0])
summary_df[["day", "sample_num"]] = summary_df[["day", "sample_num"]].apply(pd.to_numeric)
# Get V J gene assignments by taking the first assignment, collapsing alleles.
summary_df["V-GENE"] = [x.split(", or ")[0].split("*")[0].replace("Homsap ", "") for x in summary_df["V-GENE and allele"]]
summary_df["J-GENE"] = [x.split(", or ")[0].split("*")[0].replace("Homsap ", "") for x in summary_df["J-GENE and allele"]]
#
# Get clonotypes for known mAbs
mAb_df = pd.read_csv("/nfs/users/nfs_b/bb9/workspace/rotation2/team115_lustre/data/existing_bcr_data/MalariaSamplesBen/IgBlast_bnAbs.csv")
mAb_df["V-GENE"] = [x.split(",")[0].split("*")[0] for x in mAb_df["V gene"]]
mAb_df["J-GENE"] = [x.split(",")[0].split("*")[0] for x in mAb_df["J gene"]]
mAb_df["VJ-GENE"] = mAb_df["V-GENE"] + "." + mAb_df["J-GENE"]
mAb_df["CDR3_aa"] = mAb_df["CDR3 amino acid seq"]
mAb_df["CDR3_len"] = mAb_df["CDR3 amino acid seq"].map(len)
mAb_df["clonotype"] = mAb_df[clonotype_filter].apply(clonotype_format, axis=1)
#
# Get clonotypes for existing data
summary_df["clonotype"] = summary_df[clonotype_filter].apply(clonotype_format, axis=1)
summary_df["VJ-GENE"] = summary_df["V-GENE"] + "." + summary_df["J-GENE"]
#
# Read in sample information
sample_info_df = pd.read_excel("/nfs/users/nfs_b/bb9/workspace/rotation2/team115_lustre/data/existing_bcr_data/MalariaSamplesBen/Malaria_Samples_SeqInfo.xlsx")
#
# Merge in patient number and cell type
summary_df = pd.merge(summary_df, sample_info_df[["Tag Index", "patient_code", "cell_type"]], 
         how="left", left_on="sample_num", right_on="Tag Index")
#
# Merge correspondance to known mAbs
# summary_df = pd.merge(summary_df, mAb_df[["Ab.Name", "clonotype"]], how="left", on="clonotype")
#
# Check if V-gene is in the set of V-genes used in known mAbs
# summary_df["V-GENE_known_usage"] = np.where(summary_df["V-GENE"].isin(mAb_df["V-GENE"]), summary_df["V-GENE"], "NaN")
#
# Mark clones with multiple isotypes
summary_df["isotypes"] = summary_df["digest"].map(lambda x: x.split("/"))
summary_df["n_isotypes"] = summary_df["isotypes"].map(len)
#
# Calculate mutational freq per base pair
summary_df["V-REGION_unmut_len"] = summary_df["V-REGION identity nt"].map(lambda x: int(x.split()[0].split('/')[0]))
summary_df["J-REGION_unmut_len"] = summary_df["J-REGION identity nt"].map(lambda x: int(x.split()[0].split('/')[0]))
summary_df["VJ-REGION_unmut_len"] = summary_df["V-REGION_unmut_len"] + summary_df["J-REGION_unmut_len"]
#
summary_df["V-REGION_len"] = summary_df["V-REGION identity nt"].map(lambda x: int(x.split()[0].split('/')[1]))
summary_df["J-REGION_len"] = summary_df["J-REGION identity nt"].map(lambda x: int(x.split()[0].split('/')[1]))
summary_df["VJ-REGION_len"] = summary_df["V-REGION_len"] + summary_df["J-REGION_len"]
#
summary_df["mut_freq_v"] = (
    (summary_df["V-REGION_len"] - summary_df["V-REGION_unmut_len"])
    /summary_df["V-REGION_len"]
)
summary_df["mut_freq_j"] = (
    (summary_df["J-REGION_len"] - summary_df["J-REGION_unmut_len"])
    /summary_df["J-REGION_len"]
)
summary_df["mut_freq_vj"] = (
    (summary_df["VJ-REGION_len"] - summary_df["V-REGION_unmut_len"] - summary_df["J-REGION_unmut_len"])
    / summary_df["VJ-REGION_len"]
)
# Change to categorical patient code
summary_df["patient_code"] = summary_df["patient_code"].astype("category")

#
# Write out counts, read in normalised counts from edgeR
#
# Write out summary_df
# summary_df.to_csv("../team115_lustre/1_analyse_clonotypes/summary.csv", index=False)
#
# Biological unit: clonotype
# count_clonotype = get_clonotype_freq(summary_df, "clonotype", groups=["sample_num"])
# count_clonotype.to_csv("../team115_lustre/1_analyse_clonotypes/count_clonotype.csv")
#
# Biological unit: vj gene combo
# count_vj = get_clonotype_freq(summary_df, "VJ-GENE", groups=["sample_num"])
# count_vj.to_csv("../team115_lustre/1_analyse_clonotypes/count_vj.csv")
#
# Biological unit: v gene
# count_v = get_clonotype_freq(summary_df, "V-GENE", groups=["sample_num"])
# count_v.to_csv("../team115_lustre/1_analyse_clonotypes/count_v.csv")
#
# Write out sample information csv
# sample_info_df.to_csv("../team115_lustre/1_analyse_clonotypes/sample_info.csv", index=False)
#
# Read back in normalised counts
# count_clonotype_cpm_normed = pd.read_csv("./.output/count_clonotype.cpm_normed.csv", header=0, index_col=0)

#
# Two timepoint analysis
# d0 MBC vs d140 MBC
#
rep1 = summary_df.query("day == 0 & cell_type == 'MBC'")
rep2 = summary_df.query("day == 140 & cell_type == 'MBC'")
patient_codes = sorted(summary_df["patient_code"].unique())
# analysis_level = "VJ-GENE"
analysis_level = "clonotype"
# analysis_level = "CDR3_aa"

#
# Specific expansions of clonotypes
#
#
# Detect increases in clonotype proportion
# Z-test for difference in proportions
#
rep1_freq = get_clonotype_freq(rep1, analysis_level, ["patient_code"])
rep2_freq = get_clonotype_freq(rep2, analysis_level, ["patient_code"])
#
rep1_prop = rep1_freq.apply(lambda x: x/sum(x))
rep2_prop = rep2_freq.apply(lambda x: x/sum(x))
#
rep_freq_joined = rep1_freq.join(rep2_freq, how='outer', lsuffix="_rep1", rsuffix="_rep2").fillna(0)
rep_prop_joined = rep1_prop.join(rep2_prop, how='outer', lsuffix="_rep1", rsuffix="_rep2").fillna(0)
#

# TODO
# Look at distribution of counts to determine possible low counts threshold
# rep_freq_joined_long = pd.melt(rep_freq_joined.apply(np.sqrt))
# fig, ax = plt.subplots()
# for col in rep_freq_joined.columns:
    # sns.kdeplot(np.sqrt(rep_freq_joined[col]), ax=ax)

# Detect expanded clonotypes for each patient,
# then perform multipletesting correction on all samples simultaneously
clonotypes_expanded_dict = dict()
for patient_code in patient_codes: 
    x1 = rep_freq_joined[str(patient_code) + "_rep1"]
    x2 = rep_freq_joined[str(patient_code) + "_rep2"]
    zs, ps = prop_test_ztest(
        x1, x2,
        alternative="two-sided"
    )
    clonotypes_expanded_dict[patient_code] = pd.DataFrame({"clonotype": rep_freq_joined.index, 
                                                           "z": zs,
                                                           "p": ps,
                                                           "rep1_freq": x1,
                                                           "rep2_freq": x2,
                                                           "rep1_prop": x1/sum(x1),
                                                           "rep2_prop": x2/sum(x2),
                                                           "known_mAb": rep_freq_joined.index.map(lambda x: list(mAb_df.loc[mAb_df[analysis_level] == x, "Ab.Name"]))
                                                           })
# Merge and perform multipletesting correction
clonotypes_expanded_df = pd.concat(
        clonotypes_expanded_dict.values(), 
        keys=list(clonotypes_expanded_dict.keys()), 
        names=["patient_code"])
clonotypes_expanded_df = clonotypes_expanded_df.reset_index(level=0)
clonotypes_expanded_df["signif"], clonotypes_expanded_df["p_corrected"], _, _ = sm.stats.multipletests(
        clonotypes_expanded_df["p"], method="holm", alpha=0.05)
# Determine expansion based on sign of z statistic
clonotypes_expanded_df["expanded"] = clonotypes_expanded_df["signif"] & (clonotypes_expanded_df['z'] < 0)
# Change patient code to categorical
clonotypes_expanded_df['patient_code'] = clonotypes_expanded_df['patient_code'].astype("category", ordered=True)
# Subset out those that are expanded in each patient
clonotypes_expanded = dict()
for patient_code in patient_codes:
    clonotypes_expanded[patient_code] = clonotypes_expanded_df.loc[clonotypes_expanded_df["patient_code"] == patient_code].query('expanded')

# TODO
if 0:
    #
    # For each clonotype, get status:
    #
    results = pd.DataFrame({
        "clonotype": rep_freq_joined.index 
    })
    results["expanded"] = results["clonotype"].apply(lambda x:
        list(filter(lambda p: x in clonotypes_expanded[p], patient_codes))
    )
    naive_reps_d0 = dict(zip(
        patient_codes,
        (set(get_naive_rep(summary_df).query("day == 0 & patient_code == {}".format(p))[analysis_level]) for p in patient_codes)
    )) 
    pbmc_reps_d0 = dict(zip(
        patient_codes,
        (set(summary_df.query("cell_type == 'PBMCs' & day == 0 & patient_code == {}".format(p))[analysis_level]) for p in patient_codes)
    )) 
    results["in_naive_d0"] = results["clonotype"].apply(lambda x:
        list(filter(lambda p: x in naive_reps_d0[p], patient_codes))
    )
    results["in_pbmc_d0"] = results["clonotype"].apply(lambda x:
        list(filter(lambda p: x in pbmc_reps_d0[p], patient_codes))
    )
    results["known_mAb_analysis_level"] = results["clonotype"].apply(lambda x:
        list(mAb_df.loc[mAb_df[analysis_level] == x, "Ab.Name"])
    )
    # results["in_MBC_d0"]
    # results["in_MBC_d140"]
    results["isotype_rep1"] = results["clonotype"].apply(lambda x:
        collections.Counter(
            itertools.chain(*
                map(
                    lambda x: dict((i, 1/len(x)) for i in x), # Split into fractional counts for multi-isotype clones
                    rep1.loc[rep1[analysis_level] == x, "isotypes"]
                )
            )
        )
    )
    results["isotype_rep2"] = results["clonotype"].apply(lambda x:
        collections.Counter(
            itertools.chain(*
                map(
                    lambda x: dict((i, 1/len(x)) for i in x), # Split into fractional counts for multi-isotype clones
                    rep2.loc[rep2[analysis_level] == x, "isotypes"]
                )
            )
        )
    )

    results[results.loc[:, ["expanded"]].apply(any, axis=1)]

#
# Do expanded clonotypes have a higher mutational frequency on average?
# Plot mutational freqs of clones, group by patient, hue by expansion
#
data = rep2
data = data.assign(expanded=
    data.apply(lambda x: x[analysis_level] in set(clonotypes_expanded[x['patient_code']][analysis_level]), axis=1)
)
fig, ax = plt.subplots()
sns.violinplot(x="mut_freq_vj", y="patient_code", hue="expanded", data=data, ax=ax, orient="h")
# sns.stripplot(x="mut_freq_vj", y="patient_code", hue="expanded", 
        # data=data[data[analysis_level].isin(mAb_df[analysis_level])], 
        # ax=ax, orient="h", split=True)
# Perform mann whitney u test for difference in mean rank
title_parts = []
for patient_code in patient_codes:
    if len(clonotypes_expanded[patient_code]):
        _, p_mw = scipy.stats.mannwhitneyu(
            data.loc[(data['patient_code'] == patient_code) & ~data['expanded'], "mut_freq_vj"],
            data.loc[(data['patient_code'] == patient_code) & data['expanded'], "mut_freq_vj"])
        rank_ne, rank_e, med_ne, med_e = get_mean_ranks(
            data.loc[(data['patient_code'] == patient_code) & ~data['expanded'], "mut_freq_vj"],
            data.loc[(data['patient_code'] == patient_code) & data['expanded'], "mut_freq_vj"])
        # title_parts.append(" ".join(str(x) for x in [patient_code, p_mw, rank_ne, rank_e, med_ne, med_e]))
        title_parts.append(" ".join(str(x) for x in ["patient: ", patient_code, "p = ", p_mw]))
ax.set_title("\n".join(title_parts))
fig.savefig("../team115_lustre/1_analyse_clonotypes/expanded_clonotypes_mut_freq.pdf")

#
# Manhattan plot, p value for expansion
# 
clonotypes_expanded_df['is_known_mAb_clonotype'] = clonotypes_expanded_df['known_mAb'].astype(bool)
#
# Determine significance threshold line
signif = clonotypes_expanded_df["signif"]
max_positive = max(clonotypes_expanded_df["p"][signif])
min_negative = min(clonotypes_expanded_df["p"][~signif])
thresh = np.mean([-np.log10(max_positive), -np.log10(min_negative)])
#
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True, gridspec_kw={'height_ratios':[1, 10]})
# Plot same dataset on two subplots
jitter_rad = 0.35
g1 = sns.stripplot(y=-np.log10(clonotypes_expanded_df["p_corrected"]), x="patient_code",
        data=clonotypes_expanded_df, ax=ax1, jitter=jitter_rad, 
        hue="is_known_mAb_clonotype", split=True)
g2 = sns.stripplot(y=-np.log10(clonotypes_expanded_df["p_corrected"]), x="patient_code",
        data=clonotypes_expanded_df, ax=ax2, jitter=jitter_rad, 
        hue="is_known_mAb_clonotype", split=True)
# Zoom in on separate parts of dataset
ax1.set_ylim([210, 220])
ax1.xaxis.set_visible(False) # Remove x axis labels from top subplot
ax1.set_yticks([210, 220]) # Keep tick spacing consistent between subplots
ax1.set_ylabel("") # Remove y axis label
ax2.set_ylim([0, 90])
ax2.legend().set_visible(False)
# Add signif line
ax2.axhline(y=thresh, color="black", linestyle="--")
fig.tight_layout()
fig.savefig("../team115_lustre/1_analyse_clonotypes/expand_p_val_manhattan.pdf")

#
# Show proportion of repertoire taken up by certain clonotypes
# e.g. known mAb clonotypes
#
# Choose set of clonotypes to plot
clonotypes_to_plot = mAb_df[analysis_level].unique()
props = rep_prop_joined[rep_prop_joined.index.isin(clonotypes_to_plot)].transpose()
# Sort samples such that patients are adjacent
props = props.loc[props.index.sort_values()]
# Sort clonotypes by max prop in repertoire
props = props.loc[:, props.apply(sum).sort_values(ascending=False).index]
# Plot stacked bars
sbg = stackedBarGraph.StackedBarGrapher()
fig, ax = plt.subplots(figsize=(10, 10))
stack_colors = plt.cm.Set2(np.linspace(0, 1, props.shape[1]))
sbg.stackedBarPlot(
    ax,
    props,
    stack_colors,
    props.index, 
    gap=0.1,
    scale=False,
    xlabel="Sample repertoire",
    ylabel="Proportion of repertoire"
)
# Remove vertical gridlines
ax.xaxis.grid(False) 
# add lines separating patients
plt.axvline(x=1.5, color="w", linestyle="-")
plt.axvline(x=3.5, color="w", linestyle="-")
# Add legend, showing samples in which expansion occurred
legends = []
i = 0
for column in props.columns: 
    label = column
    for patient_code in patient_codes:
        if label in set(clonotypes_expanded[patient_code][analysis_level]):
            label = label + ", Expanded in: " + str(patient_code)
    legends.append(matplotlib.patches.Patch(color=stack_colors[i], label=label))
    i+=1
plt.legend(handles=legends, loc='best')
#
fig.savefig("../team115_lustre/1_analyse_clonotypes/clonotype_stacked_bars.pdf")

#
# Clonotype sharing, indicative of convergent response
#
# Sharing of expanded clonotypes: venn diagram
#
fig, ax = plt.subplots(figsize=(10, 10))
# (100, 010, 110, 001, 101, 011, 111)
v = matplotlib_venn.venn3(subsets=(1, 1, 0, 1, 0, 0, 0), set_labels=patient_codes)
# for x in v.subset_labels:
    # x.set_text("0")
v.get_label_by_id('100').set_text(len(clonotypes_expanded[patient_codes[0]]))
v.get_label_by_id('010').set_text(len(clonotypes_expanded[patient_codes[1]]))
v.get_label_by_id('001').set_text(len(clonotypes_expanded[patient_codes[2]]))
#
# plt.annotate(" \n ", 
    # xy=v.get_label_by_id('100').get_position(), xytext=(0, -120), size=11,
    # ha='center', va='top', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1)
# )
plt.annotate("\n".join(list(clonotypes_expanded[patient_codes[1]][analysis_level].sort_values())), 
    xy=v.get_label_by_id('010').get_position(), xytext=(0, -120), size=11,
    ha='center', va='top', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1)
)
plt.annotate("\n".join(list(clonotypes_expanded[patient_codes[2]][analysis_level].sort_values())), 
    xy=v.get_label_by_id('001').get_position(), xytext=(0, -120), size=11,
    ha='center', va='top', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1)
)
plt.title("Venn diagram of expanded clonotypes per patient")
fig.savefig("../team115_lustre/1_analyse_clonotypes/expanded_clonotype_venn.pdf")

# TODO explore certain clonotypes

# IGHV3-33.IGHJ4.CDR3_len16

rep_freq_joined.query("clonotype == 'IGHV3-23.IGHJ4.CDR3_len16'")

rep_freq_joined.query("clonotype == 'IGHV4-39.IGHJ4.CDR3_len15'")

# IGHV4-39.IGHJ4.CDR3_len15
rep_freq_joined.query("clonotype == 'IGHV1-18.IGHJ5.CDR3_len16'")

# ab r5
'IGHV1-18.IGHJ5.CDR3_len16' in set(mAb_df['clonotype'])

rep_freq_joined.query("clonotype == 'IGHV1-18.IGHJ5.CDR3_len16'")

clonotypes_expanded_ps_combined.query("clonotype == 'IGHV1-18.IGHJ5.CDR3_len16'")

rep_freq_joined.apply(lambda x: x/sum(x)).apply(max, axis=1).sort_values().tail()

rep_freq_joined.apply(lambda x: x/sum(x)).sort_values(["1019_rep2"]).tail()

rep_freq_joined.apply(lambda x: x/sum(x)).sort_values(["2207_rep2"]).tail()

rep_freq_joined.apply(lambda x: x/sum(x)).query("clonotype == 'IGHV1-18.IGHJ5.CDR3_len16'")

# TODO Heatmap of clonotype sharing

# Ordering with timepoints grouped
# sample_names = ["LEA_S" + str(n) for n in [1, 11, 15, 3, 9, 7]]
# Ordering with patients grouped
sample_names = ["LEA_S" + str(n) for n in [1, 3, 11, 9, 15, 7]]
sample_names = ['1017_rep1', '1019_rep1', '2207_rep1', '1017_rep2', '1019_rep2', "2207_rep2"]
clonotype_cpm_normed_props = rep_freq_joined.loc[:, sample_names].apply(lambda x: x/sum(x))
# Distance measure: number of clonotypes present in both repertoires
dist = scipy.spatial.distance.pdist(
    clonotype_cpm_normed_props.transpose(), 
    lambda x, y: 
        sum((x[i]+y[i])/(sum(x)+sum(y)) for i in range(len(x)) if x[i] and y[i])
)

dist = pd.DataFrame(scipy.spatial.distance.squareform(dist))
dist.columns = sample_names
dist.index = sample_names
# Generate a triangular mask 
mask = np.ones_like(dist, dtype=np.bool)
mask[np.tril_indices_from(mask)] = False
fig = plt.figure()
g = sns.heatmap(dist, annot=True, fmt=".2f", linewidths=.5, square=True, mask=mask)
fig.savefig("../team115_lustre/1_analyse_clonotypes/clonotype_sharing_heatmap.pdf")

# TODO 
# Pies of isotype freq by mutation and expansion
# convert to Pies of isotype freq by mutation and known

# Em_clono = results[(results["expanded"].astype(bool) & np.invert(results["mutated"].astype(bool)))]
# eM_clono = results[(results["mutated"].astype(bool) & np.invert(results["expanded"].astype(bool)))]
# EM_clono = results[(results["mutated"].astype(bool) & (results["expanded"].astype(bool)))]
# em_clono = results[np.invert(results["mutated"].astype(bool) & np.invert(results["expanded"].astype(bool)))]
Em_clono = results[(results["expanded"].astype(bool) & (~results["known_mAb_analysis_level"].astype(bool)))]
eM_clono = results[(results["known_mAb_analysis_level"].astype(bool) & (~results["expanded"].astype(bool)))]
EM_clono = results[(results["known_mAb_analysis_level"].astype(bool) & (results["expanded"].astype(bool)))]
em_clono = results[(~results["known_mAb_analysis_level"].astype(bool) & (~results["expanded"].astype(bool)))]

# Em_clono_sum = functools.reduce(lambda x, y: x+y, list(Em_clono["isotype"]))
# eM_clono_sum = functools.reduce(lambda x, y: x+y, list(eM_clono["isotype"]))
# EM_clono_sum = functools.reduce(lambda x, y: x+y, list(EM_clono["isotype"]))
# em_clono_sum = functools.reduce(lambda x, y: x+y, list(em_clono["isotype"]))
Em_clono_sum = functools.reduce(lambda x, y: x+y, list(Em_clono["isotype_rep2"]))
eM_clono_sum = functools.reduce(lambda x, y: x+y, list(eM_clono["isotype_rep2"]))
EM_clono_sum = functools.reduce(lambda x, y: x+y, list(EM_clono["isotype_rep2"]))
em_clono_sum = functools.reduce(lambda x, y: x+y, list(em_clono["isotype_rep2"]))

# Keep colors consistent
isotypes = (Em_clono_sum+ eM_clono_sum+ EM_clono_sum+ em_clono_sum).keys()
cols = sns.color_palette("hls", n_colors=len(isotypes))
cols_dict = dict(zip(isotypes, cols))

f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
ax1.axis("equal")
ax2.axis("equal")
ax3.axis("equal")
ax4.axis("equal")
ax3.pie(list(Em_clono_sum.values()), labels=Em_clono_sum.keys(), colors=[cols_dict[i] for i in Em_clono_sum.keys()])
ax2.pie(list(eM_clono_sum.values()), labels=eM_clono_sum.keys(), colors=[cols_dict[i] for i in eM_clono_sum.keys()])
ax1.pie(list(EM_clono_sum.values()), labels=EM_clono_sum.keys(), colors=[cols_dict[i] for i in EM_clono_sum.keys()])
ax4.pie(list(em_clono_sum.values()), labels=em_clono_sum.keys(), colors=[cols_dict[i] for i in em_clono_sum.keys()])
ax3.set_title("Em (n = {}, nKnown = {})".format(len(Em_clono), sum(Em_clono["known_mAb_analysis_level"].astype(bool))))
ax2.set_title("eM (n = {}, nKnown = {})".format(len(eM_clono), sum(eM_clono["known_mAb_analysis_level"].astype(bool))))
ax1.set_title("EM (n = {}, nKnown = {})".format(len(EM_clono), sum(EM_clono["known_mAb_analysis_level"].astype(bool))))
ax4.set_title("em (n = {}, nKnown = {})".format(len(em_clono), sum(em_clono["known_mAb_analysis_level"].astype(bool))))

f.savefig("../team115_lustre/1_analyse_clonotypes/pie.pdf")

# TODO
# Evaluate repertoire sharing at CDR3 seq level

# TODO graph proving lack of response in 1017 is not due to lack of clonotypes in naive repertoire
# perhaps add naive reps to heatmap

# TODO are shared clonotypes the same?
# Plot p value vs amount of sharing

