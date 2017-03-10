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
        title_parts.append(" ".join(str(x) for x in [patient_code, p_mw, rank_ne, rank_e, med_ne, med_e]))
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
# Clonotype sharing, indicative of convergent response
#

# TODO now : stacked bar
# TODO Two timepoint stacked area plot of clonotype proportions (top 100)

rep_prop_joined = rep_freq_joined.apply(lambda x: x/sum(x))
# Get clonotypes with the highest proportion of repertoire taken up at d140
top_prop_clonotypes = rep_freq_joined.apply(lambda x: x/sum(x)).loc[:, ['1017_rep2', '1019_rep2', '2207_rep2']].apply(sum, axis=1).sort_values().tail(20).index


plt.stackplot([0, 140], rep_prop_joined.reindex(top_prop_clonotypes).loc[:, ["1017_rep1", "1017_rep2"]])
# plt.ylim([0, 1])

plt.stackplot([0, 140], rep_prop_joined.reindex(top_prop_clonotypes).loc[:, ["2207_rep1", "2207_rep2"]])
# plt.ylim([0, 1])

plt.stackplot([1, 2, 3], rep_prop_joined.reindex(top_prop_clonotypes).loc[:, ['1017_rep2', '1019_rep2', '2207_rep2']])
# plt.ylim([0, 1])

plt.figure()
plt.stackplot([1, 2, 3], rep_prop_joined.reindex(mAb_df["clonotype"]).dropna().loc[:, ['1017_rep2', '1019_rep2', '2207_rep2']])
# plt.ylim([0, 1])
plt.savefig("../team115_lustre/1_analyse_clonotypes/known_mab_clonotype_stacked_area.pdf")

# Trying factorplot instead
foo = rep_prop_joined.reindex(mAb_df["clonotype"]).dropna().loc[:, ['1017_rep2', '1019_rep2', '2207_rep2']]
bar = pd.melt(foo.reset_index(), id_vars=["clonotype"])
g = sns.factorplot(x="clonotype", y="value", hue="patient_code", data=bar, kind="bar", palette="muted")



# TODO explore certain clonotypes

IGHV3-33.IGHJ4.CDR3_len16

rep_freq_joined.query("clonotype == 'IGHV3-23.IGHJ4.CDR3_len16'")

rep_freq_joined.query("clonotype == 'IGHV4-39.IGHJ4.CDR3_len15'")

IGHV4-39.IGHJ4.CDR3_len15
rep_freq_joined.query("clonotype == 'IGHV1-18.IGHJ5.CDR3_len16'")

# ab r5
'IGHV1-18.IGHJ5.CDR3_len16' in set(mAb_df['clonotype'])

rep_freq_joined.query("clonotype == 'IGHV1-18.IGHJ5.CDR3_len16'")

clonotypes_expanded_ps_combined.query("clonotype == 'IGHV1-18.IGHJ5.CDR3_len16'")



rep_freq_joined.apply(lambda x: x/sum(x)).apply(max, axis=1).sort_values().tail()

rep_freq_joined.apply(lambda x: x/sum(x)).sort_values(["1019_rep2"]).tail()

rep_freq_joined.apply(lambda x: x/sum(x)).sort_values(["2207_rep2"]).tail()

rep_freq_joined.apply(lambda x: x/sum(x)).query("clonotype == 'IGHV1-18.IGHJ5.CDR3_len16'")

# TODO
# Proportions of the repertoire taken up by known clonotypes
rep_freq_joined[rep_freq_joined.index.isin(mAb_df[analysis_level])].apply(sum) / rep_freq_joined.apply(sum)
rep_freq_joined[rep_freq_joined.index.isin(mAb_df[analysis_level])].apply(sum) / rep_freq_joined.apply(sum)

rep_freq_joined[rep_freq_joined.index.isin(mAb_df[analysis_level])]

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
# workspace below

# % TODO occuptation of rep by known clonotpye abs

# Paired bar plots of isotype proportions for clonotypes that are expanded and mutated, vs clonotypes that are not

isotype_freqs = foo.groupby('did_expand').apply(lambda x: functools.reduce(lambda x, y: x+y, x["isotype"]))

# Pies of isotype freq by mutation and expansion
# TODO 
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
# Boxplots: p value of clonotype by the number of samples a clonotype is present in.
# Do this at d0 and d140

sample_names = ["LEA_S" + str(n) for n in [1, 11, 15]]
sharing_count_d0 = count_clonotype_cpm_normed.loc[:, sample_names].apply(np.count_nonzero, axis=1)
sample_names = ["LEA_S" + str(n) for n in [3, 9, 7]]
sharing_count_d140 = count_clonotype_cpm_normed.loc[:, sample_names].apply(np.count_nonzero, axis=1)
foo = results.merge(
    pd.DataFrame({'sharing_count_d0': sharing_count_d0, 'sharing_count_d140': sharing_count_d140}),
    left_on='clonotype',
    right_index=True
)
foo = results.merge(
    clonotype_ps_combined.groupby('clonotype').apply(lambda x: np.min(x['p_expand'])).to_frame('mean_p_expand'),
    left_on='clonotype',
    right_index=True
)

def get_num_containing_samples(clonotype, rep, analysis_level="clonotype"):
    return len(rep.loc[rep[analysis_level] == clonotype, "patient_code"].unique())
sample_names = ["LEA_S" + str(n) for n in [1, 11, 15]]
sharing_count_d0 = np.array(count_clonotype_cpm_normed.loc[:, sample_names].index.map(lambda x: get_num_containing_samples(x, rep1, analysis_level)))
sample_names = ["LEA_S" + str(n) for n in [3, 9, 7]]
sharing_count_d140 = np.array(count_clonotype_cpm_normed.loc[:, sample_names].index.map(lambda x: get_num_containing_samples(x, rep1, analysis_level)))
sharing_count_combined = pd.DataFrame({'sharing_count_d0': sharing_count_d0, 'sharing_count_d140': sharing_count_d140})
sharing_count_combined.index = count_clonotype_cpm_normed.index
foo = results.merge(
    sharing_count_combined,
    left_on='clonotype',
    right_index=True
)
foo = foo.merge(
    clonotype_ps_combined.groupby('clonotype').apply(lambda x: np.mean(x['p_expand'])).to_frame('mean_p_expand'),
    left_on='clonotype',
    right_index=True
)
foo['did_expand'] = foo['expanded'].astype(bool)

sns.swarmplot(x="sharing_count_d0", y="mean_p_expand", data=foo, hue='did_expand', split=True)
sns.swarmplot(x="sharing_count_d140", y="mean_p_expand", data=foo, hue='did_expand', split=True)

# TODO
# Evaluate repertoire sharing at CDR3 seq level

# TODO graph proving lack of response in 1017 is not due to lack of clonotypes in naive repertoire
# perhaps add naive reps to heatmap

# TODO are shared clonotypes the same?
# Plot p value vs amount of sharing

# TODO PCA
# import sklearn.decomposition 
# pca = sklearn.decomposition.PCA(n_components=2)
# pca.fit(foo.transpose())
# bar = pca.transform(foo.transpose())

# TODO
# reorganisation required below

# In[156]:

def get_clonotype_mean_mut_freq(summary_df, var="mut_freq_vj"):
    '''For each clonotype present, get mean mutational freq per base pair
    '''
    return summary_df.groupby("clonotype")[var].mean()
#

def get_clonotype_mut_freqs(summary_df, var="mut_freq_vj"):
    '''For each clonotype present, get the mutational freqs of the clonotypes
    '''
    return summary_df.groupby("clonotype", "patient_code")[var]

def get_clonotype_mut_freqs(df, clonotype_field="clonotype", groups=None, prop=False, fillna=True):
    if groups is None:
        freqs = df.groupby([clonotype_field]).apply(len)
    else:
        freqs = df.groupby([clonotype_field] + groups).apply(len).unstack()
    if fillna:
        freqs = freqs.fillna(0)
    if prop:
        if len(freqs.shape) > 1:
            return freqs.apply(lambda x: x/sum(x)) 
        else:
            return freqs/sum(freqs) 
    else:
        return freqs

rep1_mut_rate = get_clonotype_mean_mut_freq(rep1)
rep2_mut_rate = get_clonotype_mean_mut_freq(rep2)

rep1_freq = get_clonotype_freq(rep1, "VJ-GENE", ["patient_code"])
rep2_freq = get_clonotype_freq(rep2, "VJ-GENE", ["patient_code"])
#
rep_freq_joined = rep1_freq.join(rep2_freq, how='outer', lsuffix="_rep1", rsuffix="_rep2").fillna(0)



# Filter by those detected in the naive rep of individuals?

#
# Diversity of a single repertoire
# Hill curves
#
# Calculate Hill curve
def calc_hill_curve(counts, max_q=15, step=0.1):
    '''
    '''
    counts = np.array([c for c in counts if not c == 0])
    rel_abund = counts/sum(counts)
    def calc_hill_num(rel_abund, q):
        # The Shannon index
        if q == 1.0:
            return math.exp(-sum(p*math.log(p) for p in rel_abund))
        # Other Hill numbers
        else:
            return math.pow(sum(math.pow(p, q) for p in rel_abund), (1/(1-q)))
    qs = np.arange(0, max_q, step)
    return (qs, [calc_hill_num(rel_abund, q) for q in qs])
#
# Plot a Hill curve for each repertoire, coding by timepoint and cell type
colors = dict(zip(sample_info_df["day"].unique(), ['red', 'green', 'blue']))
fig, axes = plt.subplots(ncols=3, sharey=True, sharex=True, figsize=(16, 10))
fig.suptitle("Repertoire Hill diversity curves")
for subplot_n, cell_type in enumerate(sample_info_df["cell_type"].unique()):
    axes[subplot_n].grid(True)
    axes[subplot_n].set_title(cell_type)
    for i in range(sample_info_df.shape[0]):
        if cell_type == sample_info_df.ix[i, "cell_type"]:
            sample_name = sample_info_df.ix[i, "Sample_name"]
            sample_num = sample_info_df.ix[i, "Tag Index"]
            day = sample_info_df.ix[i, "day"]
            repertoire = count_clonotype_cpm_normed[sample_name]
            axes[subplot_n].semilogy(calc_hill_curve(repertoire)[0], calc_hill_curve(repertoire)[1], color=colors[day], label=day)
            if subplot_n == 0: axes[subplot_n].set_ylabel("Diversity")
            if subplot_n == 1: axes[subplot_n].set_xlabel("q") 
            axes[subplot_n].legend()
fig.savefig('./.output/repertoire_hill_curves.svg')

#
# paired t test on arcsin transformed proportions of rep per clonotype
# 
foo = []
for i in range(len(rep_prop_arcsin_joined)):
    foo.append(stats.ttest_rel(rep_prop_arcsin_joined.iloc[i, :3], rep_prop_arcsin_joined.iloc[i, 3:]).pvalue)
sm.stats.multipletests(foo, method="fdr_bh")
#
# For clonotypes within the known mab set only
foo = []
for i in range(len(rep_prop_joined)):
    if (rep_prop_joined.index[i] in mAb_df["V-GENE"].values):
        foo.append(stats.ttest_rel(rep_prop_arcsin_joined.iloc[i, :3], rep_prop_arcsin_joined.iloc[i, 3:]).pvalue)

#
# boxplots of prop of repertoire per v gene
#
fig, axes = plt.subplots()
foo = pd.melt(rep1_prop.reset_index(), id_vars=['V-GENE'])
bar = pd.melt(rep2_prop.reset_index(), id_vars=['V-GENE'])
baz = pd.concat([foo, bar], keys=["foo", "bar"]).reset_index()
g = sns.boxplot(x="V-GENE", y="value", hue="level_0", data=baz, palette="PRGn")
#  g = sns.swarmplot(x="V-GENE", y="value", hue="level_0", data=baz, palette="PRGn", split=True)
for item in g.get_xticklabels():
    item.set_rotation(60)
g

#
# Hill diversity, 2 reps
#
fig, axes = plt.subplots()
rep1_clonotype_freq = get_clonotype_freq(rep1)
rep2_clonotype_freq = get_clonotype_freq(rep2)
rep1_hill_nums = calc_hill_curve(rep1_clonotype_freq.loc[:, 1017])
rep2_hill_nums = calc_hill_curve(rep2_clonotype_freq.loc[:, 1017])
axes.semilogy(rep1_hill_nums[0], rep1_hill_nums[1], label="rep1")
axes.semilogy(rep2_hill_nums[0], rep2_hill_nums[1], label="rep2")
axes.legend()

#
# Sharing of clonotypes: venn diagram
#
fig, axes = plt.subplots()
rep1_clonotypes = set(rep1_clonotype_freq.keys())
rep2_clonotypes = set(rep2_clonotype_freq.keys())
matplotlib_venn.venn2(subsets=(len(rep1_clonotypes - rep2_clonotypes), len(rep1_clonotypes & rep2_clonotypes), len(rep2_clonotypes - rep1_clonotypes)), set_labels=('rep1', 'rep2'))

#
# Find those clonotypes with increased Mutational frequencies
#
#
# Scatterplot of clonotype freq
fig, axes = plt.subplots()
df = pd.concat({"rep1": rep1_mut_rate, "rep2": rep2_mut_rate}, axis=1)
g = sns.jointplot("rep1", "rep2", data=df, xlim=(-0.0005, 0.0025), ylim=(-0.0005, 0.0025))
x0, x1 = g.ax_joint.get_xlim()
y0, y1 = g.ax_joint.get_ylim()
lims = [max(x0, y0), min(x1, y1)]
g.ax_joint.plot(lims, lims, ':k') 
#
# Violin plot of mean mut freq

df_complete = df[df["rep1"].notnull() & df["rep2"].notnull()]

#  shared rep: clonotypes that appear in all 3 patients at d140
shared_rep = rep2.groupby(["clonotype"]).apply(lambda x: len(x["patient_code"].unique()) == 3)

df_complete = df_complete.assign(shared=df_complete.index.isin(shared_rep[shared_rep].index))
df_complete_wide = pd.melt(df_complete, id_vars=["shared"])

fig, axes = plt.subplots()
sns.violinplot(x="variable", y="value", data=df_complete_wide, hue="shared", split=False)

#
# Mutational frequencies of shared clonotypes
#
shared = set(rep1["clonotype"]) & set(rep2["clonotype"])
ps = []
cs = []
#
with np.errstate(divide='raise'):
    for c in shared:
        rep1_mut_rates = np.arcsin(np.sqrt(rep1.loc[rep1["clonotype"] == c, "mut_freq_vj"]))
        rep2_mut_rates = np.arcsin(np.sqrt(rep2.loc[rep2["clonotype"] == c, "mut_freq_vj"]))
        try:
            ps.append(sm.stats.ttest_ind(rep1_mut_rates, rep2_mut_rates, alternative="smaller")[1])
        except FloatingPointError:
            ps.append(np.nan)
        cs.append(c)
bar = sm.stats.multipletests(ps, method="fdr_bh", alpha=0.05)
np.array(cs)[bar[0]]


#
# Regression of mutational frequency vs clonotype
#

#
# Isotype distributions of shared clonotypes
#

#
# Difference in proportion of clonotypes
#
#  import rpy2
#  from rpy2.robjects.packages import importr
#  from rpy2.robjects.vectors import StrVector
#  import rpy2.robjects as robjects
#  import rpy2.robjects.packages as rpackages
#  utils = rpackages.importr('utils')
#  utils.chooseCRANmirror(ind=1)
#  utils.install_packages(StrVector(['magrittr']))
#  IMGTStatClonotype = importr('IMGTStatClonotype')

# Comparison of expanded clonotypes to known sequences

# TODO
# unfinished below here

#
# Heatmap: sharing of unique clonotypes between individuals by time and cell type
#
foo = summary_df.query("cell_type == 'MBC'").groupby(["day", "patient_code"]).apply(lambda x: set(x["clonotype"]))
bar = np.reshape([len(x.intersection(y)) for x in foo for y in foo], newshape=(len(foo), len(foo)))
baz = (foo.keys().levels[0][foo.keys().labels[0]].astype("str")
    + " "
    + foo.keys().levels[1][foo.keys().labels[1]].astype("str")
    #  + " "
    #  + foo.keys().levels[2][foo.keys().labels[2]].astype("str")
)

fig, ax = plt.subplots()
ax.set_xticks(np.arange(22) + 0.5, minor=False)
ax.set_yticks(np.arange(22) + 0.5, minor=False)
ax.set_xticklabels(baz, minor=False)
ax.set_yticklabels(baz, minor=False)
ax.invert_yaxis()
ax.xaxis.tick_top()
#  plt.xticks(rotation=90)
heatmap = ax.pcolor(bar, cmap=plt.cm.Blues, alpha=0.8)

# In[120]:

adj = pd.DataFrame(bar)
adj.columns = baz
adj.to_csv("../team115_lustre/1_analyse_clonotypes/adj.csv", index=False)
adj

# In[84]:

# Plot clonotype abundances in each cell type, per patient, over time

for patient_code in ["1019", "2207"]:
# for patient_code in summary_df["patient_code"].unique():
    # Get clonotype count for each unique clonotype at each timepoint
    foo = summary_df.query("patient_code == {} & cell_type == 'MBC'".format(patient_code)).groupby(["day", "clonotype"]).apply(lambda x: x["norm_factor"].sum()).unstack().transpose().fillna(0)
    # Sort by freq at day 140
    foo.sort_values([140], ascending=True, inplace=True)
    days = np.array([0, 63, 140])
    fig, ax = plt.subplots()
    ax.stackplot(days, foo, baseline="sym", linewidth=0) 
    fig.set_size_inches(8, 8)
    plt.title(patient_code)
    plt.show()

# In[50]:

# Find V-genes that increased in abundance from day 0


def var_by_indiv_plots(summary_df, var="V-GENE", cell_type=None, day=None):
    # Filter out relevant cell types and days
    if cell_type and all((x in ["plasma", "MBC", "PBMCs"]) for x in cell_type): 
        summary_df = summary_df[summary_df["cell_type"].isin(cell_type)]
    if day and all((x in [0, 63, 140]) for x in day): 
        summary_df = summary_df[summary_df["day"].isin(day)]
        
    # Calculate normalised abundance of V genes by patient at day 0
    v_d0 = summary_df.groupby([var, "patient_code"]).sum()["norm_factor"].unstack().fillna(0)
    # Convert to proportion of repertoire
    v_d0_prop = v_d0.apply(lambda x: x/np.sum(x), axis=0)

    # Order V genes by total proportion of repertoire represented and plot
    v_d0_total_prop_idx = v_d0_prop.sum(axis=1).sort_values(ascending=False).index
    v_d0.loc[v_d0_total_prop_idx, :].plot(kind="bar", stacked=True, legend=None)
    # Normalise by total proportion of repertoire represented and replot
    # v_d0.loc[v_d0_total_prop_idx, :].apply(lambda x: x/np.sum(x), axis=1).plot(kind="bar", stacked=True, legend=None)

    # Order V genes by coefficent of variation and plot
    v_d0_cv_idx = v_d0.apply(lambda x: np.std(x)/np.mean(x), axis=1).sort_values().index
    v_d0.loc[v_d0_cv_idx, :].apply(lambda x: x/np.sum(x), axis=1).plot(kind="bar", stacked=True, legend=None)

#
# Day 0: distribution of clonotypes between individuals
#

var_by_indiv_plots(summary_df, var="V-GENE", day=[0])
# var_by_indiv_plots(summary_df, var="V-GENE", cell_type=["plasma"], day=[0]) 
# var_by_indiv_plots(summary_df, var="V-GENE", cell_type=["MBC"], day=[0]) 
# var_by_indiv_plots(summary_df, var="V-GENE", cell_type=["PBMCs"], day=[0]) 

# In[282]:

#
# Day 140: characterisation of long term memory
# In the memory compartment:
#     - increase in abundance from d0 to post prime: d28 > d0
#     - recalled at post boost (day 63): d63 > d0
#     - still present at d140: d140 > d0
#     - relatively high mutational freq (i.e. is the above group higher than expected)
#
# Compare against the naive repertoire:
#     IgD/IgM unmutated sequences taken from PBMC samples

def get_v_gene_freq(summary_df):
    '''For each V gene present, get normalised usage
    '''
    return summary_df.groupby("V-GENE")["norm_factor"].sum()

def get_isotype_freq(summary_df):
    '''For each isotype present, get normalised freq
    '''
    isotype_freqs = collections.defaultdict(float)
    for isotypes, norm_factors in zip(summary_df["isotypes"], summary_df["norm_factor"]):
        for isotype in isotypes:
            isotype_freqs[isotype] += norm_factors
    return pd.Series(isotype_freqs)

def get_clonotype_freq(summary_df):
    '''For each clonotype present, get normalised abundance
    '''
    return summary_df.groupby("clonotype")["norm_factor"].sum()

# Find V-genes that increased in abundance from day 0

# Get naive rep from day 0 
naive_d0 = get_naive_rep(summary_df).query("day == 0")

# Get plasmoblast rep from day 63
plasma_d63 = summary_df.query("day == 63 & cell_type == 'plasma'")

# Get memory rep from day 140
mbc_d140 = summary_df.query("day == 140 & cell_type == 'MBC'")

# In[324]:

naive_d0

# In[296]:

foo = pd.concat([get_v_gene_freq(naive_d0), 
           get_v_gene_freq(plasma_d63),
           get_v_gene_freq(mbc_d140)], axis=1)
foo.columns = ["naive d0", "plasma d63", "MBC d140"] 
print(foo.shape)
foo.apply(lambda x: x/np.sum(x)).plot(kind="bar")

# In[281]:

get_v_gene_freq(summary_df.query("day == 140 & cell_type == 'MBC'"))

# In[278]:

# How many clonotypes and V-genes overlap with those from known mAbs?
print(summary_df.query("day == 140 & cell_type == 'MBC'")["Ab.Name"].value_counts(dropna=False))
print(summary_df.query("day == 140 & cell_type == 'MBC'")["V-GENE_known_usage"].value_counts(dropna=False))

# In[297]:

summary_df

# In[304]:

print(summary_df["V-GENE"].unique().size)
print(summary_df["clonotype"].unique().size)

# In[322]:

bar = summary_df.groupby("sample_num").apply(lambda x: x.groupby("clonotype")["norm_factor"].sum())[1]
sorted(bar)

# In[323]:

baz = summary_df.groupby("sample_num").apply(lambda x: x["clonotype"].value_counts()[1])
sorted(baz)

