
import pandas as pd
import numpy as np
from matplotlib import colors as mcolors
import matplotlib.pyplot as plt
import glob
import re
import math
from scipy import stats
from collections import Counter, defaultdict
from functools import reduce
import seaborn as sns
import statsmodels.api as sm
from matplotlib_venn import venn2

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
summary_df["CDR3_len"] = summary_df["CDR3-IMGT length"]
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
mAb_df["CDR3_len"] = mAb_df["CDR3 amino acid seq"].map(len)
mAb_df["clonotype"] = mAb_df[clonotype_filter].apply(clonotype_format, axis=1)
# Get clonotypes for existing data
summary_df["clonotype"] = summary_df[clonotype_filter].apply(clonotype_format, axis=1)
summary_df["VJ-GENE"] = summary_df["V-GENE"] + "." + summary_df["J-GENE"]
#
# Read in sample information
sample_info_df = pd.read_excel("/nfs/users/nfs_b/bb9/workspace/rotation2/team115_lustre/data/existing_bcr_data/MalariaSamplesBen/Malaria_Samples_SeqInfo.xlsx")
sample_info_df["patient_code"] = sample_info_df["patient_code"].astype("category")
#
# Calculate normalisation factor per sample.
# If a sample has an above average number of reads, 
# the contribution to normalised clonotype abundance by that sample should be lower.
sample_info_df["norm_factor"] = np.mean(sample_info_df["Reads, count"]) / sample_info_df["Reads, count"]
#
# Merge in patient number and cell type
summary_df = pd.merge(summary_df, sample_info_df[["Tag Index", "patient_code", "cell_type", "norm_factor"]], 
         how="left", left_on="sample_num", right_on="Tag Index")
#
# Merge correspondance to known mAbs
summary_df = pd.merge(summary_df, mAb_df[["Ab.Name", "clonotype"]], 
         how="left", on="clonotype")
#
# Check if V-gene is in the set of V-genes used in known mAbs
summary_df["V-GENE_known_usage"] = np.where(summary_df["V-GENE"].isin(mAb_df["V-GENE"]), summary_df["V-GENE"], "NaN")
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
summary_df["mut_freq_per_bp_v"] = (
    (summary_df["V-REGION_len"] - summary_df["V-REGION_unmut_len"])
    /summary_df["V-REGION_len"]
    /summary_df["V-REGION_len"]
)
summary_df["mut_freq_per_bp_j"] = (
    (summary_df["J-REGION_len"] - summary_df["J-REGION_unmut_len"])
    /summary_df["J-REGION_len"]
    /summary_df["J-REGION_len"]
)
summary_df["mut_freq_per_bp_vj"] = (
    (summary_df["VJ-REGION_len"] - summary_df["V-REGION_unmut_len"] - summary_df["J-REGION_unmut_len"])
    / summary_df["VJ-REGION_len"]
    / summary_df["VJ-REGION_len"]
)
#
# Biological unit: clonotype
count_clonotype = summary_df.groupby(["sample_num", "clonotype"]).apply(len).unstack().transpose().fillna(0)
count_clonotype.to_csv("../team115_lustre/1_analyse_clonotypes/count_clonotype.csv")
#
# Biological unit: vj gene combo
count_vj = summary_df.groupby(["sample_num", "VJ-GENE"]).apply(len).unstack().transpose().fillna(0)
count_vj.to_csv("../team115_lustre/1_analyse_clonotypes/count_vj.csv")
#
# Biological unit: v gene
count_v = summary_df.groupby(["sample_num", "V-GENE"]).apply(len).unstack().transpose().fillna(0)
count_v.to_csv("../team115_lustre/1_analyse_clonotypes/count_v.csv")
#
# Write out sample information csv
sample_info_df.to_csv("../team115_lustre/1_analyse_clonotypes/sample_info.csv", index=False)
#
# Read back in normalised counts
count_clonotype_cpm_normed = pd.read_csv("./.output/count_clonotype.cpm_normed.csv", header=0, index_col=0)

#
# Diversity of a single repertoire
# Hill curves
#
#
# Calculate Hill curve
def calc_hill_curve(counts, max_q=15, step=0.1):
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
# # Questions to ask
# 
# How is immunity to the vaccine candidate achieved? Due to antibodies... which we detect via:
# 
# - V-gene usage
# - clonotype abundance
# - mutational frequency
# 
# Things:
# 
# - V(D)J gene frequencies, Ig isotype usage, and BCR clone size, stratified by timepoint and patient
# - Look for expansion of clonotypes in day 63 and 140 vs. day 0
# - Compare expanded clonotypes to known Ab sequences
# - Look for differences between patients at day 0
# - Calculate mutational frequencies of clonotypes
# - Area/bar/parallel coords chart showing expansion of clonotypes (known V genes?) over time
# - 1 long stream plot showing abundances of clones over time in the 3 individuals
# 
# ## Comparison across individuals at pre-prime (day 0)
# 
# RH5 and Duffy samples
# 
# - Distribution of clonotypes between individuals
# - How diverse is the repertoire? Gini index?
# - Run comparison across cell types
# 
# ## Clonotypes at early post prime (day 28)
# 
# Look for:
# 
# - clonotype freq
# - isotype distribution
# - mutational freq (1-% identity) vs. IMGT V gene reference
#     - linear regression of BCR mutational frequency against clonotype V gene status
#     - mutational freq stratified by isotype
#     
# In the compartments:
# 
# - For memory cells vs day 0
# - For plasmablasts vs naive repertoire (IgD/IgM unmutated sequences taken from PBMC samples at Day 0)
# 
# - Are expanded/mutated clonotypes the same as the
#     - known anti-RH5 Ab sequences
#     - ones observed in AMA1 and Duffy trials
#     - ones observed in influenza and other infection challenge
# 
# ## Clonotypes at early post boost (day 63)
# 
# Look for:
# 
# - clonotype freq
# - isotype distribution
# - mutational freq (1-% identity) vs. IMGT V gene reference
#     - linear regression of BCR mutational frequency against clonotype V gene status
#     - mutational freq stratified by isotype
#     
# In the compartments:
# 
# - For memory cells vs day 0 and day 28
# - For plasmablasts vs plasmoblasts at day 0 and day 28
#     - activation of same clonotypes during prime and boost?
# - For plasmablasts vs naive repertoire at day 0 and day 28
# - For plasmablasts vs memory repertoire at day 0 and day 28
#     - recall of memory cells generated at the prime?
# 
# - clonotype overlap between 
# 
# ## Clonotypes at long term (day 140)
# 
# Look for:
# 
# - clonotype freq
# 
# In the compartments:
# 
# - For plasmablasts vs memory repertoire at day 63
#     - how much of the boost response is detectable in long term memory?
# 
# ## Identification of malaria specific mAbs
# 
# - Overlap with known anti RH5 mAbs (isolated from day 63), or clonotypes appearing in other infections/vaccinations in the literature
# - Can we expect cross-reactivity in clonotypes that would affect the inferences based on the BCR repertoire data? 

#
# Two timepoint analysis
#
def get_naive_rep(summary_df, negate=False):
    '''Filter for the naive repertoire
    '''
    mask = (summary_df["V-REGION identity %"] == 100
            & summary_df["digest"].map(lambda x: "IGHM" in x or "IGHG" in x)
            & summary_df["cell_type"] == "PBMCs")
    if negate:
        return summary_df.loc[np.logical_not(mask)]
    else:
        return summary_df[mask]
#
#  def get_clonotype_freq(summary_df):
    #  '''For each clonotype present, get normalised abundance
    #  '''
    #  return summary_df.groupby("clonotype")["norm_factor"].sum()
#
#
#  def get_clonotype_prop(summary_df):
    #  '''For each clonotype present, get proportion of total repertoire
    #  '''
    #  return summary_df.groupby("clonotype")["norm_factor"].sum()
#
def get_clonotype_freq(df, clonotype_field="clonotype", prop=False):
    freqs = df.groupby([clonotype_field, "patient_code"]).apply(len).unstack().fillna(0)
    if prop:
        return freqs.apply(lambda x: x/sum(x)) 
    else:
        return freqs

rep1 = summary_df.query("day == 0 & cell_type == 'MBC'")
rep2 = summary_df.query("day == 140 & cell_type == 'MBC'")
#
rep1_prop = get_clonotype_freq(rep1, "V-GENE", prop=True)
rep2_prop = get_clonotype_freq(rep2, "V-GENE", prop=True)
#
rep1_prop_arcsin = rep1_prop.apply(lambda x: 2*np.arcsin(np.sqrt(x)))
rep2_prop_arcsin = rep2_prop.apply(lambda x: 2*np.arcsin(np.sqrt(x)))
#
rep_prop_joined = rep1_prop.join(rep2_prop, how='outer', lsuffix="_1", rsuffix="_2").fillna(0)
rep_prop_arcsin_joined = rep1_prop_arcsin.join(rep2_prop_arcsin, how='outer', lsuffix="_1", rsuffix="_2").fillna(0)
#

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
venn2(subsets=(len(rep1_clonotypes - rep2_clonotypes), len(rep1_clonotypes & rep2_clonotypes), len(rep2_clonotypes - rep1_clonotypes)), set_labels=('rep1', 'rep2'))

#
# tests for difference in proportion of the repertoire between two reps
#
def prop_test_ztest(rep1_freqs, rep2_freqs):
    assert(len(rep1_freqs) == len(rep2_freqs))
    ps = []
    for i in range(len(rep1_freqs)):
        if rep1_freqs[i] or rep2_freqs[i]:
            ps.append(sm.stats.proportions_ztest([rep1_freqs[i], rep2_freqs[i]], [sum(rep1_freqs), sum(rep2_freqs)])[1])
        # Situation where both counts are 0
        else:
            ps.append(1.0)
    return ps
#
def prop_test_fisher(rep1_freqs, rep2_freqs):
    assert(len(rep1_freqs) == len(rep2_freqs))
    ps = []
    for i in range(len(rep1_freqs)):
        cont_table = np.array([[rep1_freqs[i], sum(rep1_freqs) - rep1_freqs[i]], 
                               [rep2_freqs[i], sum(rep2_freqs) - rep2_freqs[i]]])
        ps.append(stats.fisher_exact(cont_table)[1])
    return ps

rep1_freq = get_clonotype_freq(rep1, "clonotype", prop=False)
rep2_freq = get_clonotype_freq(rep2, "clonotype", prop=False)
rep_freq_joined = rep1_freq.join(rep2_freq, how='outer', lsuffix="_1", rsuffix="_2").fillna(0)
#
foo = prop_test_ztest(rep_freq_joined["1017_1"], rep_freq_joined["1017_2"])
bar = sm.stats.multipletests(foo, method="fdr_bh", alpha=0.1)
print(rep_freq_joined.index[bar[0]][rep_freq_joined.index[bar[0]].isin(mAb_df["clonotype"])])
#
foo = prop_test_ztest(rep_freq_joined["1019_1"], rep_freq_joined["1019_2"])
bar = sm.stats.multipletests(foo, method="fdr_bh", alpha=0.1)
print(rep_freq_joined.index[bar[0]][rep_freq_joined.index[bar[0]].isin(mAb_df["clonotype"])])
#
foo = prop_test_ztest(rep_freq_joined["2207_1"], rep_freq_joined["2207_2"])
bar = sm.stats.multipletests(foo, method="fdr_bh", alpha=0.1)
print(rep_freq_joined.index[bar[0]][rep_freq_joined.index[bar[0]].isin(mAb_df["clonotype"])])

#
# Find those clonotypes with increased Mutational frequencies
#
def get_clonotype_mean_mut_freq(summary_df, var="mut_freq_per_bp_vj"):
    '''For each clonotype present, get mean mutational freq per base pair
    '''
    return summary_df.groupby("clonotype")[var].mean()
#
def get_clonotype_mut_freqs(summary_df, var="mut_freq_per_bp_vj"):
    '''For each clonotype present, get mean mutational freq per base pair
    '''
    return summary_df.groupby("clonotype", "patient_code")[var]

rep1_mut_rate = get_clonotype_mean_mut_freq(rep1)
rep2_mut_rate = get_clonotype_mean_mut_freq(rep2)
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
sns.violinplot(x="variable", y="value", data=df_complete_wide, hue="shared")

#
# Mutational frequencies of shared clonotypes
#
shared = set(rep1["clonotype"]) & set(rep2["clonotype"])
ps = []
cs = []
#
with np.errstate(divide='raise'):
    for c in shared:
        rep1_mut_rates = np.arcsin(np.sqrt(rep1.loc[rep1["clonotype"] == c, "mut_freq_per_bp_vj"]))
        rep2_mut_rates = np.arcsin(np.sqrt(rep2.loc[rep2["clonotype"] == c, "mut_freq_per_bp_vj"]))
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

# Get naive rep from day 0 
# naive_d0 = get_naive_rep(summary_df).query("day == 0")
get_naive_rep(summary_df)

# In[156]:

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
    isotype_freqs = defaultdict(float)
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

