
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

# Set interactive plotting mode
plt.ion()

# Define clonotype field names
# - V gene, J gene, CDR3 length
clonotype_filter = ["V-GENE", "J-GENE", "CDR3_len"]
def clonotype_format(x): return "{}.{}.CDR3_len{}".format(*x)

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

# Biological unit: clonotype
count_clonotype =  summary_df.groupby(["sample_num", "clonotype"]).apply(len).unstack().transpose().fillna(0)
count_clonotype.to_csv("../team115_lustre/1_analyse_clonotypes/count_clonotype.csv")
#
# Biological unit: vj gene combo
count_vj =  summary_df.groupby(["sample_num", "VJ-GENE"]).apply(len).unstack().transpose().fillna(0)
count_vj.to_csv("../team115_lustre/1_analyse_clonotypes/count_vj.csv")
#
# Biological unit: v gene
count_v =  summary_df.groupby(["sample_num", "V-GENE"]).apply(len).unstack().transpose().fillna(0)
count_v.to_csv("../team115_lustre/1_analyse_clonotypes/count_v.csv")
#
# Write out sample information csv
sample_info_df.to_csv("../team115_lustre/1_analyse_clonotypes/sample_info.csv", index=False)

#
# Diversity of a single repertoire
# Hill curves
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

# Plot a Hill curve for each repertoire, coding by timepoint and cell type
colors = dict(zip(sample_info_df["day"].unique(), ['red', 'green', 'blue']))
fig, axes = plt.subplots(ncols=3, sharey=True, sharex=True, figsize=(16, 10))
fig.suptitle("Repertoire Hill diversity curves")
for subplot_n, cell_type in enumerate(sample_info_df["cell_type"].unique()):
    axes[subplot_n].grid(True)
    axes[subplot_n].set_title(cell_type)
    for i in range(sample_info_df.shape[0]):
        if cell_type == sample_info_df.ix[i, "cell_type"]:
            sample_num = sample_info_df.ix[i, "Tag Index"]
            day = sample_info_df.ix[i, "day"]
            repertoire = count_clonotype[sample_num]
            axes[subplot_n].semilogy(calc_hill_curve(repertoire)[0], calc_hill_curve(repertoire)[1], color=colors[day], label=day)
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

# In[124]:

# Generate counts matrices: biological unit by samples


# In[123]:

# Heatmap: sharing of unique clonotypes between individuals by time and cell type

# foo = summary_df.groupby(["patient_code", "cell_type", "day"]).apply(lambda x: set(x["clonotype"]))

foo = summary_df.query("cell_type == 'MBC'").groupby(["day", "patient_code"]).apply(lambda x: set(x["clonotype"]))

bar = np.reshape([len(x.intersection(y)) for x in foo for y in foo], newshape=(len(foo), len(foo)))

fig, ax = plt.subplots()
heatmap = ax.pcolor(bar, cmap=plt.cm.Blues, alpha=0.8)

fig.set_size_inches(11, 11)

baz = (foo.keys().levels[0][foo.keys().labels[0]].astype("str")
 + " "
 + foo.keys().levels[1][foo.keys().labels[1]].astype("str")
#  + " "
#  + foo.keys().levels[2][foo.keys().labels[2]].astype("str")
)

ax.set_xticks(np.arange(22) + 0.5, minor=False)
ax.set_yticks(np.arange(22) + 0.5, minor=False)

ax.set_xticklabels(baz, minor=False)
ax.set_yticklabels(baz, minor=False)

# want a more natural, table-like display
ax.invert_yaxis()
ax.xaxis.tick_top()

plt.xticks(rotation=90)

plt.show()


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

def get_naive_rep(summary_df):
    '''Filter for the naive repertoire
    '''
    return summary_df.loc[(summary_df["V-REGION identity %"] == 100) 
                          & (summary_df["digest"].map(lambda x: "IGHM" in x or "IGHG" in x))
                          & (summary_df["cell_type"] == "PBMCs")]

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

def get_naive_rep(summary_df):
    '''Filter for the naive repertoire
    '''
    return summary_df.loc[(summary_df["V-REGION identity %"] == 100) 
                          & (summary_df["digest"].map(lambda x: "IGHM" in x or "IGHG" in x))
                          & (summary_df["cell_type"] == "PBMCs")]

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

