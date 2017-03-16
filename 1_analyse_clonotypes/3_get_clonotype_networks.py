# Generation of clonotype networks
import math 
import editdistance
import itertools
import pandas as pd
import igraph
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt
import random
import pygraphviz
from process_clonotype_data_helpers import *
plt.ion()
random.seed(42)

class Clonotype():
    '''
    Clonotypes have V, J, CDR3 seq, freq
    '''

    def __init__(self, name, v, j, cdr3, freq):
        self.name = name
        self.v = v
        self.j = j
        self.cdr3 = cdr3
        self.freq = freq

    def dist(self, c2):
        '''Get distance between two clonotypes

        Distance is inf if not the same VJ-gene.
        Distance is otherwise the edit distance between CDR seqs
        '''
        if self.v == c2.v and self.j == c2.j:
            return editdistance.eval(self.cdr3, c2.cdr3)
        else:
            return math.inf

    def __str__(self):
        return("{}_{}".format(self.name, self.freq))

    def __repr__(self):
        return("({}: {})".format(self.name, self.freq))

# Build fully connected graph of clonotypes,
# then drop all edges with distance > 1

summary_df = pd.read_csv("../team115_lustre/1_analyse_clonotypes/summary.csv")
clonotype_filter = ["V-GENE", "J-GENE", "CDR3_aa"]
def clonotype_format(x): return "{}.{}.CDR3_{}".format(*x)
summary_df["clonotype_aa"] = summary_df[clonotype_filter].apply(clonotype_format, axis=1)
rep = summary_df.query("day == 140 & cell_type == 'MBC' & patient_code == 1019")
rep_freqs = get_clonotype_freq(rep, "clonotype_aa")
rep_props = rep_freqs/sum(rep_freqs)

G = nx.Graph()
for c in rep_freqs.index:
    v, j, cdr3 = c.split(".")
    cdr3 = cdr3.split("_")[1]
    freq = rep_props[rep_freqs.index == c][0]
    G.add_node(Clonotype(c, v, j, cdr3, freq), freq=freq)
for c1, c2 in itertools.combinations(G.nodes(), 2):
    weight = c1.dist(c2)
    if weight == 1:
        G.add_edge(c1, c2, weight=weight)

# Remove singleton nodes
Gsub = G.copy()
for n in Gsub.nodes():
    if not len(Gsub.neighbors(n)):
        Gsub.remove_node(n)
Gsub = nx.nx_agraph.to_agraph(Gsub)

# pos = nx.spring_layout(Gsub,
        # k=4/math.sqrt(len(Gsub)),
        # iterations=60
        # )
# nx.draw_networkx_nodes(
        # Gsub, 
        # pos, 
        # node_size=[x[1]['freq']*2 for x in Gsub.nodes(data=True)],
        # alpha=0.6
        # )
# nx.draw_networkx_edges(
        # Gsub, 
        # pos)

Gsub.layout("fdp")
Gsub.draw("../team115_lustre/1_analyse_clonotypes/clonotype_networks.png")

for n in Gsub.nodes():
    n.name

