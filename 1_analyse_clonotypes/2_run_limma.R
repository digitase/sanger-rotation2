
# source("https://bioconductor.org/biocLite.R")
# biocLite("limma")
# biocLite("edgeR")

library(limma)
library(reshape2)
library(ggplot2)
library(dtplyr)
library(data.table)
library(edgeR)
library(lme4)
library(circlize)

# Read in counts matrix
counts.clonotype <- fread("../team115_lustre/1_analyse_clonotypes/count_clonotype.csv", header = T)
# counts.clonotype <- fread("../team115_lustre/1_analyse_clonotypes/count_vj.csv", header = T)
clonotype.names <- unlist(counts.clonotype[, 1, with=F])
counts.clonotype <- counts.clonotype[, -1, with=F]
rownames(counts.clonotype) <- clonotype.names

# Read in sample matrix
sampleInfo <- fread("../team115_lustre/1_analyse_clonotypes/sample_info.csv")
sampleInfo <- sampleInfo[order(sampleInfo$"Tag Index")]
sampleInfo <- sampleInfo[, patient_code := as.factor(patient_code)]
sampleInfo <- sampleInfo[, day := as.factor(day)]
sampleInfo <- sampleInfo[sampleInfo$"Tag Index" != "20", ]
sampleInfo[, group := paste(day, cell_type, sep=".")]
sampleInfo[, id := paste(Sample_name, patient_code, day, cell_type, sep=".")]
colnames(counts.clonotype) <- sampleInfo$id

# Create DGEList and normalise counts by library size
dge.clonotype <- DGEList(counts=counts.clonotype)
rownames(dge.clonotype) <- rownames(counts.clonotype)
dge.clonotype <- calcNormFactors(dge.clonotype)

# Filter out clonotypes with low counts
dge.clonotype.filtered <- dge.clonotype[rowSums(dge.clonotype$counts >= 2) >= 2,  , keep.lib.sizes=F]
plotMDS(dge.clonotype.filtered, labels=paste(sampleInfo$day, sampleInfo$cell_type, sep="_"), col=as.numeric(sampleInfo$patient_code))

# Since we need to make comparisons both within and between subjects, it is necessary to treat patient_code as a random effect.
design <- model.matrix(~ 0 + group, sampleInfo)
# When the library sizes are quite variable between samples, then the voom
# approach is theoretically more powerful than limma-trend.  In this approach,
# the voom transformation is applied to the normalized and filtered DGEList
# object.
v.clonotype <- voom(dge.clonotype.filtered, design, plot=T)
# Compensate for blocking term: estimate the correlation between measurements made on the same subject
cor.clonotype <- duplicateCorrelation(v.clonotype, design, block=sampleInfo$patient_code)

# Fit lm
fit.clonotype <- lmFit(v.clonotype, design, block=sampleInfo$patient_code, correlation=cor.clonotype$consensus.correlation)

# Setup and compute contrasts
cm.clonotype <- makeContrasts(
                              group63.MBC-group0.MBC,
                              group63.plasma-group0.MBC,
                              group140.MBC-group0.PBMCs,
                              group140.MBC-group0.MBC,
                              levels=design)
fit2.clonotype <- contrasts.fit(fit.clonotype, cm.clonotype)
# Perform shrinkage
fit2.clonotype <- eBayes(fit2.clonotype)
# Test for DE
summary(decideTests(fit2.clonotype, p.value=0.05))

topTable(fit2.clonotype, coef="group140.MBC - group0.MBC")

counts.clonotype[rownames(counts.clonotype) == "IGHV6-1.IGHJ5"]


#
# Paired t-test for proportions
#

For each

counts.clonotype[, sampleInfo[day == 0 & cell_type == "MBC", Sample_name], with=F]

counts.clonotype[, sampleInfo[day == 140 & cell_type == "MBC", Sample_name], with=F]

pairwise.prop.test(x, n, p.adjust.method = p.adjust.methods, ...)


cpm(dge.clonotype)



#
# Diversity of repertoire
#

summary(aov(gini_i ~ cell_type*day, sampleInfo))

summary(lmer(gini_i ~ cell_type+day + (1 | patient_code), sampleInfo))

(mixed(gini_i ~ cell_type*day + (1|patient_code), sampleInfo))

#
# 
#

# Generate adjacency matrix representing strength of relationships between samples
# Here, the strength of relationship is the number of unique clonotypes in the overlap

adj <- as.matrix(read.csv("../team115_lustre/1_analyse_clonotypes/adj.csv"))
rownames(adj) <- colnames(adj)
chordDiagram(adj, symmetric=F, self.link=1)

