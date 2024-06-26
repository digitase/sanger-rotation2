#
# Attempt to use DE software to look for increased clonotype props
# Ended up just using the MDS plots.
#
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
# counts.clonotype <- fread("../team115_lustre/1_analyse_clonotypes/count_clonotype.csv", header = T)
counts.clonotype <- fread("../team115_lustre/1_analyse_clonotypes/count_vj.csv", header = T)
# counts.clonotype <- fread("../team115_lustre/1_analyse_clonotypes/count_v.csv", header = T)
clonotype.names <- unlist(counts.clonotype[, 1, with=F])
counts.clonotype <- counts.clonotype[, -1, with=F]
rownames(counts.clonotype) <- clonotype.names
#
# Read in sample matrix
sampleInfo <- fread("../team115_lustre/1_analyse_clonotypes/sample_info.csv")
sampleInfo <- sampleInfo[order(sampleInfo$"Tag Index")]
sampleInfo <- sampleInfo[, patient_code := as.factor(patient_code)]
sampleInfo <- sampleInfo[, day := as.factor(day)]
sampleInfo[, group := paste(day, cell_type, sep=".")]
sampleInfo[, id := paste(Sample_name, patient_code, day, cell_type, sep=".")]
colnames(counts.clonotype) <- sampleInfo$Sample_name
#
# Create DGEList and normalise counts by library size
dge.clonotype <- DGEList(counts=counts.clonotype)
rownames(dge.clonotype) <- rownames(counts.clonotype)
dge.clonotype <- calcNormFactors(dge.clonotype)
# Write out normalised counts
write.csv(cpm(dge.clonotype, log=F, normalised.lib.sizes=T), "../team115_lustre/1_analyse_clonotypes/count_clonotype.cpm_normed.csv")
colnames(counts.clonotype) <- sampleInfo$id
colnames(dge.clonotype) <- sampleInfo$id
#
# Filter out clonotypes with low counts
dge.clonotype.filtered <- dge.clonotype[rowSums(dge.clonotype$counts >= 2) >= 0,  , keep.lib.sizes=F]
# dge.clonotype.filtered <- dge.clonotype[colnames(dge.clonotype.filtered) != "LEA_S20.1019.63.PBMCs", ]

# Filter for samples of interest
# toKeep <- T
toKeep <- colnames(dge.clonotype.filtered) %in% unlist(sampleInfo[sampleInfo[, cell_type == "MBC" & (day == 0 | day == 140)], "id", with=F])
dge.clonotype.filtered <- dge.clonotype.filtered[, toKeep]

pdf('../team115_lustre/1_analyse_clonotypes/sample_mds_edgeR.pdf', width=6.5, height=6.0)
plotMDS(dge.clonotype.filtered, 
        # labels=paste(sampleInfo$day, sampleInfo$cell_type, sep="_")[toKeep], 
        labels=paste("day_", sampleInfo$day, sep="")[toKeep], 
        col=c("darkgrey", "blue", "orange")[as.numeric(sampleInfo$patient_code)[toKeep]],
        # col=as.numeric(sampleInfo$patient_code)[toKeep],
        gene.selection="common"
)
legend("top", legend=c("Patient 1", "Patient 2", "Patient 3"), fill=c("darkgrey", "blue", "orange"))
# legend("bottomright", legend=sampleInfo$patient_code[toKeep], fill=as.numeric(sampleInfo$patient_code)[toKeep])
#
arrows(x0=c(-0.20, 0.4), y0=c(-0.45, 0.2), x1=c(-1.1, 0.65), y1=c(0.3, 0.55), col=c("orange", "blue"))
dev.off()

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
#
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

topTable(fit2.clonotype, coef="group140.MBC - group0.MBC", n=20)

topTable(fit2.clonotype, coef="group140.MBC - group0.PBMCs", n=20)

getCounts <- function(clonotype) {
    print( (dge.clonotype$counts)[rownames(dge.clonotype) == clonotype][grepl(".0.PBMCs", colnames(dge.clonotype), fixed=T)] )
    print( cpm(dge.clonotype)[rownames(dge.clonotype) == clonotype][grepl(".0.PBMCs", colnames(dge.clonotype), fixed=T)] )
    print( (dge.clonotype$counts)[rownames(dge.clonotype) == clonotype][grepl(".140.PBMCs", colnames(dge.clonotype), fixed=T)] )
    print( cpm(dge.clonotype)[rownames(dge.clonotype) == clonotype][grepl(".140.PBMCs", colnames(dge.clonotype), fixed=T)] )
}

getCounts("IGHV3-7.IGHJ4.CDR3_len23")
getCounts("IGHV4-38-2.IGHJ5.CDR3_len8")
getCounts("IGHV4-28.IGHJ5.CDR3_len11")

# TODO
# unfinished below
stopifnot(F)

#
# Paired t-test for proportions
#

survivors <- matrix(c(159,1286,142,2029), ncol=2)

prop.test(survivors)
chisq.test(survivors)

#
# Barnard test - difference in proportions
# multinomial for nonfixed margins, too slow.
#

library(Exact)

exact.test(survivors)
exact.test(c(159,1286,142,2029), ncol=2)

prop.test(c(159, 12), c(1286, 2029))

foo=matrix(c(159,142,1286,2029), ncol=2)
exact.test(foo, model="multinomial")

#
# Diversity of repertoire
#

summary(aov(gini_i ~ cell_type*day, sampleInfo))

summary(lmer(gini_i ~ cell_type+day + (1 | patient_code), sampleInfo))

(mixed(gini_i ~ cell_type*day + (1|patient_code), sampleInfo))

# Generate adjacency matrix representing strength of relationships between samples
# Here, the strength of relationship is the number of unique clonotypes in the overlap

adj <- as.matrix(read.csv("../team115_lustre/1_analyse_clonotypes/adj.csv"))
rownames(adj) <- colnames(adj)
chordDiagram(adj, symmetric=F, self.link=1)

