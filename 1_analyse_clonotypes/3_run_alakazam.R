#
# Direct sourcing of alakazam modules to generate hill curves
#
# Necessary due to failure of installation
# cc1plus: error: unrecognised command line option ‘-std=c++11’
#
library(lazyeval)
library(ggplot2)
# library(dtplyr)
library(dplyr)
library(grid)
library(data.table)

source("~/packages/alakazam/R/Diversity.R")
source("~/packages/alakazam/R/Core.R")
source("~/packages/alakazam/R/Classes.R")

# Convert summary format from IMGT to alakazam format
summary2Alakazam <- function(df) {
    data.frame(
        SEQUENCE_ID=df$"Sequence ID",
        SEQUENCE_IMGT=NA,
        V_CALL=df$"V-GENE and allele",
        V_CALL_GENOTYPED=df$"V-GENE",
        D_CALL=df$"D-GENE and allele",
        J_CALL=df$"J-GENE and allele",
        JUNCTION=df$"AA JUNCTION",
        JUNCTION_LENGTH=NA,
        NP1_LENGTH=NA,
        NP2_LENGTH=NA,
        SAMPLE=df$"sample_num",
        DAY=df$"day",
        ISOTYPE=df$"isotypes",
        DUPCOUNT=1,
        CLONE=df$"VJ-GENE"
    )
}

# Read in repertoires
summaryDf <- fread("../team115_lustre/1_analyse_clonotypes/summary.csv")

pdf("../team115_lustre/1_analyse_clonotypes/hill_curves_alakazam.pdf", width=3, height=3)
for (p in unique(summaryDf$"patient_code")) {

    rep1 <- summaryDf[cell_type == "MBC" & (day == 0 | day == 140) & patient_code == p]
    rep1Db <- summary2Alakazam(rep1)
    sample_div <- rarefyDiversity(rep1Db, "DAY", min_q=0, max_q=32, step_q=0.05, ci=0.95, min_n=5, nboot=200)
    # Plot a log-log (log_q=TRUE, log_d=TRUE) plot of sample diversity
    # sample_main <- paste0("Patient code: ", p, "\nSample diversity (no. clones with > 30 reads = ", sample_div@n, ")")
    # sample_main <- paste0("Patient code: ", p)
    # No sample names for poster
    sample_main <- ""
    sample_colors <- c("0"="seagreen", "140"="steelblue")
    curve = plotDiversityCurve(sample_div, colors=sample_colors, main_title=sample_main, legend_title="Day", log_q=F, log_d=T, ylim=c(NA, 2**7.5),
                       silent = T,
                       legend.position = c(0.80, 0.80)
                       )

    # Adjust the theming
    curve = curve + 
        labs(x="Diversity order (q)", y=expression('Diversity ('^q * D * ')')) +
        theme(
            plot.title = element_text(size=12),
            axis.title = element_text(size=12),
            axis.text.x = element_text(size=11),
            axis.text.y = element_text(size=11),
            legend.key.height = unit(0.3, "in"),
            legend.background = element_rect(fill="transparent"),
            legend.title = element_text(size=12),
            legend.text = element_text(size=12)
        )

    # Remove legend except for last plot, for poster
    # if (p != 2207) {
        # curve = curve + theme(legend.position="none")
    # }

    # Try to reduce margins
    par(mar=c(0,0,0,0), oma=c(0,0,0,0)) 
    plot(curve, xaxs = "i", yaxs = "i")

}
dev.off()

# Vignette
if(F) {
    load("~/packages/alakazam/data/ExampleDb.rda")
    str(ExampleDb)
    # df, required columns in ExampleDb for rarefyDiversity function
    #     <group>: grouping factor for which to draw multiple curves e.g. timepoint, isotype
    #     <clone>: clonotype colname
    #     <copy>: name of the data column containing copy numbers for each sequence. If copy=NULL
                # (the default), then clone abundance is determined by the number of sequences.
                # If a copy column is specified, then clone abundances is determined by the sum
                # of copy numbers within each clonal group
    #
    # Compare diversity curve across values in the "SAMPLE" column
    # q ranges from 0 (min_q=0) to 32 (max_q=32) in 0.05 incriments (step_q=0.05)
    # A 95% confidence interval will be calculated (ci=0.95)
    # 2000 resampling realizations are performed (nboot=200)
    sample_div <- rarefyDiversity(ExampleDb, "SAMPLE", min_q=0, max_q=32, step_q=0.05, ci=0.95, nboot=200)
    #
    # Plot a log-log (log_q=TRUE, log_d=TRUE) plot of sample diversity
    # Indicate number of sequences resampled from each group in the title
    sample_main <- paste0("Sample diversity (n=", sample_div@n, ")")
    sample_colors <- c("-1h"="seagreen", "+7d"="steelblue")
    plotDiversityCurve(sample_div, colors=sample_colors, main_title=sample_main, legend_title="Sample", log_q=F, log_d=TRUE)
}

