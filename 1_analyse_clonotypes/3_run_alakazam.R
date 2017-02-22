#
# Direct sourcing of alakazam modules
#
# Necessary due to failure of installation
# cc1plus: error: unrecognised command line option ‘-std=c++11’
#
library(dplyr)
library(lazyeval)
library(ggplot2)

source("~/packages/alakazam/R/Diversity.R")
source("~/packages/alakazam/R/Core.R")
source("~/packages/alakazam/R/Classes.R")

# Vignette
if(F) {
    load("~/packages/alakazam/data/ExampleDb.rda")
    head(ExampleDb)
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



