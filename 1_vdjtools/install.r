#args <- c("/nfs/users/nfs_b/bb9/packages/vdjtools-1.1.4/Rpackages/", "reshape2", "FField", "reshape", "gplots", "gridExtra", "circlize", "ggplot2", "grid", "VennDiagram", "RColorBrewer", "ape", "MASS", "plotrix", "scales")
.libPaths("/nfs/users/nfs_b/bb9/packages/vdjtools-1.1.4/Rpackages/")
# This simple script automates R dependency installation for VDJtools
args    <- commandArgs(TRUE)
libPath <- args[1]

# iterate over dependencies
for (i in 2:length(args)) {
   # install packages to local library, use 0-Cloud mirror
   # note that local library is set as default search path path by RUtil
   install.packages(args[i], lib=libPath, repos="http://cran.rstudio.com/")
}
