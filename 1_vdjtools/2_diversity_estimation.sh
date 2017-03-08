#!/usr/local/bin/bash

JAVA="/nfs/users/nfs_b/bb9/packages/jre1.8.0_121/bin/java"
VDJTOOLS="/nfs/users/nfs_b/bb9/packages/vdjtools-1.1.4/vdjtools-1.1.4.jar"
PATH=/software/R-3.3.0/bin/:$PATH

OUTDIR="/nfs/users/nfs_b/bb9/workspace/rotation2/1_vdjtools/.output"

Rscript --version

# Rarefaction curves
"$JAVA" -Xmx2G -jar "$VDJTOOLS" RarefactionPlot \
    "$OUTDIR/converted/"*3_Nt-sequences.txt \
    "$OUTDIR"

