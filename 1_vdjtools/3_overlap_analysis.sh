#!/usr/local/bin/bash

JAVA="/nfs/users/nfs_b/bb9/packages/jre1.8.0_121/bin/java"
VDJTOOLS="/nfs/users/nfs_b/bb9/packages/vdjtools-1.1.4/vdjtools-1.1.4.jar"
PATH=/software/R-3.3.0/bin/:$PATH
OUTDIR="/nfs/users/nfs_b/bb9/workspace/rotation2/1_vdjtools/.output"
Rscript --version

function overlap_pair {
    f1="$1"
    f2="$2"
    prefix="$3"
    "$JAVA" -Xmx2G -jar "$VDJTOOLS" OverlapPair \
        --intersect-type aaVJ \
        --top 20 \
        --plot --plot-area-v2 --plot-type pdf \
        "$f1" "$f2" "$prefix"
}

overlap_pair \
        "$OUTDIR/converted/isotyper_Fully_reduced_LEA_S1.vdjtools.3_Nt-sequences.txt" \
        "$OUTDIR/converted/isotyper_Fully_reduced_LEA_S3.vdjtools.3_Nt-sequences.txt" \
        "$OUTDIR/1017"

overlap_pair \
        "$OUTDIR/converted/isotyper_Fully_reduced_LEA_S11.vdjtools.3_Nt-sequences.txt" \
        "$OUTDIR/converted/isotyper_Fully_reduced_LEA_S9.vdjtools.3_Nt-sequences.txt" \
        "$OUTDIR/1019"

overlap_pair \
        "$OUTDIR/converted/isotyper_Fully_reduced_LEA_S15.vdjtools.3_Nt-sequences.txt" \
        "$OUTDIR/converted/isotyper_Fully_reduced_LEA_S7.vdjtools.3_Nt-sequences.txt" \
        "$OUTDIR/2207"

#Performf an all-versus-all pairwise overlap for a list of samples and computes a set of repertoire similarity measures. 
#
#Pairwise overlap circos plot. Count, frequency and diversity panels correspond
#to the read count, frequency (both non-symmetric) and the total number of
#clonotypes that are shared between samples. Pairwise overlaps are stacked, i.e.
#segment arc length is not equal to sample size.  
"$JAVA" -Xmx2G -jar "$VDJTOOLS" CalcPairwiseDistances \
    --intersect-type aaVJ \
    --plot \
    -m "$OUTDIR/converted/metadata.txt" \
    "$OUTDIR/pairwise"

#This routine provides additional cluster analysis (hierarchical clustering),
#multi-dimensional scaling (MDS) and plotting for CalcPairwiseDistances output.
"$JAVA" -Xmx2G -jar "$VDJTOOLS" ClusterSamples \
    --intersect-type aaVJ \
    --plot \
    "$OUTDIR/pairwise"

#This routine performs an all-vs-all intersection between an ordered list of samples for clonotype tracking purposes. 
"$JAVA" -Xmx2G -jar "$VDJTOOLS" TrackClonotypes \
    --intersect-type aaVJ \
    --plot \
    -m "$OUTDIR/converted/metadata.txt" \
    "$OUTDIR/track"

