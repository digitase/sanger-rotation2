#!/usr/local/bin/bash

# Converts datasets from imgthighvquest format to VDJtools format. 

JAVA="/nfs/users/nfs_b/bb9/packages/jre1.8.0_121/bin/java"
VDJTOOLS="/nfs/users/nfs_b/bb9/packages/vdjtools-1.1.4/vdjtools-1.1.4.jar"

OUTDIR="../team115_lustre/1_vdjtools/converted/"

mkdir -p "$OUTDIR"
shopt -s globstar

for i in "../team115_lustre/data/2017-03-01_high_vquest_results/"**/3_Nt-sequences.txt; do
    prefix="$OUTDIR/$(basename $(dirname $i)).vdjtools"
    echo $prefix
    "$JAVA" -Xmx2G -jar "$VDJTOOLS" Convert -S imgthighvquest "$i" "$prefix"
done
    
