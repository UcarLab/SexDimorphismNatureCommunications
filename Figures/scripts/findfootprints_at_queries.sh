#!/bin/bash

# for a given query (peaks of interest), finds overlap with pre-calculated footprints, by sample
# output placed in motif directory for respective query

fpdir="/Volumes/emback/temp/aging2018_footprints" # location of footprints
cd /Volumes/emext/Dropbox/Jax/PBMC/DifferentialAnalysis/manuscripts/aging-02/figures/browser.shots/queries # location of queries

for qid in `seq 1 8`; do
	bed=query${qid}.qry
	odir=query${qid}/motifs
	[[ ! -d  "$odir" ]] && mkdir "$odir"

	for f in $(find $fpdir -name "*_pooledPIQcalls.bed"); do 
		samp=$(basename $f | sed s/\_.*//g);
		bedtools intersect -wo -a $f -b $bed > ${odir}/${samp}_query${qid}_footprints.bed
	done
done
