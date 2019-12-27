#!/bin/bash

# for a given query (peaks of interest), finds motifs overlaps from a database of motifs, using Homer
# output placed in motif directory for respective query
# Homer install path hardcoded in call

for qid in `seq 1 8`; do
	cd /Volumes/emext/Dropbox/Jax/PBMC/DifferentialAnalysis/manuscripts/aging-02/figures/browser.shots/queries/query${qid}
	odir=motifs
	[[ ! -d  "${odir}" ]] && mkdir "${odir}"
	perl -lane 'print $F[0],"\t",$F[1],"\t",$F[2],"\t",$.,"_",$F[3],"_",$F[4]' ../query${qid}.qry > ${odir}/query${qid}.bed
	findMotifsGenome.pl ${odir}/query${qid}.bed hg19 ${odir} -bg ${odir}/query${qid}.bed -find ~/seqanal/homer/data/knownTFs/known.motifs > ${odir}/query${qid}_motifs.txt
done
