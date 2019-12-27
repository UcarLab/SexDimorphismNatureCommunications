#!/bin/bash

module load bedtools/2.17.0

outdir=NULL
indir=NULL # Path to individual sample MACS2 peak files

cd ${indir};

bedtools multiinter -i $(find . -name "*narrowPeaks.bed" -maxdepth 1) > ${outdir}/narrowPeaks_multiinter.txt;
bedtools merge -i ${outdir}/narrowPeaks_multiinter.txt > ${outdir}/consensus_narrowPeaks.tmp;
bedtools annotate -counts -i ${outdir}/consensus_narrowPeaks.tmp -files $(find . -maxdepth 1 -name "*narrowPeaks.bed")	| sortBed | cut -f4- | sed 's/[1-9]/1/g' | sed 's/	/+/g' | sed 's/+$//g' | bc | paste ${outdir}/consensus_narrowPeaks.tmp - > ${outdir}/consensus_narrowPeaks.bed;

rm ${outdir}/consensus_narrowPeaks.tmp;
