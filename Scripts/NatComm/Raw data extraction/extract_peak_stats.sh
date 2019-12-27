#!/bin/bash


#PBS -t 1-100   # Parallel processing of 100 peak files. Thread stored in variable $PBS_ARRAYID

module load bedtools/2.17.0

# This script retrieves peak signals and scores from peak files at consensus peaks (i.e. peak stats).
# Also includes peak width and overlapping width for each peak.
# Input is first three columns (chr/start/end) of consensus and each BED file originally used to build the consensus
# Output is a peakStats file for each original sample

indir=NULL # Path to peak files
outdir=NULL
[[ ! -d "${outdir}" ]] && mkdir "${outdir}"
psdir=NULL # Path to store peak stats files
[[ ! -d "${psdir}" ]] && mkdir "${psdir}"
cd ${indir}

consensus="/path/to/consensus/consensus_narrowPeaks.bed";	
f=$(find . -name "*.narrowPeak" | head -$PBS_ARRAYID | tail -1)
cut -f1-3 ${consensus} \
	| perl -ne 'print if $.>1' | intersectBed -wao -a - -b $f | cut -f1-3,5-8,10-13 \
	| perl -lane 'print $F[0],"\t",$F[1],"\t",$F[2],"\t",$F[4]-$F[3],"\t",$F[5],"\t",$F[6],"\t",$F[7],"\t",$F[8],"\t",$F[9],"\t",$F[10]' > ${outdir}/$(basename $f .narrowPeak)_narrowPeaks.peakStats

cp ${outdir}/$(basename $f .narrowPeak)_narrowPeaks.peakStats ${psdir}
