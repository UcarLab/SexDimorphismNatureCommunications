#!/bin/bash

#PBS -t 1-100  # Parallel processing of 100 BAM files. Thread stored in variable $PBS_ARRAYID

module load MACS/2.1.0.20151222
module load bedtools/2.17.0
# Also assumes Homer is installed in home directory

indir=NULL  # Directory where BAM files are stored

outbeds=NULL # Path to store BED files converted from BAM files
outpeaks=NULL # Path to store MACS2 peak files
outtags=NULL # Path to store Homer tags
outbgs=NULL # Path to store bedGraph files

cd ${indir}
# Locating and preparing BAM file
bamname=$(find ${indir} -name "*sorted.bam" -maxdepth 1 | head -$PBS_ARRAYID | tail -1);
bamToBed -i ${bamname} > ${outbeds}/$(basename ${bamname} .bam).bed;

# Peak calling
export thisn=$(basename ${bamname} .bam);
macs2 callpeak -t ${bamname} -f BAMPE -n $(basename ${bamname} .bam) -g 'mm' --outdir ${outpeaks} --nomodel;
perl -lane 'print $F[0],"\t",$F[1],"\t",$F[2],"\t$ENV{thisn}\t",$F[4],"\t",$F[5]' ${outpeaks}/${thisn}_peaks.narrowPeak  | bedtools sort > ${outpeaks}/${thisn}_peaks.bed;

# Building bedGraph for sample. This is used later to calculate accessibility statistics on an arbitrary window in the genome (i.e., beyond peaks)
samp=$(basename ${bamname} | sed 's/\_/	/g' | perl -lane 'print $F[0],"_",$F[1]');
sampname=$(basename ${bamname} .bam);
[[ ! -d "${outtags}/${samp}" ]] && mkdir "${outtags}/${samp}";
makeTagDirectory ${outtags}/${samp} ${bamname} -genome hg19 -single -fragLength 150;
makeUCSCfile ${outtags}/${samp} -fragLength 150 > ${outbgs}/${sampname}.bedGraph;
