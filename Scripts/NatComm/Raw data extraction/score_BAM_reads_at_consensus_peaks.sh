#!/bin/bash

module load bedtools/2.17.0
module load python/2.7.3

bamdir=NULL # Directory containing BAM files
outdir=NULL # Output directory
consensus=NULL # Path to consensus file, in BED format
samplist=NULL # Path to text file with paths to sample peak files in BED format, one sample per line
cd ${outdir}

>bamdirs.tmp
for samp in $(cat $samplist); do
	dirname $(find ${bamdir} -name "*${samp}_*PBMC*bam") >> bamdirs.tmp
done
bamdirs=$(cat bamdirs.tmp | sort -u);
for d in $bamdirs; do
	thisn=$(basename $d);
	
	>bamnames0.tmp
	for samp in $(cat samplist); do
		find ${d} -name "*${samp}_*PBMC*bam" >> bamnames0.tmp
	done
	bamnames=$(cat bamnames0.tmp);

	>bamnames1.tmp;
	for f in ${bamnames};
	do
	  basename $f | sed 's/\_/	/g' | perl -lane 'print $F[0],"_",$F[1]' >> bamnames1.tmp;
	done;
	sed ':a;N;$!ba;s/\n/\t/g' bamnames1.tmp > bamnames2.tmp;
	
	sed ':a;N;$!ba;s/\n/\t/g' bamnames1.tmp | perl -ne 'print "chr\tstart\tend\tnhits\t",$_' > ${outdir}/peaks_rawcounts_${thisn}.txt;
	bedtools multicov -D -bams ${bamnames} -bed ${consensus} >> ${outdir}/peaks_rawcounts_${thisn}.txt;
done
rm "bamnames[0-2].tmp"