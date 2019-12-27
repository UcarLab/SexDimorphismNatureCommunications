#!/bin/bash

# Submit from queue
# This reads bedGraph files and selects coordinates based on a query consisting of peaks
# Peak files are in BED format, consisting of a chr, start, end, and label tab-delimited columns
# bedGraph files produced by Homer annotatePeaks.pl
# Output file contains 8 columns: first four are coordinates and standardized accessibility value for intervals overlapping query peaks plus flanks, last for are query peak coordinates plus labels
# Output file is named after the original bedGraph, adding a _query suffix, and a .txt extension. A queried bedGraph can be reconstituted using the first 4 columns of this file
# Syntax: bedGraph_query.sh <input directory> <output directory> <query file> <minimum flank size in bp around peak ends> <sequence order of bedGraph file as returned by $(find -name *bedGraph) run in indir>

if [[ "$#" -lt 4 ]]
then
    echo "$(basename $0) [inDir] [outDir] [queryFile] [flankSize] [FID]"  1>&2
    echo "   [inDir]: path to directory that contains bedGraph files. bedGraph file names assumed to end in suffix bedGraph.gz" 1>&2
    echo "   [outDir]: path to directory where query results will be placed. A new directory is created with a name based on the query file name." 1>&2
    echo "   [queryFile]: path to query file. Assumed to have extension .qry" 1>&2
    echo "   [FID]: Number index for bedGraph file to be processed in a single run" 1>&2
    exit 1
fi

indir=$(echo $1 | sed 's:/$::g')
outdir=$(echo $2 | sed 's:/$::g')
qry=$(echo $3 | sed 's:/$::g')
FID=$(echo $4)

qrydir=$(basename $qry .qry)
[[ ! -d "$qrydir" ]] && mkdir "$qrydir"

cd $outdir
oneqry=$(basename $qry .qry)$FID.qry
onebg=$(find $indir -maxdepth 1 -name "*bedGraph.gz" | head -$FID | tail -1)
tmpbg=$(basename $onebg .gz).tmp
outbg=$(basename $onebg .gz)

zcat $onebg > ${qrydir}/${tmpbg}

perl -ne 'print if $.>1' ${qrydir}/${tmpbg} > ${qrydir}/${outbg}
rm ${qrydir}/${tmpbg}

# Intersect with working bedGraph
outqry=$(basename $outbg | sed s/\.bedGraph/_query\.txt/g)
bedtools intersect -wo -a ${qrydir}/${outbg} -b $oneqry | cut -f1-9 > ${qrydir}/${outqry}