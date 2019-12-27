#!/bin/bash

# Adds Homer annotations to consensus peaks extracted form raw counts matrix

cd /path/to/data/narrowPeaks
data=narrow_raw_whitelisted_filtered.txt

bed=$(basename $data .txt).bed
perl -ne 'print if $.>1' $data | cut -f1-3 > $bed
fullanno=$(basename $data .txt)_annofull.txt
anno=$(basename $data .txt)_annotated
head -1 $data | cut -f1-3 > h.txt

annotatePeaks.pl $bed hg19 > $fullanno
cut -f2-4,8,10,16,19 $fullanno | sed 's/ //g' | sed s/\(.*\)//g | sed 's/\.[0-9]	/	/g' | sed 's/\...	/	/g' > ${anno}.tmp
head -1 ${anno}.tmp | sed 's/^Chr	Start	End/chr	start	end/g' > h.txt
perl -ne 'print if $.>1' ${anno}.tmp | sortBed | perl -lane 'print $F[0],"\t",$F[1]-1,"\t",$F[2],"\t",$F[3],"\t",$F[4],"\t",$F[5],"\t",$F[6]' | cat h.txt - > ${anno}.txt
rm *tmp
rm h.txt
