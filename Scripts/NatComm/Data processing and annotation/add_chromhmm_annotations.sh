#!/bin/bash

# Creates annotation files based on chromHMM states from various sources for consensus peaks extracted from raw counts matrix
# ChromHMM annotations extracted from segments BED files obtained from the sources described in the manuscript.
# Output consists of an annotation file for each profiled cell type, scored on a common set of consensus peaks. Multiple annotations per
# peak are possible and resolved in the R-based pipeline (see ATACseq_data_load.R script)

cd /path/to/data/narrowPeaks/chromHMM
data=/path/to/data/narrowPeaks/narrow_raw_whitelisted_filtered.txt

# Roadmap Epigenome
# 18-state model
for hmmanno_segments in /path/to/Datasets/Roadmap/ChromHMM_18/E029_18_core_K27ac_segments.bed \
		/path/to/Datasets/Roadmap/ChromHMM_18/E032_18_core_K27ac_segments.bed \
		/path/to/Datasets/Roadmap/ChromHMM_18/E037_18_core_K27ac_segments.bed \
		/path/to/Datasets/Roadmap/ChromHMM_18/E038_18_core_K27ac_segments.bed \
		/path/to/Datasets/Roadmap/ChromHMM_18/E041_18_core_K27ac_segments.bed \
		/path/to/Datasets/Roadmap/ChromHMM_18/E042_18_core_K27ac_segments.bed \
		/path/to/Datasets/Roadmap/ChromHMM_18/E044_18_core_K27ac_segments.bed \
		/path/to/Datasets/Roadmap/ChromHMM_18/E046_18_core_K27ac_segments.bed \
		/path/to/Datasets/Roadmap/ChromHMM_18/E047_18_core_K27ac_segments.bed \
		/path/to/Datasets/Roadmap/ChromHMM_18/E048_18_core_K27ac_segments.bed \
		/path/to/Datasets/Roadmap/ChromHMM_18/E050_18_core_K27ac_segments.bed \
		/path/to/Datasets/Roadmap/ChromHMM_18/E062_18_core_K27ac_segments.bed; do
	cs=$(basename $hmmanno_segments | sed s/\_.*//g)
	bed=$(basename $data .txt)_${cs}.bed
	perl -ne 'print if $.>1' $data | cut -f1-3 > $bed
	fullanno=$(basename $data .txt)_${cs}_annofull.txt
	anno=$(basename $data .txt)_${cs}_annotated
	head -1 $data | cut -f1-3 | sed 's/$/	chromHMMstate/g' > h.txt
	bedtools intersect -wo -a $bed -b $hmmanno_segments | cut -f1-3,7 | sortBed | uniq > $(basename $data .txt)_${cs}_chromHMM.tmp
	sed s/E18/Quies/g $(basename $data .txt)_${cs}_chromHMM.tmp \
		| sed s/E17/ReprPCWk/g | sed s/E16/ReprPC/g | sed s/E15/EnhBiv/g | sed s/E14/TssBiv/g | sed s/E13/Het/g | sed s/E12/ZNF_Rpts/g \
		| sed s/E11/EnhWk/g | sed s/E10/EnhA2/g | sed s/E9/EnhA1/g | sed s/E8/EnhG2/g | sed s/E7/EnhG1/g | sed s/E6/TxWk/g | sed s/E5/Tx/g \
		| sed s/E4/TssFlnkD/g | sed s/E3/TssFlnkU/g | sed s/E2/TssFlnk/g | sed s/E1/TssA/g | cat h.txt - \
		> $(basename $data .txt)_${cs}_chromHMM.txt
	rm *tmp
	rm h.txt
	rm $bed
done
# 15-state model
for hmmanno_segments in /path/to/Datasets/Roadmap/ChromHMM_15/E030_15_coreMarks_segments.bed \
		/path/to/Datasets/Roadmap/ChromHMM_15/E035_15_coreMarks_segments.bed; do
	cs=$(basename $hmmanno_segments | sed s/\_.*//g)
	bed=$(basename $data .txt)_${cs}.bed
	perl -ne 'print if $.>1' $data | cut -f1-3 > $bed
	fullanno=$(basename $data .txt)_${cs}_annofull.txt
	anno=$(basename $data .txt)_${cs}_annotated
	head -1 $data | cut -f1-3 | sed 's/$/	chromHMMstate/g' > h.txt
	bedtools intersect -wo -a $bed -b $hmmanno_segments | cut -f1-3,7 | sortBed | uniq > $(basename $data .txt)_${cs}_chromHMM.tmp
	sed s/E15/Quies/g $(basename $data .txt)_${cs}_chromHMM.tmp \
		| sed s/E14/ReprPCWk/g | sed s/E13/ReprPC/g | sed s/E12/EnhBiv/g | sed s/E11/BivTssFlnk/g | sed s/E10/TssBiv/g | sed s/E9/Het/g \
		| sed s/E8/ZNF_Rpts/g | sed s/E7/Enh/g | sed s/E6/EnhG/g | sed s/E5/TxWk/g | sed s/E4/Tx/g | sed s/E3/TxFlnk/g | sed s/E2/TssAFlnk/g \
		| sed s/E1/TssA/g | cat h.txt - \
		> $(basename $data .txt)_${cs}_chromHMM.txt
	rm *tmp
	rm h.txt
	rm $bed
done

# Blueprint epigenome
for hmmanno_segments in /path/to/Datasets/Blueprint/ChromHMM_12_hg19/S002R5H1_12_Blueprint_release_082014_segments.bed \
		/path/to/Datasets/Blueprint/ChromHMM_12_hg19/S00622H1_12_Blueprint_release_082014_segments.bed \
		/path/to/Datasets/Blueprint/ChromHMM_12_hg19/S001S7H2_12_Blueprint_release_082014_segments.bed \
		/path/to/Datasets/Blueprint/ChromHMM_12_hg19/S004BTH2_12_Blueprint_release_082014_segments.bed; do
	cs=$(basename $hmmanno_segments | sed s/\_.*//g)
	bed=$(basename $data .txt)_${cs}.bed
	perl -ne 'print if $.>1' $data | cut -f1-3 > $bed
	fullanno=$(basename $data .txt)_${cs}_annofull.txt
	anno=$(basename $data .txt)_${cs}_annotated
	head -1 $data | cut -f1-3 | sed 's/$/	chromHMMstate/g' > h.txt
	bedtools intersect -wo -a $bed -b $hmmanno_segments | cut -f1-3,7 | sortBed | uniq > $(basename $data .txt)_${cs}_chromHMM.tmp
	sed s/E12/TssA2/g $(basename $data .txt)_${cs}_chromHMM.tmp \
		| sed s/E11/TssA1/g | sed s/E10/TssADistal/g | sed s/E9/EnhA1/g | sed s/E8/EnhWk/g | sed s/E7/EnhG1/g | sed s/E6/TxWk/g \
		| sed s/E5/Tx/g | sed s/E4/Het/g | sed s/E3/Quies/g | sed s/E2/ReprPCWk/g | sed s/E1/ReprPC/g | cat h.txt - \
		> $(basename $data .txt)_${cs}_chromHMM.txt
	rm *tmp
	rm h.txt
	rm $bed
done
# CLL reference epigenome
for hmmanno_segments in /path/to/Datasets/CLL_ref_epigenome/ChromHMM_12_hg19/csMBC1_12_segments.bed \
		/path/to/Datasets/CLL_ref_epigenome/ChromHMM_12_hg19/ncsMBC_12_segments.bed \
		/path/to/Datasets/CLL_ref_epigenome/ChromHMM_12_hg19/GCBC3_12_segments.bed \
		/path/to/Datasets/CLL_ref_epigenome/ChromHMM_12_hg19/NBCB1_12_segments.bed \
		/path/to/Datasets/CLL_ref_epigenome/ChromHMM_12_hg19/PCT2_12_segments.bed; do
	cs=$(basename $hmmanno_segments | sed s/\_.*//g)
	bed=$(basename $data .txt)_${cs}.bed
	perl -ne 'print if $.>1' $data | cut -f1-3 > $bed
	fullanno=$(basename $data .txt)_${cs}_annofull.txt
	anno=$(basename $data .txt)_${cs}_annotated
	head -1 $data | cut -f1-3 | sed 's/$/	chromHMMstate/g' > h.txt
	bedtools intersect -wo -a $bed -b $hmmanno_segments | cut -f1-3,7 | sortBed | uniq > $(basename $data .txt)_${cs}_chromHMM.tmp
	sed s/E12/Repr2/g $(basename $data .txt)_${cs}_chromHMM.tmp \
		| sed s/E11/Quies/g | sed s/E10/Repr1/g | sed s/E9/TxElong/g | sed s/E8/TxWk/g | sed s/E7/Tx/g | sed s/E6/EnhA2/g \
		| sed s/E5/EnhWk/g | sed s/E4/TssBiv/g | sed s/E3/TssWk/g | sed s/E2/EnhA1/g | sed s/E1/TssA/g | cat h.txt - \
		> $(basename $data .txt)_${cs}_chromHMM.txt
	rm *tmp
	rm h.txt
	rm $bed
done
