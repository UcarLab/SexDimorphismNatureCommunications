#!/bin/bash

#PBS -w /home/emarquez
#PBS -l nodes=1:ppn=1,mem=8gb,walltime=6:00:00 
#PBS -m a
#PBS -t 1-120
#PBS -M Eladio.Marquez@jax.org 
#PBS -e /home/emarquez/scripts/bgqueryerr
#PBS -o /home/emarquez/scripts/bgqueryout

# calls script to extract signal data from bedGraph files for a set of query peaks

module load bedtools/2.17.0

flank=250000 # flanking distance in bp to extract bedGraph signal relative to each query peak

# one BED-formatted file per query. Query files are named query<n>.qry, with n being any number (1-8 in the loop below)
# qry files have 5 columns: chr|start|end|target gene|relevant cell type for peak (relevant perhaps because it is cell-specific or because its connection to the target gene is specific to this population)

for qno in `seq 1 8`; do 

	bgdir=/projects/emarquez/bedGraph/PBMC  # repository of all PBMC bedGraph files, precalculated
	qrydir=/fastscratch/emarquez/bedGraph/queries  # output/work directory
	oqry=${qrydir}/query${qno}.qry
	qry=${qrydir}/$(basename $oqry .qry)${PBS_ARRAYID}.qry

	cd $qrydir
	# Query chromosome sizes for hg19
	[[ ! -f "hg19.genome" ]] && mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo" | sed s/\|//g | perl -ne 'print if $.>1' > hg19.genome

	# First extend query peaks by flank and then run query using this modified query file

	bedtools slop -i $oqry -g hg19.genome -b $flank > $qry  # used to add flanks prior to query bedGraphs

	bash ~/scripts/bedGraph_query.sh $bgdir $qrydir $oqry $PBS_ARRAYID

# 	rm $qry
done