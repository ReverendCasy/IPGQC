#!/bin/bash

if [ ! -d "IPG_complete_stats" ]; then
mkdir IPG_complete_stats
fi

for id in $(cut -f2 Freq10CompleteGenome.tsv | grep -w -f - CompleteGenomeTable.tsv | cut -f1)
do
echo $id
accs=( $(esearch -db assembly -query $id | \
 esummary | \
 xtract -pattern DocumentSummary -element Id Genbank RefSeq | \
 grep $id | \
 cut -f2,3) )

for acc in ${accs[@]}
do
echo $acc
prefix=$(echo $acc | cut -d'_' -f1)
esearch -db ipg -query $acc | \
 esummary | \
 xtract -pattern DocumentSummary -element Id > IPG_complete_stats/$id\_$prefix.txt
done

done

