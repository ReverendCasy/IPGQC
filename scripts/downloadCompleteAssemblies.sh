#!/bin/bash

if not [ ! -d "complete_assemblies" ]; then
mkdir complete_assemblies
fi

for id in $(cut -f2 Freq10CompleteGenome.tsv | grep -w -f - CompleteGenomeTable.tsv | cut -f1)
do
esearch -db assembly -query "$id" | \
 esummary | \
 xtract -pattern DocumentSummary -element FtpPath_GenBank | \
 grep 'GCA' | \
 awk -F'/' '{print $0"/"$NF"__genomic.fna.gz"}' | \
 xargs -n1 wget -P complete_assemblies
done
