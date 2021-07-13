#!/bin/bash
esearch -db assembly -query 'txid2[Organism:exp] AND \
(latest[filter] AND complete genome[filter] AND all[filter] NOT anomalous[filter])' | \
esummary -db assembly | xtract -pattern DocumentSummary -element Id SpeciesName > \
CompleteGenomeTable.tsv
