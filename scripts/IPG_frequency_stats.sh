#!/usr/bin/zsh

for line in $(cat $1); do
esearch -db ipg -query $line | \
 esummary | \
 xtract -pattern DocumentSummary \
        -element IPG Title AssemblyCount >> $2
done
