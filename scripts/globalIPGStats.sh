#!/usr/bin/zsh

for file in $(ls $1); do
species=$(echo $file | cut -d'_' -f1,2)
assembly=$(echo $file | cut -d'_' -f3)
source=$(echo $file | cut -d'_' -f4 | cut -d'.' -f1)
echo "$species\t$assembly\t$source\t$(wc -l $1/$file | cut -d' ' -f1)"
done
