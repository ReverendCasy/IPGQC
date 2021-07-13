#!/usr/bin/zsh

for species in $(ls $1 | cut -d'_' -f1,2 | sort | uniq); do
nspecies=$(echo $species | sed 's/_/ /g')
echo "$nspecies\t$(ls $1 | grep -c $species)" 
done
