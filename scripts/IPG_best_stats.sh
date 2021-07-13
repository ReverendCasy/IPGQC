#!/usr/bin/zsh

if [ ! -d IPG_best_stats ]; then
mkdir IPG_best_stats
fi

for prefix in $(ls $1 | cut -d'_' -f1 | sort | uniq); do
prefix=$(echo $prefix | sed 's/$/_/')
if [ $(ls $1 | v="$prefix" awk 'index($0, ENVIRON["v"])==1' | wc -l) -eq 1 ]; then
cp $1/$(ls $1 | v="$prefix" awk 'index($0, ENVIRON["v"])==1') IPG_best_stats
echo "$prefix: sole file"
else
#for file in ${$(ls $1 | grep "$prefix")[@]}; do
#if [ $(echo $file | grep -c 'GCF') -eq 1 ]; then
#cp $1/$file IPG_best_stats/
#fi
gcf=$(ls $1 | v="$prefix" awk 'index($0, ENVIRON["v"])==1'| grep 'GCF')
gca=$(ls $1 | v="$prefix" awk 'index($0, ENVIRON["v"])==1' | grep 'GCA')
echo "$gcf : $(wc -l $1/$gcf | cut -d' ' -f1), $gca : $(wc -l $1/$gca | cut -d' ' -f1)"
if [ $(wc -l $1/$gcf | cut -d' ' -f1) -eq 0 ] && [ $(wc -l $1/$gca | cut -d' ' -f1) -eq 0 ]; then
echo "Empty files!"
else
if [ $(wc -l $1/$gcf | cut -d' ' -f1) -gt $(wc -l $1/$gca | cut -d' ' -f1) ]; then
cp $1/$gcf IPG_best_stats
echo "$gcf wins"
else
cp $1/$gca IPG_best_stats
echo "$gca wins"
fi
#done
fi
fi
done
