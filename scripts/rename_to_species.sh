#!/usr/bin/zsh

specs="${(@f)$(cut -f2 Freq10CompleteGenome.tsv)}"
#echo $specs

for line in "${(@f)$(cat CompleteGenomeTable.tsv)}"; do
acc=$(echo $line | cut -f1)
file=$(ls IPG_best_stats | v="$acc" awk 'index($0, ENVIRON["v"])==1')
#echo $file
#echo $acc
species=$(echo $line | cut -f2 | cut -d' ' -f1,2)
#echo $species
if [ ${specs[(i)$species]} -le ${#specs} ]; then
echo $acc
echo $species
newfile=$(echo $species | sed 's/\s/_/g')
newfile+=$(echo $file | sed 's/^/_/g')
echo $file
echo $newfile
mv IPG_best_stats/$file IPG_best_stats/$newfile
fi
done
