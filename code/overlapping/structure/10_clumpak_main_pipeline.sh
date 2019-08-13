#!/bin/bash
#main pipeline in CLUMPAK
#/Applications/CLUMPAK/ should contain downloaded folders from http://clumpak.tau.ac.il/
#I have to run CLUMPAK from the app path instead of my project folder because CLUMPAK calls
#other scripts within the app folder from relative paths
h=$(pwd)
cd /Applications/CLUMPAK

find "$h"/data/intermediate/shared_clumpak -name *.zip | while read str_zipped
do
dataset=$(echo $str_zipped | sed 's/^.*\///g' | sed 's/.zip//')
pop_file=$(echo "$dataset"_populations_file)
perl CLUMPAK.pl --id 1 \
  --dir "$h"/data/intermediate/shared_clumpak/$dataset \
  --file $str_zipped \
  --indtopop "$h"/data/intermediate/shared_clumpak/$pop_file
done
cd -

#for some reason in the distribution I have for CLUMPAK the populations_file
#has to be passed to --indtopop and not to --podtopop (as the manual says).
