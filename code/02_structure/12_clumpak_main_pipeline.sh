#!/bin/bash
#main pipeline in CLUMPAK
#/Applications/CLUMPAK/ should contain downloaded folders from http://clumpak.tau.ac.il/
#I have to run CLUMPAK from the app path instead of my project folder because CLUMPAK calls
#other scripts within the app folder from relative paths
h=$(pwd)
cd /Applications/CLUMPAK

find "$h"/data/intermediate/clumpak -name *.zip | while read str_zipped
do
dataset=$(echo $str_zipped | sed 's/^.*\///g' | sed 's/.zip//')
pop_file=$(echo "$dataset"_populations_file)
perl CLUMPAK.pl --id 1 \
  --dir "$h"/data/intermediate/clumpak/$dataset \
  --file $str_zipped \
  --podtopop "$h"/data/intermediate/clumpak/$pop_file
done
cd -
