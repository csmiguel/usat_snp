#!/bin/bash
#main pipeline in CLUMPAK
#/Applications/CLUMPAK/ should contain downloaded folders from http://clumpak.tau.ac.il/
#I have to run CLUMPAK from the app path instead of my project folder because CLUMPAK calls
#other scripts within the app folder from relative paths
h=$(pwd)
cd /Applications/CLUMPAK

perl CLUMPAK.pl --id 1 \
  --dir "$h"/data/intermediate/clumpak_comparison \
  --file "$h"/data/intermediate/clumpak_comparison.zip \
done
cd -
