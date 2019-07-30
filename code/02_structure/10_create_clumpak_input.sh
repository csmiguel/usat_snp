#!/bin/bash
#put all your results in one folder ‘results_folder’
#Go to the command line terminal and type:
#1. prepare input file. It needs the below structure:
# Structure results zipped into a single file
# stresults.zip
# 	k1.zip
#		k1/
#			project_data_k1_run10_f
#			project_data_k1_run1_f
#			...
#	k2.zip
#		k2/
#			...
#1. 1. create folder STRUCTURE
rm -rf mkdir data/intermediate/clumpak
mkdir data/intermediate/clumpak

find data/final -name run_* | while read dataset
#loop 1: iteration across datasets
do
d=$(echo $dataset | sed 's/^.*\///g')
mkdir data/intermediate/clumpak/$d
  for k in {1..8}
  #loop 2: iteration across K
  do
  outdir=$(echo data/intermediate/clumpak/$d/K$k)
  mkdir $outdir
  #list with runs for Kx within dataset
    find $dataset -name *K"$k"_*_f | while read str_run
    #loop 3: iteration across runs within K
    do
      #store repetition id in variable
    rep=$(echo $str_run | sed 's/.*\(rep\)/\1/')
      #create output name for str_run
    new=$(echo $outdir/K"$k"_"$rep")
    cp $str_run $new
    done
  done
  cd data/intermediate/clumpak/ && zip -r "$d".zip "$d" -x "*.DS_Store" && cd -
  rm -r data/intermediate/clumpak/$d
done
