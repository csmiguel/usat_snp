rm -rf mkdir data/intermediate/clumpak_comparison
mkdir data/intermediate/clumpak_comparison

for k in {2..8}
#loop 1: iteration across K
do
outdir=$(echo data/intermediate/clumpak_comparison/K$k)
mkdir $outdir
#list with runs for Kx within dataset
  find data/final/run_usat_pelo -name *K"$k"_*_f | while read str_run
  #loop 2: iteration across runs within K
  do
    #store repetition id in variable
  rep=$(echo $str_run | sed 's/.*\(rep\)/\1/')
    #create output name for str_run
  new=$(echo $outdir/usat_K"$k"_"$rep")
  cp $str_run $new
  done
  find data/final/run_dart_pelo -name *K"$k"_*_f | while read str_run
  #loop 2: iteration across runs within K
  do
    #store repetition id in variable
  rep=$(echo $str_run | sed 's/.*\(rep\)/\1/')
    #create output name for str_run
  new=$(echo $outdir/dart_K"$k"_"$rep")
  cp $str_run $new
  done
done
cd data/intermediate/ && zip -r clumpak_comparison.zip clumpak_comparison -x "*.DS_Store" && cd -
rm -r data/intermediate/clumpak_comparison
