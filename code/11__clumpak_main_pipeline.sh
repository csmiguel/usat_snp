#main pipeline
#/Applications/CLUMPAK/ should contain downloaded folders from http://clumpak.tau.ac.il/
find data/intermediate/clumpak -name *.zip | while read str_zipped
do
outdir=(echo $str_zipped | sed 's/^.*\///g')
perl /Applications/CLUMPAK/CLUMPAK.pl --id 1 \
--dir data/intermediate/clumpak/"$outdir" \
--file str_zipped
done
