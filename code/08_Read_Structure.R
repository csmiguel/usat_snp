path <- "/Users/miguelcamacho/Dropbox/Research/CiBIO/NGC/dartseq_genotyping/structure/dartseq/run1"

read str stlog  cat K6_rep4.stlog | grep [00]: > lk.txt
h <- read.table("data/intermediate/lk.txt", skip =1)
plot(h$X0.068, type = "l")
