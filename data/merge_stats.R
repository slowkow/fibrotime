library(data.table)
library(dplyr)
library(readr)

d1 <- as.data.frame(fread("rnaseq-data-1.tsv.gz"))
d2 <- as.data.frame(fread("rnaseq-data-2.tsv.gz"))

d <- inner_join(
  d1,
  d2[,c(1,9:ncol(d2))],
  "ensembl_id"
)

write_tsv(d, "rnaseq-data.tsv.gz")
