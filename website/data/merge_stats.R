library(data.table)
library(dplyr)
library(readr)

d1 <- as.data.frame(fread("rnaseq-data-1.tsv.gz"))
d1_cols <- as.character(d1[2,])
d1_headers <- as.character(d1[1,])
for (i in which(as.character(d1[1,]) != "")) {
  if (i != 1) {
    ix <- seq(i, i + 4)
    d1_cols[ix] <- sprintf("%s_%s", d1_headers[i], d1_cols[ix])
  }
}
d1 <- d1[3:nrow(d1),]
colnames(d1) <- d1_cols

d2 <- as.data.frame(fread("rnaseq-data-2.tsv.gz"))
d2_cols <- as.character(d2[2,])
d2_headers <- as.character(d2[1,])
for (i in which(as.character(d2[1,]) != "")) {
  if (i != 1) {
    ix <- seq(i, i + 4)
    d2_cols[ix] <- sprintf("%s_%s", d2_headers[i], d2_cols[ix])
  }
}
d2 <- d2[3:nrow(d2),]
colnames(d2) <- d2_cols

d <- inner_join(
  d1,
  d2[,c(1,9:ncol(d2))],
  "ensembl_id"
)

write_tsv(d, "rnaseq-data.tsv.gz")
