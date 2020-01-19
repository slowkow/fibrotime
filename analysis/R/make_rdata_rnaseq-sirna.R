#!/usr/bin/env Rscript
# make_rdata_rnaseq-sirna.R
# 2017-07-11

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  data.table,
  dplyr,
  pbapply,
  readr,
  readxl,
  stringr
)

# STAR RSEM created by Broad Technology Labs
# ----------------------------------------------------------------------------

rsem_file <- "data/rnaseq-sirna_rsem_tpm.rds"

if (!file.exists(rsem_file)) {
  filenames <- Sys.glob(
    "/data/srlab/slowikow/fibroblast_timecourse/data-raw/rnaseq-sirna/SN*/btl/projects/SSF/SmartSeq2/SSF*/SSF*_SmartSeqAnalysisV1.1_1/*/RSEM/*.genes.results"
  )
  dat <- rbindlist(pblapply(filenames, function(filename) {
    x        <- fread(filename)
    x$sample <- strsplit(basename(filename), "\\.")[[1]][1]
    return(x)
  }))
  tpm           <- as.data.frame(dcast(dat, gene_id ~ sample, value.var = "TPM"))
  rownames(tpm) <- tpm$gene_id
  tpm$gene_id   <- NULL
  tpm           <- as.matrix(tpm)
  # Drop version numbers from gene ids.
  rownames(tpm) <- sapply(strsplit(rownames(tpm), "\\."), "[", 1)
  saveRDS(tpm, rsem_file)
}

# Metadata
# ----------------------------------------------------------------------------

# 000008520469 is Plate 1
meta1 <- read_excel(
  path = "data-raw/rnaseq-sirna/2017-05-08_sample-map.xlsx",
  sheet = "Sample List-1"
)
meta1$Position <- sprintf("%s_%s%s",
  "000008520469",
  substr(meta1$Position, 1, 1),
  sprintf("%02i", as.integer(substr(meta1$Position, 2, 3)))
)
# 000008519369 is Plate 2
meta2 <- read_excel(
  path = "data-raw/rnaseq-sirna/2017-05-08_sample-map.xlsx",
  sheet = "Sample List-2"
)
meta2$Position <- sprintf("%s_%s%s",
  "000008519369",
  substr(meta2$Position, 1, 1),
  sprintf("%02i", as.integer(substr(meta2$Position, 2, 3)))
)
meta <- rbind(meta1, meta2)
rm(meta1, meta2)

colnames(meta) <- c(
  "sample", "info", "buffer", "volume", "concentration",
  "species", "celltype"
)
meta <- cbind(meta, str_split_fixed(meta$info, ", ", 4))
colnames(meta) <- c(
  "sample", "info", "buffer", "volume", "concentration",
  "species", "celltype", "donor", "sirna", "time", "stimulation"
)
levels(meta$stimulation) <- c("None", "TNF", "TNF_IL17")
meta$time          <- as.integer(str_replace(meta$time, "HR", ""))
meta$timefactor    <- factor(meta$time)
meta$buffer        <- NULL
meta$volume        <- NULL
meta$concentration <- NULL
meta$info          <- NULL
meta$species       <- NULL
meta$celltype      <- NULL
# Mask donor names
mask <- c(
  "RA1357"  = "RA5",
  "RA16425" = "RA6",
  "RA355"   = "RA7",
  "RA3916"  = "RA8"
)
meta$donor <- mask[meta$donor]
saveRDS(meta, "data/rnaseq-sirna-meta.rds")

# Kallisto Ensembl 89
# ----------------------------------------------------------------------------

abundance <- fread("gunzip -c data-raw/rnaseq-sirna/kallisto_ensembl89/abundance.tsv.gz")

kallisto_mat <- function(abundance, value.var = 'est_counts') {
  retval <- dcast(abundance, target_id ~ sample, value.var = value.var, fill = 0)
  retval <- as.data.frame(retval)
  rownames(retval) <- retval$target_id
  retval$target_id <- NULL
  retval
}

countst <- kallisto_mat(abundance, 'est_counts')
tpmt    <- kallisto_mat(abundance, 'tpm')

d <- fread(
  "gunzip -c data-raw/ensembl89/Homo_sapiens.GRCh38.cdna.all.filtered.tsv.gz"
)
d <- unique(d[,c("ensembl_gene_id", "ensembl_transcript_id")])
transcript_to_gene <- split(d$ensembl_gene_id, d$ensembl_transcript_id)
transcript_to_gene <- unlist(transcript_to_gene)

x <- as.data.table(data.frame(gene = transcript_to_gene[rownames(countst)], countst))
x <- x[, lapply(.SD, sum), by = gene, .SDcols = colnames(x)[2:ncol(x)]]
counts <- as.data.frame(x[,2:ncol(x)])
rownames(counts) <- x$gene
# Drop the version number.
rownames(counts) <- str_split_fixed(rownames(counts), "\\.", 2)[,1]
rownames(countst) <- str_split_fixed(rownames(countst), "\\.", 2)[,1]
colnames(counts) <- colnames(countst)

x <- as.data.table(data.frame(gene = transcript_to_gene[rownames(tpmt)], tpmt))
x <- x[, lapply(.SD, sum), by = gene, .SDcols = colnames(x)[2:ncol(x)]]
tpm <- as.data.frame(x[,2:ncol(x)])
rownames(tpm) <- x$gene
# Drop the version number.
rownames(tpm) <- str_split_fixed(rownames(tpm), "\\.", 2)[,1]
rownames(tpmt) <- str_split_fixed(rownames(tpmt), "\\.", 2)[,1]
colnames(tpm) <- colnames(tpmt)

# #' @param dat A numeric matrix or data.frame.
# #' @param xs A vector of groups (e.g. gene names).
# #' @return A data.table with the aggregated mean for each group.
# #' @seealso stats::aggregate
# sum_by <- function(dat, xs) {
#   # Convert to data.table.
#   dat <- data.table(dat)
#   # Append the vector of group names as an extra column.
#   dat$agg_var <- xs
#   # Melt the data.table so all values are in one column called "value".
#   dat <- melt(dat, id.vars = "agg_var")
#   # Cast the data.table back into the original shape.
#   dat <- dcast.data.table(
#     dat, agg_var ~ variable, value.var = "value",
#     fun.aggregate = sum, na.rm = TRUE
#   )
#   rs <- dat$agg_var
#   dat <- as.data.frame(dat)
#   rownames(dat) <- rs
#   dat$agg_var <- NULL
#   return(dat)
# }
# 
# gene_ids <- unlist(transcript_to_gene[rownames(tpmt)])
# tpm      <- sum_by(tpmt, gene_ids)
# tpmt     <- as.data.frame(tpmt)
# 
# # Drop the version numbers from the Ensembl identifiers.
# rownames(tpmt) <- str_split_fixed(rownames(tpmt), "\\.", 2)[,1]
# rownames(tpm)  <- str_split_fixed(rownames(tpm), "\\.", 2)[,1]

stopifnot(all(sort(meta$sample) == sort(colnames(counts))))
stopifnot(all(sort(meta$sample) == sort(colnames(countst))))
stopifnot(all(sort(meta$sample) == sort(colnames(tpm))))
stopifnot(all(sort(meta$sample) == sort(colnames(tpmt))))

counts  <- counts[,meta$sample]
countst <- countst[,meta$sample]
tpm     <- tpm[,meta$sample]
tpmt    <- tpmt[,meta$sample]

stopifnot(all.equal(colnames(tpm), meta$sample))
stopifnot(all.equal(colnames(tpmt), meta$sample))

myprint <- function(...) {
  print(substitute(...))
  print(...)
  cat("\n")
}

myprint(counts[1:5,1:5])
myprint(countst[1:5,1:5])

myprint(tpm[1:5,1:5])
myprint(tpmt[1:5,1:5])

myprint(meta[1:5,])

# Save the data
# ----------------------------------------------------------------------------
save(list = c("counts", "countst", "tpm", "tpmt", "meta"), file = "data/rnaseq-sirna_kallisto_ensembl89.rda")

# Write supplementary files for NCBI GEO
# ----------------------------------------------------------------------------

write_matrix <- function(x, path) {
  x <- data.frame(
    ID_REF = rownames(x),
    signif(x),
    stringsAsFactors = FALSE
  )
  rownames(x) <- NULL
  readr::write_tsv(x, path)
}

out_dir <- "share/NCBI-GEO/rnaseq-data-2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_matrix(counts, sprintf("%s/rnaseq_gene_counts.tsv.gz", out_dir))
write_matrix(countst, sprintf("%s/rnaseq_transcript_counts.tsv.gz", out_dir))

write_matrix(tpm, sprintf("%s/rnaseq_gene_tpm.tsv.gz", out_dir))
write_matrix(tpmt, sprintf("%s/rnaseq_transcript_tpm.tsv.gz", out_dir))

meta_cols <- c(
  "sample",
  "donor",
  "stimulation",
  "sirna",
  "time"
)
x <- meta[,meta_cols]
readr::write_tsv(x, sprintf("%s/rnaseq_metadata.tsv.gz", out_dir))

