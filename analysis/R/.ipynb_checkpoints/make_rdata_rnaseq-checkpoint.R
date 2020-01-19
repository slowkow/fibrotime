#!/usr/bin/env Rscript
# make_rdata_rnaseq.R
# 2017-07-11

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    data.table,
    dplyr,
    readr,
    readxl,
    reshape2,
    stringr
)

abundance <- fread("gunzip -c data-raw/rnaseq/kallisto_ensembl89/abundance.tsv.gz")
head(abundance)
##                              sample          target_id length eff_length
## 1: Pl1A0113_10-B-2_OA452TNF1__IL171  ENST00000358097.8   1655   1473.000
## 2: Pl1A0113_10-B-2_OA452TNF1__IL171  ENST00000610593.4   1363   1181.000
## 3: Pl1A0113_10-B-2_OA452TNF1__IL171  ENST00000534447.5   5330   5148.000
## 4: Pl1A0113_10-B-2_OA452TNF1__IL171 ENST00000299598.11   5818   5636.000
## 5: Pl1A0113_10-B-2_OA452TNF1__IL171  ENST00000621118.4   1112    930.003
## 6: Pl1A0113_10-B-2_OA452TNF1__IL171  ENST00000457194.6   1984   1802.000
##     est_counts         tpm
## 1: 2.56476e+00 1.04088e+00
## 2: 1.98146e+00 1.00297e+00
## 3: 1.15941e+01 1.34634e+00
## 4: 1.01725e+00 1.07898e-01
## 5: 2.80320e-05 1.80188e-05
## 6: 4.76865e+00 1.58196e+00

kallisto_mat <- function(abundance, value.var = 'est_counts') {
  retval <- dcast(abundance, target_id ~ sample, value.var = value.var, fill = 0)
  retval <- as.data.frame(retval)
  rownames(retval) <- retval$target_id
  retval$target_id <- NULL
  retval
}

countst <- kallisto_mat(abundance, 'est_counts')
tpmt    <- kallisto_mat(abundance, 'tpm')

# Extract metadata from the sample names.
m2             <- data.frame(str_split_fixed(colnames(tpmt), "_", 5), stringsAsFactors = FALSE)
m2$X4          <- NULL
m2             <- m2[order(m2[,2]),]
colnames(m2)   <- c("X1", "ID", "X3", "X5")
m2$Cell_Line   <- substr(m2[,3], 1, 5)
m2$Stimulation <- paste(substr(m2[,3], 6, 10), m2$X5)

# Ensure this is a unique identifier for each sample.
stopifnot(all(table(m2$ID) == 1))

colnames(tpmt) <- sapply(str_split(colnames(tpmt), "_"), "[", 2)
tpmt           <- tpmt[,order(colnames(tpmt))]

colnames(countst) <- sapply(str_split(colnames(countst), "_"), "[", 2)
countst           <- countst[,order(colnames(countst))]

stopifnot(all.equal(m2$ID, colnames(tpmt)))

# Read metadata from the plate design table.
m <- read_csv("data-raw/rnaseq/plate_design_final.csv")
# Delete duplicated column.
all(m[,3] == m[,16])
m[,3] <- NULL
colnames(m) <- c(
  "Plate", "Well", "ID", "Cell_Line", "Stimulation", "Time", 
  "Concentration", "RNA_Prep", "Disease", "Conc_High", "Unused_Well", 
  "BLANK", "Plate#", "Target_Well", "Sample", "Conc_ng_ul", "75ng_ul", 
  "Add_H2O_to_15ul"
)
# Delete duplicated column.
m[['Plate#']] <- NULL
# Order by ID.
m <- m[!is.na(m$ID),]
m <- m[order(m$ID),]
m$Row <- substr(m$Well, 1, 1)
m$Col <- substr(m$Well, 2, 10)

# Sanity checks.
stopifnot(all.equal(m2$Cell_Line, m$Cell_Line))
table(m2$Stimulation, m$Stimulation)
stopifnot(all.equal(m$ID, colnames(tpmt)))
stopifnot(all.equal(m$ID, colnames(countst)))

m$Sample <- sprintf("S%s", m$Sample)
# Set column names to a short identifier.
colnames(tpmt) <- m$Sample
colnames(countst) <- m$Sample

# Sort by the identifier.
m <- m[order(m$Sample), ]

tpmt <- tpmt[, order(colnames(tpmt))]
stopifnot(all.equal(m$Sample, colnames(tpmt)))

countst <- countst[, order(colnames(countst))]
stopifnot(all.equal(m$Sample, colnames(countst)))

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

x <- as.data.table(data.frame(gene = transcript_to_gene[rownames(tpmt)], tpmt))
x <- x[, lapply(.SD, sum), by = gene, .SDcols = colnames(x)[2:ncol(x)]]
tpm <- as.data.frame(x[,2:ncol(x)])
rownames(tpm) <- x$gene
# Drop the version number.
rownames(tpm) <- str_split_fixed(rownames(tpm), "\\.", 2)[,1]
rownames(tpmt) <- str_split_fixed(rownames(tpmt), "\\.", 2)[,1]

##' @param dat A numeric matrix or data.frame.
##' @param xs A vector of groups (e.g. gene names).
##' @return A data.table with the aggregated mean for each group.
##' @seealso stats::aggregate
#sum_by <- function(dat, xs) {
#  # Convert to data.table.
#  dat <- data.table(dat)
#  # Append the vector of group names as an extra column.
#  dat$agg_var <- xs
#  # Melt the data.table so all values are in one column called "value".
#  dat <- melt(dat, id.vars = "agg_var")
#  # Cast the data.table back into the original shape.
#  dat <- dcast.data.table(
#    dat, agg_var ~ variable, value.var = "value",
#    fun.aggregate = sum, na.rm = TRUE
#  )
#  rs <- dat$agg_var
#  dat <- as.data.frame(dat)
#  rownames(dat) <- rs
#  dat$agg_var <- NULL
#  return(dat)
#}
#
#gene_ids <- unlist(transcript_to_gene[rownames(tpmt)])
#tpm <- sum_by(tpmt, gene_ids)
#
#tpmt <- as.data.frame(tpmt)
#
## Drop the version numbers from the Ensembl identifiers.
#rownames(tpmt) <- str_split_fixed(rownames(tpmt), "\\.", 2)[,1]
#rownames(tpm)  <- str_split_fixed(rownames(tpm), "\\.", 2)[,1]

# Edit the metadata dataframe
# ----------------------------------------------------------------------------
m$BLANK <- NULL
m$TimeStim <- factor(paste(m$Time, m$Stimulation))

m$TNF <- as.integer(str_detect(m$Stimulation, "TNF"))
m$IL17 <- 0
m$IL17[str_detect(m$Stimulation, "IL17 \\(1\\)")] <- 1
m$IL17[str_detect(m$Stimulation, "IL17 \\(10\\)")] <- 10

doses <- c(
  "None"                = 0,
  "TNF (1)"             = 0,
  "TNF (1) + IL17 (1)"  = 1,
  "TNF (1) + IL17 (10)" = 10
)
m$Dose <- doses[m$Stimulation]
m$RA <- as.logical(m$Disease == "RA")

m$DoseFactor   <- factor(m$Dose)
m$DoseExposure <- m$Dose > 0
m$TimeFactor   <- factor(m$Time)

# Donor information
d <- readxl::read_excel("data-raw/rnaseq/Cell_Information.xlsx", skip = 1)
d <- d[,c(1,3,4,5,6)]
colnames(d) <- c("Cell_Line", "Joint", "Sex", "Age", "RF")
d$Disease <- "OA"
d$Disease[1:5] <- "RA"
d <- d[d$Cell_Line %in% c("552", "553", "554", "704", "502", "504", "452"),]
d$Cell_Line <- paste(d$Disease, d$Cell_Line, sep = "")

m <- cbind.data.frame(m, d[match(m$Cell_Line, d$Cell_Line),-1])
m$path <- NULL

# Mask cell line identifiers.
m$Cell_Line[m$Cell_Line == "RA552"] <- "RA1"
m$Cell_Line[m$Cell_Line == "RA553"] <- "RA2"
m$Cell_Line[m$Cell_Line == "RA554"] <- "RA3"
m$Cell_Line[m$Cell_Line == "RA704"] <- "RA4"

m$Cell_Line[m$Cell_Line == "OA502"] <- "OA1"
m$Cell_Line[m$Cell_Line == "OA504"] <- "OA2"
m$Cell_Line[m$Cell_Line == "OA452"] <- "OA3"

myprint <- function(...) {
  print(substitute(...))
  print(...)
  cat("\n")
}

myprint(counts[1:5,1:5])
myprint(countst[1:5,1:5])

myprint(tpm[1:5,1:5])
myprint(tpmt[1:5,1:5])

myprint(m[1:5,])

# Save the data
# ----------------------------------------------------------------------------
# save(list = c("tpm", "tpmt", "m"), file = "data/rnaseq_kallisto_ensembl89.rda")
save(list = c("counts", "countst", "tpm", "tpmt", "m"), file = "data/rnaseq_kallisto_ensembl89.rda")

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

out_dir <- "share/NCBI-GEO/rnaseq-data-1"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_matrix(counts, sprintf("%s/rnaseq_gene_counts.tsv.gz", out_dir))
write_matrix(countst, sprintf("%s/rnaseq_transcript_counts.tsv.gz", out_dir))

write_matrix(tpm, sprintf("%s/rnaseq_gene_tpm.tsv.gz", out_dir))
write_matrix(tpmt, sprintf("%s/rnaseq_transcript_tpm.tsv.gz", out_dir))

meta_cols <- c(
  "Sample",
  "Plate",
  "Well",
  "Cell_Line",
  "Stimulation",
  "Dose",
  "Time",
  "Disease",
  "Sex"
)
x <- m[,meta_cols]
readr::write_tsv(x, sprintf("%s/rnaseq_metadata.tsv.gz", out_dir))

