# Libraries -------------------------------------------------------------------

if (!require("pacman")) install.packages("pacman")

# Plotting
pacman::p_load(
    RColorBrewer,
    cowplot,
    ggplot2,
    ggrepel,
    gridExtra,
    pheatmap,
    scales,
    superheat,
    viridis,
    wiggleplotr,
    patchwork
)

# Reading and manipulating
pacman::p_load(
    GenomicRanges,
    GenomicFeatures,
    data.table,
    dplyr,
    magrittr,
    matrixStats,
    readr,
    readxl,
    reshape2,
    rlist,
    stringr
)

# Parallel computing
pacman::p_load(
    doMC,
    foreach,
    parallel,
    pbapply
)

# Modeling
pacman::p_load(
    NMF,
    limma,
    lme4,
    qvalue,
    mixtools
)

#

# Globals
# ----------------------------------------------------------------------------

load_message <- function(file, ...) {
  
  load(file, ...)
}

font_size <- 20

source("R/meta_colors.R")

sirnas <- c("CUX1", "ELF3", "LIFR", "STAT3", "STAT4")


# MSigDB
# ----------------------------------------------------------------------------

file_gene_sets <- "data/gene_sets.rda"

message(sprintf("load(%s)", file_gene_sets))
load(file_gene_sets)

# Dose and time RNA-seq data
# ----------------------------------------------------------------------------

file_rnaseq <- "data/rnaseq_kallisto_ensembl89.rda"
message(sprintf("load(%s)", file_rnaseq))
load(file_rnaseq)
log2tpm  <- log2(tpm + 1)
log2tpmt <- log2(tpmt + 1)
# Drop unexpressed genes.
ix_genes <- apply(tpm, 1, function(x) sum(x > 0) >= 3)
tpm      <- tpm[ix_genes,]
log2tpm  <- log2tpm[ix_genes,]
# Drop unexpressed transcripts.
ix_transcripts <- apply(tpmt, 1, function(x) sum(x > 0) >= 3)
tpmt           <- tpmt[ix_transcripts,]

# siRNA RNA-seq data
# ----------------------------------------------------------------------------

file_sirna <- "data/rnaseq-sirna_kallisto_ensembl89.rda"
message(sprintf("load(%s)", file_sirna))
sirna <- new.env()
load(file_sirna, envir = sirna)
sirna$meta$dosefactor <- factor(1 * (sirna$meta$stimulation == "TNF_IL17"))
levels(sirna$meta$stimulation) <- c("None", "TNF (1)", "TNF (1) + IL17 (1)")
sirna$log2tpm <- log2(sirna$tpm + 1)
stopifnot(all(sirna$meta$sample == colnames(sirna$log2tpm)))
# Drop unexpressed genes.
ix_genes      <- apply(sirna$log2tpm, 1, function(x) sum(x > 0) >= 3)
sirna$log2tpm <- sirna$log2tpm[ix_genes,]
sirna$tpm     <- sirna$tpm[ix_genes,]
# Drop unexpressed transcripts.
ix_transcripts <- apply(sirna$tpmt, 1, function(x) sum(x > 0) >= 3)
sirna$tpmt     <- sirna$tpmt[ix_transcripts,]
                        
#

# Gene symbols ---------------------------------------------------------------

file_gene_names <- "data/ensembl89/Homo_sapiens.GRCh38.89.gene_names.tsv.gz"
message("fread('gunzip -c %s')", file_gene_names)
gene_names <- fread(sprintf("gunzip -c %s", file_gene_names))
x                 <- unique(gene_names[,c("gene_id", "gene_name")])
ensembl_to_symbol <- unlist(split(x$gene_name, x$gene_id))
symbol_to_ensembl <- unlist(split(x$gene_id, x$gene_name))
ensembl_ids       <- union(rownames(tpm), rownames(sirna$tpm))
gene_symbols      <- ensembl_to_symbol[ensembl_ids]
