# Jupyter functions -----------------------------------------------------------

show_plot <- function(obj, name = NULL, width = NULL, height = NULL, units = "in", res = 120, save_pdf = TRUE) {
  out_dir <- "notebooks/figures"
  dir.create(out_dir, showWarnings = FALSE)

  if (is.null(width)) {
    width <- getOption("repr.plot.width", 3)
  }
  if (is.null(height)) {
    height <- getOption("repr.plot.height", 3)
  }
  if (is.null(name)) {
    name <- substr(tempfile("", tmpdir = ""), 2, 6)
  }

  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

  #filename <- sprintf("%s/%s_%s", out_dir, timestamp, name)
  filename <- sprintf("%s/%s", out_dir, name)
  file_png <- sprintf("%s.png", filename)
  file_pdf <- sprintf("%s.pdf", filename)

  message(sprintf("Writing %s", file_png))
  ggsave(
    filename = file_png,
    plot = obj, width = width, height = height, units = units, dpi = res
  )

  if (save_pdf) {
    message(sprintf("Writing %s", file_pdf))
    # pdf.options(encoding = 'ISOLatin2')
    ggsave(
      filename = file_pdf,
      plot = obj,
      device = cairo_pdf,
      width = width, height = height, units = units, dpi = res
    )
  }

  IRdisplay::display_png(file = file_png)
}

# Pure functions --------------------------------------------------------------

corner <- function(x) {
  nr <- min(5, nrow(x))
  nc <- min(5, ncol(x))
  x[1:nr, 1:nc]
}

# Read about Bessel's correction: https://en.wikipedia.org/wiki/Bessel%27s_correction
# The sd() function already uses (n-1) as the denominator.
# So, here we also use (n-1) as the denominator.
sem <- function(x) sd(x) / sqrt(length(x) - 1)

# This function can read bed narrowPeak files from ENCODE.
read_encode_bed <- function(file) {
  retval <- readr::read_tsv(
    file = file,
    col_names = FALSE,
    col_types = 'ciicicdddi'
  )
  colnames(retval) <- c(
    "chrom", "chromStart", "chromEnd",
    "name", "score", "strand",
    "signalValue", "pValue", "qValue",
    "peak"
  )
  GenomicRanges::GRanges(
    seqnames = retval$chrom,
    ranges = IRanges::IRanges(
      start = retval$chromStart,
      end = retval$chromEnd,
      names = retval$name
    ),
    strand = stringr::str_replace(retval$strand, "\\.", "*"),
    score = retval$score,
    signalValue = retval$signalValue,
    pValue = retval$pValue,
    qValue = retval$qValue
  )
}

sigfig <- function(vec, digits = 2){
  retval <- gsub("\\.$", "",
    formatC(
      x      = signif(x =vec, digits = digits),
      digits = digits,
      format = "fg",
      flag   = "#"
    )
  )
  ix_small <- which(vec < 1e-3)
  if (length(ix_small) != 0L) {
    retval[ix_small] <- gsub("\\.$", "",
      formatC(
        x      = signif(x = vec[ix_small], digits = digits),
        digits = digits,
        format = "g",
        flag   = "#"
      )
    )
  }
  # Convert 1e-06 to 1e-6
  retval <- gsub("e-0(\\d)$", "e-\\1", retval)
  retval
}

scale_rows <- function(x, ...) t(scale(t(x), ...))

sort_hclust <- function(...) { as.hclust(dendsort(as.dendrogram(...))) }

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

trim_image <- function(filename) {
  new_filename <- file.path(
    dirname(filename),
    sprintf("convert_trim_%s", basename(filename))
  )
  command <- sprintf(
    'convert "%s" -trim "%s" && mv -f "%s" "%s"',
    filename, new_filename,
    new_filename, filename
  )
  system(command)
}

optimize_png <- function(filename) {
  opt_filename <- sprintf(
    "%s-fs8.png", substr(filename, 1, nchar(filename) - 4)
  )
  command <- sprintf(
    'pngquant --ext -fs8.png -- "%s" && mv -f "%s" "%s"',
    filename, opt_filename, filename
  )
  system(command)
}

ggsave_optimize_png <- function(filename, overwrite = TRUE, ...) {
  if (!overwrite && file.exists(filename)) {
    message("Not overwriting ", filename)
  } else {
    message(filename)
  }
  ggplot2::ggsave(filename, ...)
  n <- nchar(filename)
  if (substr(filename, n - 3, n) == ".png") {
    optimize_png(filename)
  }
}

sort_hclust <- function(...) as.hclust(dendsort::dendsort(as.dendrogram(...)))

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

low_density <- function(x, y, n = 100, q = 0.05) {
  kde <- MASS::kde2d(x, y, n = n)
  kx <- cut(x, kde$x, labels = FALSE, include.lowest = TRUE)
  ky <- cut(y, kde$y, labels = FALSE, include.lowest = TRUE)
  kz <- sapply(seq_along(kx), function(i) kde$z[kx[i], ky[i]])
  kz < quantile(kz, q)
}

extreme_idx <- function(xs, low = 0.005, high = 0.995) {
  qs <- quantile(xs, probs = c(low, high))
  xs < qs[1] | xs > qs[2]
}

# Given an Ensembl ID like "ENSG00.1", return "ENSG00".
drop_version <- function(xs) stringr::str_split_fixed(xs, "\\.", 2)[,1]

do_fisher <- function(ids, universe, gene_sets) {
  retval <- as.data.frame(t(sapply(gene_sets, function(gene_set) {
    x <- sum(ids %in% gene_set)
    n <- length(ids)
    X <- sum(universe %in% gene_set)
    N <- length(universe)
    if (x > 0 && X > 0) {
      fish <- fisher.test(
        x = matrix(c(x, n - x, X - x, N - X - (n - x)), 2),
        alternative = "greater"
      )
      retvec <- c(
        "x" = x, "n" = n, "X" = X, "N" = N,
        "enrichment" = (x / n) / (X / N),
        "orlow" = fish$conf.int[1],
        "orhigh" = fish$conf.int[2],
        "or" = unname(fish$estimate),
        "pval" = fish$p.value
      )
    } else {
      retvec <- c(
        "x" = x, "n" = n, "X" = X, "N" = N,
        "enrichment" = 1,
        "orlow" = 1,
        "orhigh" = 1,
        "oddsratio" = 1,
        "pval" = 1
      )
    }
    retvec
  })))
  retval$Name <- names(gene_sets)
  retval <- retval[order(retval$pval),]
  retval$qval <- p.adjust(retval$pval, method = "fdr")
  retval
}

#

# Plotting functions ----------------------------------------------------------

pow_trans <- function(power = 2) {
  force(power)
  trans <- function(x) x ^ power
  inv <- function(x) x ^ (1 / power)
  scales::trans_new(paste0("power-", format(power)), trans, inv)
}
               
theme_clean <- function(base_size) {
  theme_classic(base_size = base_size) %+replace% theme(
    panel.border = element_rect(fill = NA, size = 0.2),
    axis.line.y = element_line(colour = "black", size = 0.2),
    axis.line.x = element_line(colour = "black", size = 0.2),
    axis.ticks = element_line(colour = "black", size = 0.2),
    strip.background = element_blank()
    #strip.background = element_rect(size = 0.25)
  )
}
  
# Overwrite default draw_colnames in the pheatmap package.©
# Thanks to Josh O'Brien at http://stackoverflow.com/questions/15505607
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

# https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r
get_breaks <- function(test, paletteLength) {
  myColor <- colorRampPalette(c("red", "white", "blue"))(paletteLength)
  # use floor and ceiling to deal with even/odd length palette lengths
  c(
    seq(min(test), 0, length.out = ceiling(paletteLength / 2) + 1), 
    seq(max(test) / paletteLength, max(test), length.out = floor(paletteLength / 2))
  )
}

plot_geneset_by_stimulation <- function(geneset, meta, log2tpm, ci = 0.5) {
  # Take the columns we need for this plot.
  dat <- meta[,c("Stimulation", "Time")]
  # Try to accept Ensembl IDs or HGNC symbols.
  if (!grepl("^ENSG", geneset[1])) {
    # Use the global variable with mappings from Ensembl to Symbol.
    geneset <- names(gene_symbols[which(gene_symbols %in% geneset)])
    stopifnot(length(geneset) > 1)
  }
  geneset <- geneset[geneset %in% rownames(log2tpm)]  
  dat <- cbind.data.frame(dat, t(log2tpm[geneset,]))
  # Add extra rows, only for the purpose of plotting the 0h time point.
  x <- dat[dat$Time == 0,]
  dat_zero <- rbind(x, x, x)
  stims <- unique(dat$Stimulation)
  stims <- stims[stims != "None"]
  dat_zero$Stimulation <- rep(stims, each = nrow(x))
  dat <- rbind(dat, dat_zero)
  dat <- subset(dat, Stimulation != "None")
  # Convert Time to a factor.
  dat$Time <- factor(dat$Time)
  dat_mean <- reshape2::melt(dat, id.vars = c("Stimulation", "Time")) %>%
    group_by(Stimulation, Time, variable) %>%
    summarise(
      mean = mean(value)
    ) %>% 
    group_by(Stimulation, Time) %>%
    summarise(
      median = median(mean),
      mean = mean(mean),
      low = mean(mean) - sem(mean),
      high = mean(mean) + sem(mean)
      # low = quantile(mean, probs = (1 - ci) / 2),
      # high = quantile(mean, probs = 1 - ((1 - ci) / 2))
    )
  ggplot(
    data = dat_mean,
    mapping = aes(
      x = Time,
      # y = median,
      y = mean,
      ymin = low,
      ymax = high,
      group = Stimulation,
      color = Stimulation
    )
  ) +
  geom_line(size = 1) +
  geom_ribbon(aes(fill = Stimulation), alpha = 0.15, color = NA) +
  scale_color_manual(values = brewer.pal(9, "YlOrRd")[c(4,7,9)]) +
  scale_fill_manual(values = brewer.pal(9, "YlOrRd")[c(4,7,9)]) +
  annotate(
    geom  = "text",
    label = sprintf("Summary of %s genes", scales::comma(length(geneset))),
    size = 5,
    x     = Inf,
    y     = -Inf,
    vjust = -1,
    hjust = 1.1
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  #theme_classic(base_size = font_size) +
  labs(
    x = "Hours",
    y = bquote("Log"[2]~"TPM")
  )
  # theme(
  #   panel.grid.major = element_line(color = 'grey60', size = 0.1),
  #   panel.grid.minor = element_blank()
  # )
}

plot_gene_sirna <- function(
  gene      = NULL,
  trend     = FALSE,
  save_plot = TRUE,
  show_plot = TRUE,
  facet_by  = "stimulation ~ sirna",
  color_by  = "donor",
  filename  = "figures/sirna/genes/%s.png"
) {
  # Use the global variable that contains metadata.
  dat <- sirna$meta
  # Try to accept Ensembl IDs or HGNC symbols.
  if (!grepl("^ENSG", gene)) {
    # Use the global variable with mappings from Ensembl to Symbol.
    gene <- names(which(gene_symbols == gene))
    stopifnot(length(gene) == 1)
  }
  dat[['Gene']] <- as.matrix(sirna$log2tpm)[gene,]
  # Add extra rows, only for the purpose of plotting the 0h time point.
  x <- dat[dat$time == 0,]
  dat_zero <- rbind(x, x)
  stims <- unique(dat$stimulation)
  stims <- stims[stims != "None"]
  dat_zero$stimulation <- rep(stims, each = nrow(x))
  dat <- rbind(dat, dat_zero)
  dat <- subset(dat, stimulation != "None")
  # Compute the averages.  
  dat_mean <- aggregate(Gene ~ time + stimulation, dat, mean)
  # Make a plot.
  p <- ggplot()
  if (trend == "lm") {
    p <- p + stat_smooth(
      data = dat,
      aes(x = time, y = Gene, group = stimulation),
      color = "grey50", method = "lm", size = 2
    )
  } else if (trend == "mean") {
    p <- p + geom_line(
      data = dat_mean,
      aes(x = time, y = Gene, group = stimulation),
      size = 1.5, color = "black"
    )
  }
  if (!is.null(facet_by)) {
    p <- p +
      facet_grid(as.formula(facet_by))
  }
  p <- p + geom_line(
      data = dat,
      aes_string(
        x = "time", y = "Gene",
        color = color_by, group = "donor"
      ),
      size = 1
    ) +
    scale_color_manual(values = brewer.pal(8, "Dark2")) +
    scale_x_continuous(breaks = sort(unique(dat$time))) +
    theme_bw(base_size = font_size) +
    labs(
      title = gene_symbols[gene],
      x = "Hours",
      y = bquote("Log"[2]~"TPM")
    ) +
    theme(
      panel.grid.major = element_line(color = 'grey60', size = 0.1),
      panel.grid.minor = element_blank()
      # panel.grid.major.x = element_blank(),
      # panel.grid.minor = element_line(color = 'grey50', size = 0.1),
      # panel.grid.minor.x = element_blank()
    )
  if (save_plot) {
    dir.create(dirname(filename), showWarnings = FALSE)
    ggsave_optimize_png(
      filename = sprintf(filename, gene_symbols[gene]),
      plot = p,
      # width = 14,
      # height = 8,
      width = 14,
      height = 6,
      units = "in"
    )
  }
  if (show_plot) {
    return(p)
  }
  return(NULL)
}

plot_gene <- function(
  gene = NULL,
  # trend = "lm",
  trend = FALSE,
  save_plot = TRUE,
  facet_by = "Disease ~ Stimulation",
  color_by = "Cell_Line",
  filename = "figures/genes/%s.pdf"
) {
  # Use the global variable that contains metadata.
  dat <- m
  
  # Try to accept Ensembl IDs or HGNC symbols.
  if (!grepl("^ENSG", gene)) {
    # Use the global variable with mappings from Ensembl to Symbol.
    gene <- names(which(gene_symbols == gene))
    stopifnot(length(gene) == 1)
  }
  dat[['Gene']] <- as.matrix(log2tpm)[gene,]
  
  # Add extra rows, only for the purpose of plotting the 0h time point.
  x <- dat[dat$Time == 0,]
  dat_zero <- rbind(x, x, x)
  stims <- unique(dat$Stimulation)
  stims <- stims[stims != "None"]
  dat_zero$Stimulation <- rep(stims, each = nrow(x))
  dat <- rbind(dat, dat_zero)
  dat <- subset(dat, Stimulation != "None")
  
  dat_mean <- aggregate(Gene ~ Disease + Time + Stimulation, dat, mean)
  
  p <- ggplot()
  
  if (trend == "lm") {
    p <- p + stat_smooth(
      data = dat,
      aes(x = Time, y = Gene, group = paste(Disease, Stimulation)),
      color = "grey50", method = "lm", size = 2
    )
  } else if (trend == "mean") {
    p <- p + geom_line(
      data = dat_mean,
      aes(x = Time, y = Gene, group = Disease),
      size = 1.5, color = "black"
    )
  }
  if (!is.null(facet_by)) {
    p <- p +
    # facet_wrap(Disease ~ Stimulation, ncol = 3) +
    # facet_grid(Disease ~ Stimulation) +
    facet_grid(as.formula(facet_by))
  }
  p <- p + geom_line(
    # data = subset(dat, Stimulation != "None"),
    data = dat,
    aes_string(
      x = "Time", y = "Gene",
      # color = Cell_Line, group = Cell_Line
      color = color_by, group = "Cell_Line"
    ),
    # shape = 21, size = 4
    size = 1
  ) +
    scale_color_manual(values = brewer.pal(8, "Dark2")) +
    scale_x_continuous(breaks = sort(unique(m$Time))) +
    # scale_y_continuous(labels = function(x) round(2 ^ x)) +
    # theme_cowplot(font_size = font_size) +
    theme_bw(base_size = font_size) +
    labs(
      title = gene_symbols[gene],
      x = "Hours",
      y = bquote("Log"[2]~"TPM")
    ) +
    theme(
      panel.grid.major = element_line(color = 'grey60', size = 0.1),
      panel.grid.minor = element_blank()
      # panel.grid.major.x = element_blank(),
      # panel.grid.minor = element_line(color = 'grey50', size = 0.1),
      # panel.grid.minor.x = element_blank()
    )
  
  if (save_plot) {
    dir.create(dirname(filename), showWarnings = FALSE)
    ggsave(
      filename = sprintf(filename, gene_symbols[gene]),
      plot = p,
      # width = 14,
      # height = 8,
      width = 12,
      height = 6,
      units = "in"
    )
  }
  
  return(p)
}

plot_gene_mean <- function(
  gene = NULL,
  save_plot = TRUE,
  color_by = "Disease",
  filename = "figures/mean-genes/%s.pdf"
) {
  # Use the global variable that contains metadata.
  dat <- m
  
  # Try to accept Ensembl IDs or HGNC symbols.
  if (!grepl("^ENSG", gene)) {
    # Use the global variable with mappings from Ensembl to Symbol.
    gene <- names(which(gene_symbols == gene))
    stopifnot(length(gene) == 1)
  }
  dat[['Gene']] <- as.matrix(log2tpm)[gene,]
  
  # Add extra rows, only for the purpose of plotting the 0h time point.
  x <- dat[dat$Time == 0,]
  dat_zero <- rbind(x, x, x)
  stims <- unique(dat$Stimulation)
  stims <- stims[stims != "None"]
  dat_zero$Stimulation <- rep(stims, each = nrow(x))
  dat <- rbind(dat, dat_zero)
  dat <- subset(dat, Stimulation != "None")
  
  # dat_mean <- aggregate(Gene ~ Disease + Time + Stimulation, dat, mean)
  
  dat_mean <- group_by(dat, Disease, Time, Stimulation) %>%
    summarize(
      Mean = median(Gene),
      # Mean = mean(Gene),
      SD = sd(Gene),
      Q5 = median(Gene) - mad(Gene),
      Q95 = median(Gene) + mad(Gene)
      # Q5 = quantile(Gene, 0.025),
      # Q95 = quantile(Gene, 0.975)
    )
  
  p <- ggplot()
  
  if (color_by == "Stimulation") {
    plot_width <- 10
    p <- p +
      geom_ribbon(
        data = dat_mean,
        mapping = aes(
          x = Time, ymin = Q5, ymax = Q95, fill = Stimulation
        ),
        alpha = 0.15
      ) +
      # geom_line(
      #   data = dat,
      #   mapping = aes(x = Time, y = Gene, Group = Cell_Line, color = Stimulation),
      #   size = 0.5, alpha = 0.5
      # ) +
      geom_line(
        data = dat_mean,
        mapping = aes(x = Time, y = Mean, color = Stimulation),
        size = 1
      ) +
      facet_grid(~ Disease) +
      # scale_color_manual(values = brewer.pal(8, "Set2")) +
      scale_color_manual(values = meta_colors$Stimulation) +
      scale_fill_manual(values = meta_colors$Stimulation) +
      scale_x_continuous(breaks = sort(unique(m$Time))) +
      # scale_y_continuous(labels = function(x) round(2 ^ x)) +
      # theme_cowplot(font_size = font_size) +
      theme_bw(base_size = font_size) +
      labs(
        title = gene_symbols[gene],
        x = "Hours",
        y = bquote("Log"[2]~"TPM")
      ) +
      theme(
        panel.grid.major = element_line(color = 'grey60', size = 0.1),
        panel.grid.minor = element_blank()
        # legend.position = c(0.15, 0.85)
      )
  } else {
    plot_width <- 15
    p <- p +
      geom_ribbon(
        data = dat_mean,
        mapping = aes(
          x = Time, ymin = Q5, ymax = Q95, fill = Disease
        ),
        alpha = 0.15
      ) +
      # geom_line(
      #   data = dat,
      #   mapping = aes(x = Time, y = Gene, Group = Cell_Line, color = Disease),
      #   size = 0.5, alpha = 0.5
      # ) +
      geom_line(
        data = dat_mean,
        mapping = aes(x = Time, y = Mean, color = Disease),
        size = 1
      ) +
      facet_grid(~ Stimulation) +
      scale_color_manual(values = brewer.pal(10, "Paired")[c(10, 8)]) +
      scale_fill_manual(values = brewer.pal(10, "Paired")[c(10, 8)]) +
      scale_x_continuous(breaks = sort(unique(m$Time))) +
      # scale_y_continuous(labels = function(x) round(2 ^ x)) +
      # theme_cowplot(font_size = font_size) +
      theme_bw(base_size = font_size) +
      labs(
        title = gene_symbols[gene],
        x = "Hours",
        y = bquote("Log"[2]~"TPM")
      ) +
      theme(
        panel.grid.major = element_line(color = 'grey60', size = 0.1),
        panel.grid.minor = element_blank(),
        legend.position = c(0.05, 0.85)
      )
  }
  
  if (save_plot) {
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
    ggsave(
      filename = sprintf(filename, gene_symbols[gene]),
      plot = p,
      width = plot_width,
      height = 6,
      units = "in"
    )
  }
  
  return(p)
}

plot_transcript <- function(
  transcript = NULL,
  trend = "lm",
  save_plot = TRUE,
  transcript_to_gene = NULL,
  filename = "figures/transcripts/%s.pdf"
) {
  # Drop the version number of the transcript: "ENST001.1" -> "ENST001"
  transcript <- drop_version(transcript)
  
  dat <- m
  dat[['Transcript']] <- as.matrix(log2tpmt)[transcript,]
  
  # Add extra rows, only for the purpose of plotting the 0h time point.
  x <- dat[dat$Time == 0,]
  dat_zero <- rbind(x, x, x)
  stims <- unique(dat$Stimulation)
  stims <- stims[stims != "None"]
  dat_zero$Stimulation <- rep(stims, each = nrow(x))
  dat <- rbind(dat, dat_zero)
  dat <- subset(dat, Stimulation != "None")
  
  dat_mean <- aggregate(Transcript ~ Disease + Time + Stimulation, dat, mean)
  
  Title <- transcript
  if (!is.null(transcript_to_gene)) {
    tg <- subset(t2g, target_id == transcript)
    if (nrow(tg) == 1) {
      Gene <- tg[["ext_gene"]]
      Title <- sprintf("%s\n%s", Gene, transcript)
    }
  }
  
  p <- ggplot()
  
  if (trend == "lm") {
    p <- p + stat_smooth(
      data = dat,
      aes(x = Time, y = Transcript, group = paste(Disease, Stimulation)),
      color = "grey50", method = "lm", size = 2
    )
  } else if (trend == "mean") {
    p <- p + geom_line(
      data = dat_mean,
      aes(x = Time, y = Transcript, group = Disease),
      size = 1.5, color = "black"
    )
  }
  p <- p + geom_line(
    # data = subset(dat, Stimulation != "None"),
    data = dat,
    aes(x = Time, y = Transcript, color = Cell_Line, group = Cell_Line),
    # shape = 21, size = 4
    size = 1
  ) +
    # facet_wrap(Disease ~ Stimulation, ncol = 3) +
    facet_grid(Disease ~ Stimulation) +
    scale_color_manual(values = brewer.pal(8, "Dark2")) +
    scale_x_continuous(breaks = sort(unique(m$Time))) +
    # scale_y_continuous(labels = function(x) round(2 ^ x)) +
    # theme_cowplot(font_size = font_size) +
    theme_bw(base_size = font_size) +
    labs(
      title = Title,
      x = "Hours",
      y = bquote("Log"[2]~"TPM")
    ) +
    theme(
      panel.grid.major = element_line(color = 'grey60', size = 0.1),
      panel.grid.minor = element_blank()
      # panel.grid.major.x = element_blank(),
      # panel.grid.minor = element_line(color = 'grey50', size = 0.1),
      # panel.grid.minor.x = element_blank()
    )
  
  if (save_plot) {
    dir.create(dirname(filename), showWarnings = FALSE)
    ggsave(
      filename = sprintf(filename, str_replace(Title, "\n", "_")),
      plot = p,
      width = 14,
      height = 8,
      units = "in"
    )
  }
  
  return(p)
}

plot_gene_by_disease <- function(gene = NULL, meta, log2tpm, ci = 0.5) {
  # Take the columns we need for this plot.
  dat <- meta[,c("Disease", "Stimulation", "Time")]
  # Try to accept Ensembl IDs or HGNC symbols.
  if (!grepl("^ENSG", gene)) {
    # Use the global variable with mappings from Ensembl to Symbol.
    gene <- names(which(gene_symbols == gene))
    stopifnot(length(gene) == 1)
  }
  dat[['Gene']] <- as.matrix(log2tpm)[gene,]
  # Add extra rows, only for the purpose of plotting the 0h time point.
  x <- dat[dat$Time == 0,]
  dat_zero <- rbind(x, x, x)
  stims <- unique(dat$Stimulation)
  stims <- stims[stims != "None"]
  dat_zero$Stimulation <- rep(stims, each = nrow(x))
  dat <- rbind(dat, dat_zero)
  dat <- subset(dat, Stimulation != "None")
  # Convert Time to a factor.
  dat$Time <- factor(dat$Time)
  dat_mean <- dat %>%
    group_by(Time, Disease) %>%
    summarise(
      median = median(Gene),
      low = quantile(Gene, probs = (1 - ci) / 2),
      high = quantile(Gene, probs = 1 - ((1 - ci) / 2))
    )
  ggplot(
    data = dat_mean,
    mapping = aes(
      x = Time,
      y = median,
      ymin = low,
      ymax = high,
      group = Disease,
      color = Disease
    )
  ) +
    geom_line(size = 1) + geom_ribbon(aes(fill = Disease), alpha = 0.15, color = NA) +
    scale_color_manual(values = meta_colors[["Disease"]]) +
    scale_fill_manual(values = meta_colors[["Disease"]]) +
    theme_bw(base_size = font_size) +
    labs(
      title = gene_symbols[gene],
      x = "Hours",
      y = bquote("Log"[2]~"TPM")
    ) +
    theme(
      panel.grid.major = element_line(color = 'grey60', size = 0.1),
      panel.grid.minor = element_blank()
    )
}

plot_gene_by_stimulation <- function (gene = NULL, meta, log2tpm, ci = 0.5, font_size = 14)
{
    dat <- meta[, c("Stimulation", "Time")]
    if (!grepl("^ENSG", gene)) {
        gene <- names(which(gene_symbols == gene))
        stopifnot(length(gene) == 1)
    }
    dat[["Gene"]] <- as.matrix(log2tpm)[gene, ]
    x <- dat[dat$Time == 0, ]
    dat_zero <- rbind(x, x, x)
    stims <- unique(dat$Stimulation)
    stims <- stims[stims != "None"]
    dat_zero$Stimulation <- rep(stims, each = nrow(x))
    dat <- rbind(dat, dat_zero)
    dat <- subset(dat, Stimulation != "None")
    #dat$Time <- factor(dat$Time)
    dat_mean <- dat %>%
      group_by(Time, Stimulation) %>%
      summarise(
        median = median(Gene),
        mid    = mean(Gene), 
        low    = mean(Gene) - sem(Gene),
        high   = mean(Gene) + sem(Gene)
#        low    = quantile(Gene, probs = (1 - ci) / 2),
#        high   = quantile(Gene, probs = 1 - ((1 - ci) / 2))
      )
    ggplot(
      data = dat_mean,
      mapping = aes(x = Time, y = mid, ymin = low, ymax = high, group = Stimulation, color = Stimulation)
    ) + 
      geom_ribbon(aes(fill = Stimulation), alpha = 0.15, color = NA) +
      scale_color_manual(values = brewer.pal(9, "YlOrRd")[c(4, 7, 9)]) +
      scale_fill_manual(values = brewer.pal(9, "YlOrRd")[c(4, 7, 9)]) +
      scale_x_continuous(breaks = unique(dat_mean$Time), labels = c(0, rep('', 7), 24)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
      geom_line(size = 1) +
      labs(
        title = gene_symbols[gene],
        x = "Hours",
        y = bquote("Log"[2] ~ "TPM")
      ) +
      theme_clean(base_size = font_size) + 
      theme(
        panel.border = element_rect(fill = NA, size = 0.2),
        axis.text = element_text(size = font_size),
        plot.title = element_text(face = "italic", size = font_size)
      )
}

# limma -----------------------------------------------------------------------

lmCompare <- function(object, fit1, fit2) {
  # Residuals
  res1 <- limma::residuals.MArrayLM(fit1, object)
  res2 <- limma::residuals.MArrayLM(fit2, object)
  # Residual sums of squares
  rss1 <- apply(res1, 1, function(xs) sum(xs ^ 2))
  rss2 <- apply(res2, 1, function(xs) sum(xs ^ 2))
  # Degrees of freedom
  df1 <- fit2$rank - fit1$rank
  df2 <- ncol(object) - fit2$rank
  # Number of parameters in each model
  if (df1 <= 0) {
    stop("fit2 must have more parameters than fit1")
  }
  # F statistics and p values
  #fs <- (rss1 - rss2) / (rss2 / df1)
  fs <- (rss1 - rss2) / rss2
  fs <- fs * df2 / df1
  ps <- pf(q = fs, df1 = df1, df2 = df2, lower.tail = FALSE)
  data.frame(
    "fstat" = fs,
    "pval"  = ps
  )
}

plot_limma_volcano <- function(
  dat, lfc = log2(1.5), fdr = 0.05, n_text = 10, classes = NULL,
  title = "", subtitle = "",
  fdr_line = TRUE
) {
  fdrs <- p.adjust(p = dat$P.Value, method = "fdr")
  p_threshold <- -log10(dat$P.Value[tail(which(fdrs < fdr), 1)])
  
  n_genes <- sum(abs(dat$logFC) >= lfc & -log10(dat$P.Value) > p_threshold)
  
  p <- ggplot(mapping = aes(x = logFC, y = -log10(P.Value))) +
    geom_point(
      data = dat,
      size = 0.5, color = "grey40"
    ) +
    geom_point(
      data = subset(dat, abs(logFC) >= lfc & adj.P.Val <= fdr),
      size = 1, color = "red"
    ) +
    scale_x_continuous(
                       breaks = scales::pretty_breaks(8),
      labels = function(x) fractional::fractional(2 ^ x)
    ) +
    geom_vline(
      xintercept = c(-lfc, 0, lfc),
      color = "grey70"
    ) +
    labs(
      x = bquote("Log"[2]~"Fold-Change"),
      y = bquote("-Log"[10]~"P"),
      title = title,
      subtitle = sprintf(
        "%s genes at %.2f absolute fold-change and %s%% FDR",
        comma(n_genes), 2^lfc, signif(fdr * 100, 1)
      )
    )
  
  if (fdr_line) {
    p <- p +
      geom_hline(
        yintercept = p_threshold,
        color = "grey70"
      ) +
      annotate(
        geom = "text",
        x = min(dat$logFC),
        y = p_threshold,
        hjust = 0,
        vjust = -0.3,
        color = "grey50",
        label = sprintf("%s%% FDR", signif(100 * fdr, 2))
      )
  }
  
  if (!is.null(classes)) {
    # %1 percent in the tails
    q50 <- quantile(dat$logFC, c(0.005, 0.995))
    p <- p + annotate(
      geom = "text",
      fontface = "bold",
      # x = log2(c(1 / (lfc * 1.2), lfc * 1.2)),
      # x = range(dat$logFC),
      x = range(dat$logFC[dat$logFC > q50[1] & dat$logFC < q50[2]]),
      y = max(-log10(dat$P.Value)) * 1.15,
      label = classes,
      size = 5
    )
  }
  
  if (n_text > 0) {
    # text_order <- with(
    #   dat, order(-log10(adj.P.Val) + abs(logFC), decreasing = TRUE))
    
    text_dat <- subset(dat, abs(logFC) >= lfc & adj.P.Val <= fdr)
    
    text_order <- with(
      text_dat, order(-log10(P.Value) + abs(logFC), decreasing = TRUE))
    
    p <- p + 
      geom_text_repel(
        data = head(text_dat[text_order, ], n_text),
        mapping = aes(x = logFC, y = -log10(P.Value), label = Gene),
        size = 4, fontface = "italic"
      ) +
      geom_point(
        data = head(text_dat[text_order, ], n_text),
        mapping = aes(x = logFC, y = -log10(P.Value)),
        shape = 21, size = 1.25, color = "black", fill = "red"
      )
  }
  
  return(p)
}

plot_limma_ma <- function(dat) {
  ggplot() +
    geom_point(
      data = dat,
      mapping = aes(x = AveExpr, y = logFC),
      size = 0.1
    ) +
    geom_point(
      data = subset(dat, abs(logFC) >= log2(1.5) & adj.P.Val <= 0.05),
      mapping = aes(x = AveExpr, y = logFC),
      size = 1, color = "red"
    ) +
    geom_text_repel(
      data = subset(dat, abs(logFC) >= log2(4) & adj.P.Val <= 0.05),
      mapping = aes(x = AveExpr, y = logFC, label = Gene),
      size = 3, color = "red"
    ) +
    geom_hline(
      yintercept = log2(c(1.5, 1, 0.6666)),
      color = "grey70"
    ) +
    scale_y_continuous(
      labels = function(x) signif(2^x, 2),
      breaks = pretty_breaks(5)
    ) +
    labs(
      x = bquote("Mean Log"[2]~"TPM"),
      y = bquote("Fold-Change")
    )
}
