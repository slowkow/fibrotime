# Load libraries -------------------------------------

pacman::p_load(
  data.table,
  gsheet,
  janitor,
  dplyr,
  ggplot2,
  patchwork,
  ggstance,
  fractional,
  RColorBrewer,
  stringr,
  limma,
  tidyr
)

cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9",
  "#009E73", "#F0E442", "#0072B2",
  "#D55E00", "#CC79A7"
)

theme_set(theme_classic(base_size = 18) %+replace% theme(
  strip.background = element_blank(),
  # axis.line.y = element_line(colour = "black", size = 0.2),
  # axis.line.x = element_line(colour = "black", size = 0.2),
  axis.ticks   = element_line(colour = "black", size = 0.3),
  panel.border = element_rect(size = 0.3, fill = NA),
  axis.line    = element_blank(),
  plot.title   = element_text(size = 18, vjust = 2, hjust = 0.5),
  panel.grid.major = element_line(size = 0.2)
))

# Read data  -------------------------------------

chip <- read_excel("data/experiments.xlsx", sheet = "ChIP-QPCR")
chip <- clean_names(chip)
chip$date <- as.character(chip$date)
chip$sample <- sprintf("S%d", seq(nrow(chip)))
chip$time <- as.character(chip$time)
# chip$row <- substr(chip$well, 1, 1)
# chip$replicate <- rep(c(1, 2), each = 6)
head(chip, 6)

chip <- subset(chip, date %in% c("2019-10-11", "2019-10-12"))

chip <- chip[,apply(chip, 2, function(x) length(unique(x))) > 1]

table(chip$promoter, chip$date)

chip$measurement <- 1:2

compute_chip <- function(x) {
  x1 <- x[x$antibody != "Input",]
  x0 <- x[x$antibody == "Input",]
  stopifnot(all(x1$promoter == x0$promoter))
  stopifnot(all(x1$time == x0$time))
  x1$ct_deltarn <- x1$ct_delta_rn - x0$ct_delta_rn
  x1$ddct <- 2 ^ - x1$ct_delta_rn
  x1$promoter <- str_split_fixed(x1$promoter, ",", 2)[,1]
  x1
}

my_chip <- rbind(
  compute_chip(chip[chip$date == "2019-10-11",]),
  compute_chip(chip[chip$date == "2019-10-12",])
) %>%
  group_by(antibody, promoter, time, replicate) %>%
  summarize(
    ct_delta_rn = mean(ct_delta_rn),
    ddct = mean(ddct)
  )


# Plot data  -------------------------------------

my_chip$sample <- with(my_chip, order(antibody, promoter, time, replicate))

p1 <- ggplot(my_chip) +
  aes(
    x = time,
    y = 100 * ddct,
    # group = paste(promoter, replicate, antibody),
    group = sample,
    fill = antibody
  ) +
  geom_point(
    position = position_dodge(width = 1), shape = 21, size = 3
  ) +
  facet_grid( ~ promoter) +
  scale_y_log10(breaks = c(1e-4, 1e-5, 1e-6)) +
  annotation_logticks(
    base = 10, sides = "l", size = 0.3
  ) +
  scale_fill_manual(
    name = "Antibody", values = c("white", "grey50")
  ) +
  theme(
    strip.text = element_text(face = "italic", size = 16),
    plot.title = element_text(size = 16),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.position = "bottom",
    legend.text = element_text(margin = margin(r = 30, unit = "pt"))
  ) +
  labs(
    x = "Hours after TNF and IL-17A",
    y = bquote("ddC"[t])
  )

x <- my_chip %>%
  group_by(antibody, promoter, replicate) %>%
  # group_by(antibody, measurement, promoter, replicate) %>%
  summarize(
    t1 = ddct[time == 1] / ddct[time == 0]
  )
#x <- x %>% gather("time", "fold", -promoter, -replicate, -measurement, -antibody)

x %>% group_by(antibody, promoter) %>%
  summarize(fold = mean(t1), foldsd = sd(t1)) %>%
  group_by(antibody) %>% summarize(mean = mean(fold))

p2 <- x %>% group_by(antibody, promoter) %>%
  summarize(fold = mean(t1), foldsd = sd(t1)) %>%
ggplot() +
  aes(x = antibody, fill = antibody) +
  facet_grid(~ promoter) +
  geom_hline(yintercept = 1, linetype = 2, size = 0.5) +
  geom_pointrange(
    aes(y = fold, ymin = fold - foldsd, ymax = fold + foldsd)
  ) +
  geom_point(aes(y = fold), shape = 21, size = 3) +
  scale_fill_manual(
    name = "Antibody", values = c("white", "grey50"), guide = FALSE
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(2)) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    strip.text = element_text(face = "italic", size = 16),
    plot.title = element_text(size = 16),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16)
  ) +
  labs(y = "Fold recruitment", x = "Antibody")

p2 +
  (p1 + theme(strip.text = element_blank())) +
  plot_layout(ncol = 1, heights = c(1, 1.5))

ggsave("chip2.pdf", width = 5, height = 6, units = "in", dpi = 300, useDingbats = FALSE)
