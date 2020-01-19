library(OSAT)
library(readxl)
library(stringr)

d <- data.frame(read_excel("data-raw/sample_concentrations.xlsx", sheet = 2))

d$Cell_Line <- str_replace_all(d$Cell_Line, "FLS ", "")

d <- d[d$Cell_Line != "OA453",]
d$Disease <- "OA"
d$Disease[str_detect(d$Cell_Line, "RA")] <- "RA"

d$ConcHigh <- d$Concentration >= median(d$Concentration)

gs <- setup.sample(
  x = d,
  optimal = c("Stimulation", "Disease", "Cell_Line", "ConcHigh"),
  strata = "Time"
)

myChip <- new("BeadChip", nRows=8, nColumns=12, byrow=TRUE)
myPlate <- new("BeadPlate", chip=myChip, nRows=1L, nColumns=1L)

# gc <- setup.container(IlluminaBeadChip96Plate, 2, batch="plates")
gc <- setup.container(myPlate, 2, batch="plates")

set.seed(1234)
gSetup <- create.optimized.setup(
  fun = "optimal.block",
  sample = gs,
  container = gc,
  nSim = 1e3
)
QC(gSetup)

pdf("figures/plate_design.pdf", width = 14, height = 8)
plot(gSetup)
dev.off()

out <- map.to.MSA(gSetup, MSA4.plate)
out$plates <- NULL
out$wellID <- NULL

write.csv(
  x = out,
  file = "data-raw/plate_design.csv",
  row.names = FALSE
)