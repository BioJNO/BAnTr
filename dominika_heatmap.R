# Load packages ---------------------------------------------------------------
# Data transformation.
library(reshape2)
# Generate heatmap.
library(ggplot2)
library(RColorBrewer)

# Extract data ----------------------------------------------------------------
# Create a list of genes to plot.
dom_mikado <- c("mikado.chr1_part2G6880", "mikado.chr7_part1G1000",
                "mikado.chr3_part1G4610", "mikado.chr2_part1G3270",
                "mikado.chr3_part1G2852", "mikado.chr2_part2G3702",
                "mikado.chr2_part2G3978", "mikado.chr5_part2G8170",
                "mikado.chr6_part2G1532", "mikado.chr2_part2G4100",
                "mikado.chr2_part2G4110", "mikado.chr2_part2G4112",
                "mikado.chr5_part2G792", "mikado.chr5_part2G8148",
                "mikado.chr7_part2G1774", "mikado.chr7_part2G6810",
                "mikado.chr1_part2G324", "mikado.chr2_part2G76",
                "mikado.chr5_part2G700", "mikado.chr5_part1G2702",
                "mikado.chr2_part1G1990", "mikado.chr7_part2G996",
                "mikado.chr1_part1G2354", "mikado.chr4_part2G4084",
                "mikado.chr6_part2G1536", "mikado.chr7_part1G6482",
                "mikado.chr5_part2G8548", "mikado.chr1_part2G442",
                "mikado.chr7_part2G1636", "mikado.chr5_part2G854")

# Load in the expression data.
expr.dat <- read.csv("GeneExpressionR.csv")
# Pull out rows with the genes you want to plot.
plot.data <- expr.dat[expr.dat$gene.id %in% dom_mikado, ]
## Sort them by wgcna module, roughly the expression pattern.
#plot.data <- plot.data[order((plot.data$wgcna.mod)), ]
# Check to see that the list and extracted genes match,
# hope to see character(0).
# check that we didn't get anything we didn't want
check <- setdiff(plot.data$gene.id, dom_mikado)
check
# Check that we got everything we did want
check <- setdiff(dom_mikado, plot.data$gene.id)
check
# Re-order data columns so that gene and module data is at the start and the
# columns are in order of the meiotic stage of tissues:
# premeiosis -> leptotene/zygotene -> pachytene/diplotene -> anaphase/tetrad.
colnames(plot.data)
plot.data <- plot.data[,c(1,21,3:8,15:17,9:11,18:20,12:14)]

# Data transformation ---------------------------------------------------------
# Vastly different expression levels makes the relative change acrosss a set
# difficult to see.
plot.data.scaled <- plot.data

## Create a column with the sum of transcript count for each gene in all samples.
#plot.data.scaled$total <- rowSums(plot.data.scaled[3:20])

## Convert the raw count values into values relative to the total for each gene.
## The proper name for this is "feature scaling". 
#plot.data.scaled[,3:20] <- plot.data.scaled[,3:20]/plot.data.scaled$total

## Remove the total column.
#plot.data.scaled <- plot.data.scaled[,-21]
#colnames(plot.data.scaled)

# convert from wide to long format.
long.plot.data <- melt(plot.data.scaled)

# Keep the genes in current order in the plot, not alphabetical.
long.plot.data$gene.id <- as.character(long.plot.data$gene.id)
long.plot.data$gene.id <- factor(long.plot.data$gene,
                                 levels=unique(dom_mikado))

# Check if the data is normally distributed.
# First plot a simple histogram and look at the distribution.
hist(long.plot.data$value)

# Clearly left skewed. Confirm with shapiro-wilk test,
# where p > 0.05 = normal distribution
shapiro.test(long.plot.data$value)

# Transform scaled values to normalise them, or get them closer to normal.
long.plot.data$log.expr <- log10(long.plot.data$value)

# Check with a histogram again.
hist(long.plot.data$log.expr)
# Better but still skewed. 

# Generate heatmaps -----------------------------------------------------------
# Using RColorBrewer pallette. Diverging pallettes give an appreciable scale
# from a "cold" to a "hot" colour.
#
# Some diverging colour blind friendly pallettes are: BrBG, RdYlBu, RdBu, PuOr,
# PRGn, and PiYG. Gone with RdYlBu because red and blue are "heat-y" colours.
meiotic.heat <- ggplot(long.plot.data, aes(variable, gene.id)) +
  geom_tile(aes(fill = log.expr), colour = "white") +
  scale_fill_distiller(palette = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
meiotic.heat

# Save it.
ggsave("Dominika_heat_map_log_only.svg",
       height = 5, width = 6, dpi=1200)
