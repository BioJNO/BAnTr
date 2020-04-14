# Plot a heatmap of Meiotic genes in anther and meicyte transcript data

# Load required packages -----------------------------------------------------
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(plyr)
library(scales)

# Load and format data -------------------------------------------------------
expr.dat1 <- read.csv("upreg_meiocyte_pachytene_TFs.csv")
expr.dat2 <- read.csv("upreg_meiocyte_leptotene_TFs.csv")

expr.dat <- rbind(expr.dat1, expr.dat2)

# scale each column (relative change in expression of each)
# This means columns can be compared but not rows.
# It will make all values very small!
colnames(expr.dat)
expr.dat$total <- rowSums(expr.dat[5:22])
expr.dat[,5:22] <- expr.dat[,5:22]/expr.dat$total
colnames(expr.dat)

plot.data <- expr.dat[,c(2,5:22)]
colnames(plot.data)

# re-order the colums to put anther and meiocyte adjacent
plot.data <- plot.data[,c(1:7,14:16,8:10,17:19, 11:13)]

# convert to long format
m.plot.data <- melt(plot.data)
m.plot.data$gene.id <- as.character(m.plot.data$gene.id)
m.plot.data$gene.id <- factor(m.plot.data$gene.id, levels=unique(m.plot.data$gene))
# As scaling the rows makes all values small plotting the raw scaled values
# will give a very blue map. Log transformation makes the relative changes 
# visible when plotted but values will be negative.
m.plot.data$logexpr <- log(m.plot.data$value)

# Plot heatmap ---------------------------------------------------------------
TF.heat <- ggplot(m.plot.data, aes(variable, gene.id)) +
  geom_tile(aes(fill = logexpr), colour = "white") +
  scale_fill_distiller(palette = "BrBG") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
TF.heat

ggsave("meiocyte_TFs_upreg_heat_map.svg", height = 12, width = 9, dpi=1200)
