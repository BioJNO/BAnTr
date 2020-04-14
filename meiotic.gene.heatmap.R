# Plot a heatmap of Meiotic genes in anther and meicyte transcript data

# Load required packages -----------------------------------------------------
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(plyr)
library(scales)

# Load and format data -------------------------------------------------------
expr.dat <- read.csv("anther_transcripts_pass_expression.csv", na.strings = c("", " "))
expr.dat <- read.csv("TranscriptomeResultsUpdate_2.csv", na.strings = c("", " "))
# mei.list <- read.csv("meiotic.gene.list.csv")

expr.dat.meiosis <- expr.dat[!is.na(expr.dat$MeiosisGenes), ]
expr.dat.meiosis$MeiosisGenes

colnames(expr.dat)

expr.dat.meiosis.sig <- subset(expr.dat.meiosis, A.LEP.vs.A.PRE.Pval < 0.01 & A.LEP.vs.A.PRE.log2FC > 1 |
                                 A.LEP.vs.A.PRE.Pval < 0.01 & A.LEP.vs.A.PRE.log2FC < -1 |
                                 A.PAC.vs.A.LEP.Pval < 0.01 & A.PAC.vs.A.LEP.log2FC > 1 |
                                 A.PAC.vs.A.LEP.Pval < 0.01 & A.PAC.vs.A.LEP.log2FC < -1 |
                                 A.TET.vs.A.PAC.Pval < 0.01 & A.TET.vs.A.PAC.log2FC > 1 |
                                 A.TET.vs.A.PAC.Pval < 0.01 & A.TET.vs.A.PAC.log2FC < -1 |
                                 M.LEP.vs.A.LEP.Pval < 0.01 & M.LEP.vs.A.LEP.log2FC > 1 |
                                 M.LEP.vs.A.LEP.Pval < 0.01 & M.LEP.vs.A.LEP.log2FC < -1 |
                                 M.PAC.vs.A.PAC.Pval < 0.01 & M.PAC.vs.A.PAC.log2FC > 1 |
                                 M.PAC.vs.A.PAC.Pval < 0.01 & M.PAC.vs.A.PAC.log2FC < -1 |
                                 M.PAC.vs.M.LEP.Pval < 0.01 & M.PAC.vs.M.LEP.log2FC > 1 |
                                 M.PAC.vs.M.LEP.Pval < 0.01 & M.PAC.vs.M.LEP.log2FC < -1)


# 29 significant DEG in any comparison
colnames(expr.dat.meiosis.sig)
# write.csv(expr.dat.meiosis.sig, "significant_meiotic_genes.csv")

expr.dat.meiosis.sig <- read.csv("significant_meiotic_genes.csv")

# order by wgcna module
# expr.dat.meiosis <- expr.dat.meiosis[order((expr.dat.meiosis$WGCNA)), ]
expr.dat.meiosis.sig <- expr.dat.meiosis.sig[order((expr.dat.meiosis.sig$WGCNA)), ]

# subset the non-significant genes
# unsig <- setdiff(expr.dat.meiosis$BAnTr, expr.dat.meiosis.sig$BAnTr)
# expr.dat.meiosis.unsig <- expr.dat.meiosis[expr.dat.meiosis$BAnTr %in% unsig, ]

# restrict above to data columns for plotting (gene name and transcript
# read numbers)
# colnames(expr.dat.meiosis)
# plot.data <- expr.dat.meiosis[, c(8,18:35,10)]

colnames(expr.dat.meiosis.sig)
plot.data <- expr.dat.meiosis.sig[, c(9,19:36,11)]

# colnames(expr.dat.meiosis.unsig)
# plot.data <- expr.dat.meiosis.unsig[, c(8,18:35,10)]

# scale each column (relative change in expression of each)
# This means columns can be compared but not rows.
# It will make all values very small!
colnames(plot.data)
# plot.data.scaled <- plot.data
#plot.data.scaled$total <- rowSums(plot.data.scaled[3:20])
#plot.data.scaled[,3:20] <- plot.data.scaled[,3:20]/plot.data.scaled$total
#plot.data.scaled <- plot.data.scaled[,-21]
#colnames(plot.data.scaled)

# re-order the colums to put anther and meiocyte adjacent
colnames(plot.data)
plot.data <- plot.data[,c(1:7,14:16,8:10,17:19,11:13,20)]

# convert from wide to long format.
long.plot.data <- melt(plot.data)

# Keep the genes in current order in the plot, not alphabetical.
long.plot.data$gene <- as.character(long.plot.data$MeiosisGenes)
long.plot.data$gene <- factor(long.plot.data$gene,
                              levels=unique(long.plot.data$gene))
# Check if the data is normally distributed.
# First plot a simple histogram and look at the distribution.
hist(long.plot.data$value)

# Clearly left skewed. Confirm with shapiro-wilk test,
# where p > 0.05 = normal distribution
shapiro.test(long.plot.data$value)

# Transform scaled values to normalise them, or get them closer to normal.
long.plot.data$log.expr <- log1p(long.plot.data$value)

# Check with a histogram again.
hist(long.plot.data$log.expr)
# Better
shapiro.test(long.plot.data$log.expr)

# Plot heatmap ---------------------------------------------------------------
long.plot.data$WGCNA
long.plot.data <- long.plot.data[order((long.plot.data$WGCNA)), ]
meiotic.heat <- ggplot(long.plot.data, aes(variable, gene)) +
                       geom_tile(aes(fill = log.expr), colour = "white") +
                       scale_fill_distiller(palette = "BrBG") +
                       theme(axis.text.x = element_text(angle = 45, hjust = 1),
                             axis.text.y = element_text(size = 6),
                             axis.title.x = element_blank(),
                             axis.title.y = element_blank()) +
                       labs(fill = "Log1p Expression")
meiotic.heat

ggsave("meiotic_heat_map_sig_3.svg", width=17.4, height=13.65, units = "cm", dpi=600)
