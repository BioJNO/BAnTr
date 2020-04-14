# Join all data ---------------------------------------------------------------
expr.dat <- read.csv("TranscriptomeResultsUpdate_2.csv")
full_nog <- read.csv("full_nog.csv")
GOs <- read.csv("PANNZER_GOs.csv")

expr.dat <- merge(expr.dat, full_nog, by="BAnTr", all.x=TRUE)
expr.dat <- merge(expr.dat, GOs, by="BAnTr", all.x=TRUE)
expr.dat <- subset(expr.dat, ExpressionFilter != "LowExpressionFilter")

# depreciated part of script kept for example of removing duplicates
# UBIQ.all <- UBIQ.all[!duplicated(UBIQ.all$BAnTr), ]

# Grep ligases  based on GO terms ---------------------------------------------
# GO:0061630 = ubiquitin protein ligase activity
# GO:0004842 = ubiquitin-protein transferase activity
# GO:0016567 = protein ubiquitination
# GO:0034450 = ubiquitin-ubiquitin ligase activity
# GO:0000151 = ubiquitin ligase complex
# GO:0031462 = Cul2-RING ubiquitin ligase complex
# GO:0019005 = SCF ubiquitin ligase complex
# GO:0031461 = cullin-RING ubiquitin ligase complex

# GO:0004839 = ubiquitin activating enzyme activity
# GO:0019948 = SUMO activiating enzyme activity
# GO:0031510 = SUMO activating enzyme complex

# GO:0061656 = SUMO conjugating enzyme activity
# GO:0061631 = ubiquitin conjugating enzyme activity
# GO:0031371 = ubiquitin conjugating enzyme complex
# GO:0061650 = ubiquitin-like protein conjugating enzyme activity

# GO:0019788 = NEDD8 transferase activity
# GO:0016579 = protein deubiquitination
# GO:0000338 = protein deneddylation
# GO:0016926 = protein desumoylation

# GO:0061659 = ubiquitin-like protein ligase activity
# GO:0019787 = ubiquitin-like protein transferase activity
# GO:0019789 = SUMO transferase activity
# GO:0106068 = SUMO ligase complex 
# GO:0061665 = SUMO ligase activity
# GO:0140082 = SUMO-ubiquitin ligase activity
# GO:0033768 = SUMO-targeted ubiquitin ligase complex

# Without including ligase GO term we were not selecting a number of genes
# described as SINA ligases but without Ubiquitin E3 ligase GO terms! 
# This is how this problem was identified.
# 
# SINA_grep_desc <- UBIQ.all[grepl("SINA", UBIQ.all$PANNZER) |
#                            grepl("SINA", UBIQ.all$Mercator.annotation) |
#                            grepl("SINA", UBIQ.all$Mercator.Bin.Name) |
#                            grepl("SINA", UBIQ.all$eggNOG), ]
# 
# SINA_grep_GO <- go_ligase[grepl("SINA", go_ligase$PANNZER) |
#                           grepl("SINA", go_ligase$Mercator.annotation) |
#                           grepl("SINA", go_ligase$Mercator.Bin.Name) |
#                           grepl("SINA", go_ligase$eggNOG), ]
# 
# desc_not_GO <- setdiff(SINA_grep_desc$BAnTr, SINA_grep_GO$BAnTr)
# GO_not_desc <- setdiff(SINA_grep_GO$BAnTr, SINA_grep_desc$BAnTr)
# 
# missing_SINA <- SINA_grep_desc[SINA_grep_desc$BAnTr %in% desc_not_GO, ]
# write.csv(missing_SINA, "described_as_SINA_not_annotated_GO.csv")

ligase.terms <- c("GO:0061630", "GO:0004842", "GO:0000151",
                  "GO:0031462", "GO:0019005", "GO:0031461",
                  "GO:0019788", "GO:0061659", "GO:0019787",
                  "GO:0019789", "GO:0106068", "GO:0061665",
                  "GO:0140082", "GO:0033768", "GO:0016874")

go_ligase <- expr.dat[grepl(paste(ligase.terms, collapse = "|"),
                          expr.dat$GO_PANNZER), ]

# While GO term GO:0016874 pulls out E3 ligases that would be otherwise
# excluded by ubiquitin ligase specific terms it also pulls out some
# genes that I don't want.
# Identify the terms we don't want 
e3.ligase.terms <- c("GO:0061630", "GO:0004842", "GO:0000151",
                     "GO:0031462", "GO:0019005", "GO:0031461",
                     "GO:0019788", "GO:0061659", "GO:0019787",
                     "GO:0019789", "GO:0106068", "GO:0061665",
                     "GO:0140082", "GO:0033768")
e3_go_ligase <- expr.dat[grepl(paste(e3.ligase.terms, collapse = "|"),
                            expr.dat$GO_PANNZER), ]

extra_with_ligase_list <- setdiff(go_ligase$BAnTr, e3_go_ligase$BAnTr)
ligase_not_ub_GO <- go_ligase[go_ligase$BAnTr %in% extra_with_ligase_list, ]
# An extra 503 genes. Some of which are not ubiqitin or ubiquitin like ligases.
# Care is needed here as some genes with unexpected descriptions (for example
# annotated as histone methyl-transferases) may be capable of E3 ligase
# activity in addition to their canonical function.

# Exclude ligase GO terms extras which are explicily described as E3 ligases
e3.descriptors <- c("SKP", "Cullin", "F-box", "SCF",
                    "E3", "APC/C", "Anaphase promoting complex",
                    "SINA")
ligase_extra_ub <- ligase_not_ub_GO[grepl(paste(e3.descriptors, collapse = "|"),
                                         ligase_not_ub_GO$PANNZER,
                                         ignore.case = TRUE) |
                                  grepl(paste(e3.descriptors, collapse = "|"),
                                         ligase_not_ub_GO$Mercator.annotation,
                                         ignore.case = TRUE) |
                                  grepl(paste(e3.descriptors, collapse = "|"),
                                         ligase_not_ub_GO$eggNOG,
                                         ignore.case = TRUE), ]

all_ligase <- rbind(e3_go_ligase, ligase_extra_ub)

# subset tissue specific ligases
anther_ligase <- subset(all_ligase, ExpressionFilter == "AntherDataset")
embryo_ligase <- subset(all_ligase, ExpressionFilter == "EmbryoDataset")
both_ligase <- subset(all_ligase, ExpressionFilter == "BothDatasets")

# Subset stat-sig DEG ligases
ligase_anther <- subset(all_ligase, ExpressionFilter != "EmbryoDataset")
# write out and read in to stop NAs from embryo data reading columns as factors
write.csv(ligase_anther, "all_ligase_anther.csv")
ligase_anther <- read.csv("all_ligase_anther.csv")

# All stat-sig in all contrase groups.
ligase_sig <- subset(ligase_anther,
                    A.LEP.vs.A.PRE.Pval <= 0.01 & A.LEP.vs.A.PRE.log2FC >= 1 |
                    A.LEP.vs.A.PRE.Pval <= 0.01 & A.LEP.vs.A.PRE.log2FC <= -1 |
                    A.PAC.vs.A.LEP.Pval <= 0.01 & A.PAC.vs.A.LEP.log2FC >= 1 |
                    A.PAC.vs.A.LEP.Pval <= 0.01 & A.PAC.vs.A.LEP.log2FC <= -1 |
                    A.TET.vs.A.PAC.Pval <= 0.01 & A.TET.vs.A.PAC.log2FC >= 1 |
                    A.TET.vs.A.PAC.Pval <= 0.01 & A.TET.vs.A.PAC.log2FC <= -1 |
                    M.LEP.vs.A.LEP.Pval <= 0.01 & M.LEP.vs.A.LEP.log2FC >= 1 |
                    M.LEP.vs.A.LEP.Pval <= 0.01 & M.LEP.vs.A.LEP.log2FC <= -1 |
                    M.PAC.vs.A.PAC.Pval <= 0.01 & M.PAC.vs.A.PAC.log2FC >= 1 |
                    M.PAC.vs.A.PAC.Pval <= 0.01 & M.PAC.vs.A.PAC.log2FC <= -1 |
                    M.PAC.vs.M.LEP.Pval <= 0.01 & M.PAC.vs.M.LEP.log2FC >= 1 |
                    M.PAC.vs.M.LEP.Pval <= 0.01 & M.PAC.vs.M.LEP.log2FC <= -1)
write.csv(ligase_sig, "Significant_ligase.csv")

# Subset ligases enriched in meiocytes compared to anthers.
meiocyte_enriched <- subset(ligase_sig,
                            M.LEP.vs.A.LEP.Pval <= 0.01 & M.LEP.vs.A.LEP.log2FC >= 1 |
                            M.PAC.vs.A.PAC.Pval <= 0.01 & M.PAC.vs.A.PAC.log2FC >= 1)
write.csv(meiocyte_enriched, "meiocyte_enriched_ligases.csv")

# Heatmap ---------------------------------------------------------------------
library(ggplot2)
library(RColorBrewer)
library(reshape)
library(plyr)
library(scales)

# restrict above to data columns for plotting (BAnTr and transcript read numbers)
colnames(meiocyte_enriched)
plot.data <- meiocyte_enriched[, c(2,18:39)]

# re-order the colums to put anther and meiocyte adjacent
colnames(plot.data)
plot.data <- plot.data[,c(1:7,14:16,8:10,17:19,11:13, 20:23)]

# convert from wide to long format.
long.plot.data <- melt(plot.data)

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
# Better
shapiro.test(long.plot.data$log.expr)

# Plot heatmap ---------------------------------------------------------------
merged <- merge(long.plot.data, meiocyte_enriched, by="BAnTr")

merged$WGCNA
merged <- merged[order((merged$WGCNA)), ]
meiotic.heat <- ggplot(merged, aes(variable, BAnTr)) +
  geom_tile(aes(fill = log.expr), colour = "white") +
  scale_fill_distiller(palette = "BrBG") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(fill = "Log10 Expression")
meiotic.heat

ggsave("E3_ligase_meiocyte_enriched_heat.svg", width=17.4, height=22, units = "cm", dpi=600)
