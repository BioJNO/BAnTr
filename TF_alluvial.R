# expr.dat <- read.csv("expression_dat_with_annotation.csv")
# 
# # Exclude ncRNAs 
# colnames(expr.dat)
# print(levels(expr.dat$coding.ncRNA))
# coding <- subset(expr.dat, coding.ncRNA =="coding")
# 
# # replace 0 and N/A values
# print(levels(coding$COG))
# levels(coding$COG) <- c(levels(coding$COG), "None") 
# coding$COG[coding$COG=="#N/A"]  <- "None"
# coding$COG[coding$COG=="N/A"]  <- "None"
# coding$COG[coding$COG=="0"] <- "None"
# 
# 
# all_sig <- subset(coding,
#                   alep.mlep.adj.P.Val < 0.01 & alep.mlep.logFC > 1 |
#                     alep.mlep.adj.P.Val < 0.01 & alep.mlep.logFC < -1 |
#                     apac.mpac.adj.P.Val < 0.01 & apac.mpac.logFC > 1 |
#                     apac.mpac.adj.P.Val < 0.01 & apac.mpac.logFC < -1 |
#                     lep.pre.adj.P.Val < 0.01 & lep.pre.logFC > 1 |
#                     lep.pre.adj.P.Val < 0.01 & lep.pre.logFC < -1 |
#                     pac.lep.adj.P.Val < 0.01 & pac.lep.logFC > 1 |
#                     pac.lep.adj.P.Val < 0.01 & pac.lep.logFC < -1 |
#                     tet.pac.adj.P.Val < 0.01 & tet.pac.logFC > 1 |
#                     tet.pac.adj.P.Val < 0.01 & tet.pac.logFC < -1)
# 
# all_sig_up <- subset(coding,
#                   alep.mlep.adj.P.Val < 0.01 & alep.mlep.logFC > 1 |
#                     apac.mpac.adj.P.Val < 0.01 & apac.mpac.logFC > 1 |
#                     lep.pre.adj.P.Val < 0.01 & lep.pre.logFC > 1 |
#                     pac.lep.adj.P.Val < 0.01 & pac.lep.logFC > 1 |
#                     tet.pac.adj.P.Val < 0.01 & tet.pac.logFC > 1)
# all_sig_up_O <- subset(all_sig, COG == "O")
# 
# write.csv(all_sig_up_O, "ALL_upregulated_COG_O.csv")
# 
# all_sig_up_O_Ub_mercator <- all_sig_up_O[grep("ubiquitin", all_sig_up_O$Mercator.Bin.Name), ]
# all_sig_up_O_Ub_nogg <- all_sig_up_O[grep("ubiquitin", all_sig_up_O$eggNOG.Annotation), ]
# all_sig_up_O_Ub_panzer <- all_sig_up_O[grep("ubiquitin", all_sig_up_O$PANNZER.annotation), ]
# 
# all_Ub_up_sig <- rbind(all_sig_up_O_Ub_mercator, all_sig_up_O_Ub_nogg, all_sig_up_O_Ub_panzer)
# 
# write.csv(all_Ub_up_sig, "all_Ub_upreg_any_stage.csv")
# 
# 
# # e.g. pull out genes annotated as Transcription factors
#tfs.total <- coding[grepl("transcription factor",
#                                    coding$Mercator.annotation,
#                                    ignore.case = TRUE), ]
# # 1353 total
# 
# tfs.sig <- all_sig[grepl("transcription factor",
#                           all_sig$Mercator.annotation,
#                           ignore.case = TRUE), ]
# # 382 significant
# 
# # Add a new column based on the values in another ----------------------------
# # e.g. for genes annotated as transcription factors create a tf.family column
# # based on the annotation in the mercator.Function column
# #
# # Mercator annotations follow the format *.FAMILY transcription factor
# # e.g. "Cell cycle.interphase.G1 phase.DP transcription factor"
# #
# # we need a package called stringr for this 
# library(stringr)
# 
# # Make a list of plant TF families 
# tfs <- c("AP2", "ARF", "ARR-B", "B3", "BBR-BPC", "BES1", "C2H2", "C3H",
#          "CAMTA", "CO-like", "CPP", "DBB", "Dof", "E2F/DP", "EIL",
#          "ERF", "FAR1", "G2-like", "GATA", "GRAS", "GRF", "GeBP", "HB", "PHD",
#          "HB-other", "HD-ZIP", "HRT-like", "HSF", "LBD", "LFY", "LSD",
#          "MADS", "MYB", "MYB_related", "NAC", "NF-X1", "NF-YA", "NF-YB",
#          "NF-YC", "Nin-like", "RAV", "S1Fa-like", "SBP", "SRS", "STAT",
#          "TALE", "TCP", "Trihelix", "VOZ", "WOX", "WRKY", "Whirly", "YABBY",
#          "zf-HD", "bHLH", "bZIP")
# regex = paste(tfs, collapse = "|")
# 
# tfs.total$tf.family <- sapply(str_extract_all(tfs.total$Mercator.annotation, regex),
#                               function(x) paste(x, collapse = ";"))
# 
# # Look at the column to check this has worked
# tfs.total$tf.family
# colnames(tfs.total)
# 
# 
# tfs.sig$tf.family <- sapply(str_extract_all(tfs.sig$Mercator.annotation, regex),
#                               function(x) paste(x, collapse = ";"))
# 
# # Look at the column to check this has worked
# tfs.sig$tf.family
# colnames(tfs.total)
# 
# #write.csv(tfs.sig, "significant_TFs.csv")
# Some manual cleanup of TF.family columns including TFs not included in list

sig.tfs <- read.csv("significant_TFs.csv")
lep_vs_pre_up <- subset(sig.tfs,
                        lep.pre.adj.P.Val < 0.01 &
                          lep.pre.logFC > 1)
lep_vs_pre_up$comp <- "A.LepZyg-A.Pre"
lep_vs_pre_up$change <- "Up-regulated"

lep_vs_pre_down <- subset(sig.tfs,
                          lep.pre.adj.P.Val < 0.01 &
                            lep.pre.logFC < -1)
lep_vs_pre_down$comp <- "A.LepZyg-A.Pre"
lep_vs_pre_down$change <- "Down-regulated"

pac_vs_lep_up <- subset(sig.tfs,
                        pac.lep.adj.P.Val < 0.01 &
                          pac.lep.logFC > 1)
pac_vs_lep_up$comp <- "A.PacDip-A.LepZyg"
pac_vs_lep_up$change <- "Up-regulated"

pac_vs_lep_down <- subset(sig.tfs,
                          pac.lep.adj.P.Val < 0.01 &
                            pac.lep.logFC < -1)
pac_vs_lep_down$comp <- "A.PacDip-A.LepZyg"
pac_vs_lep_down$change <- "Down-regulated"

tet_vs_pac_up <- subset(sig.tfs,
                        tet.pac.adj.P.Val < 0.01 &
                          tet.pac.logFC > 1)
tet_vs_pac_up$comp <- "A.MetTet-A.PacDip"
tet_vs_pac_up$change <- "Up-regulated"

tet_vs_pac_down <- subset(sig.tfs,
                          tet.pac.adj.P.Val < 0.01 &
                          tet.pac.logFC < -1)
tet_vs_pac_down$comp <- "A.MetTet-A.PacDip"
tet_vs_pac_down$change <- "Down-regulated"


lep_vs_pre_sig <- subset(sig.tfs,
                         lep.pre.adj.P.Val < 0.01 &
                           lep.pre.logFC > 1 |
                           lep.pre.adj.P.Val < 0.01 &
                           lep.pre.logFC < -1)
colnames(sig.tfs)

GENEID <- setdiff(sig.tfs$GENEID, lep_vs_pre_sig$GENEID)
tf.family <- sig.tfs[grep(paste(GENEID, collapse = "|"), sig.tfs$GENEID), 44]
comp <- rep("A.LepZyg-A.Pre", length(GENEID))
change <- rep("Not significant", length(GENEID))

lep_vs_pre_ns <- data.frame(GENEID, comp, change, tf.family)

# Count differentially expressed transcript COGs between anthers at pachytene
# and leptotene
pac_vs_lep_sig <- subset(sig.tfs,
                         pac.lep.adj.P.Val < 0.01 &
                           pac.lep.logFC > 1 |
                           pac.lep.adj.P.Val < 0.01 &
                           pac.lep.logFC < -1)

GENEID <- setdiff(sig.tfs$GENEID, pac_vs_lep_sig$GENEID)
comp <- rep("A.PacDip-A.LepZyg", length(GENEID))
change <- rep("Not significant", length(GENEID))
tf.family <- sig.tfs[grep(paste(GENEID, collapse = "|"), sig.tfs$GENEID), 44]

pac_vs_lep_ns <- data.frame(GENEID, comp, change, tf.family)

# Count differentially expressed transcript COGs between anthers at anaphase
# and pachytene
tet_vs_pac_sig <- subset(sig.tfs,
                         tet.pac.adj.P.Val < 0.01 &
                           tet.pac.logFC > 1 |
                           tet.pac.adj.P.Val < 0.01 &
                           tet.pac.logFC < -1)

GENEID <- setdiff(sig.tfs$GENEID, tet_vs_pac_sig$GENEID)
comp <- rep("A.MetTet-A.PacDip", length(GENEID))
change <- rep("Not significant", length(GENEID))
tf.family <- sig.tfs[grep(paste(GENEID, collapse = "|"), sig.tfs$GENEID), 44]

tet_vs_pac_ns <- data.frame(GENEID, comp, change,tf.family)

all_sig_comb <- rbind(tet_vs_pac_down,
                      tet_vs_pac_up,
                      pac_vs_lep_down,
                      pac_vs_lep_up,
                      lep_vs_pre_down,
                      lep_vs_pre_up)

all_sig_comb$comp <- as.factor(all_sig_comb$comp)
print(levels(all_sig_comb$comp))
all_sig_comb$comp <- factor(all_sig_comb$comp,
                            levels(all_sig_comb$comp)[c(1,3,2)])
colnames(all_sig_comb)

plot.data <- all_sig_comb[,c(3,45,46,44)]

plot.data <- rbind(plot.data,
                   lep_vs_pre_ns,
                   pac_vs_lep_ns,
                   tet_vs_pac_ns)
library(ggalluvial)


cbPalette <- c("#E69F00",
               "#56B4E9",
               "#CC79A7")


plot.data$change <- as.factor(plot.data$change)
print(levels(plot.data$change))
plot.data$change <- factor(plot.data$change,
                           levels(plot.data$change)[c(3,2,1)])
colnames(all_sig_comb)

ggplot(plot.data,
       aes(x = comp, stratum = change, alluvium = GENEID,
           fill = change, label = change, color=change)) +
  scale_fill_manual(values = cbPalette) +
  scale_color_manual(values = cbPalette) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback") +
  geom_stratum(alpha=0.5) +
  theme(legend.position = "bottom") +
  xlab("Tissue Comparision") +
  ylab("Number of TFs") +
  theme(text = element_text(size=16))

ggsave("Sig.tfs.by.comparison.change.svg", width=17.4, units="cm", dpi=600)

ggplot(plot.data,
       aes(x = comp, stratum = change, alluvium = GENEID,
           fill = tf.family, label = change, color=tf.family)) +
  # scale_fill_manual(values = cbPalette) +
  # scale_color_manual(values = cbPalette) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback") +
  geom_stratum(alpha=0.5, aes(fill = change)) +
  theme(legend.position = "bottom") +
  xlab("Tissue Comparision") +
  ylab("Number of TFs")


library(dplyr)

# Count total each family
tf_count <- all_sig_comb %>% count(WGCNA.Module)
tf_count

write.csv(tf_count, "tf.sig.counts.bychange.and.comp.csv")
colnames(sig.tfs)
siglfc_count <- sig.tfs %>% count(WGCNA.Module)
siglfc_count

library(tidyr)
library(dplyr)


tf_gene_count <- all_sig_comb %>% count(GENEID, comp)

tf_wide <- spread(tf_gene_count, comp, n)

tf_wide[is.na(tf_wide)] <- 0
colnames(tf_wide)
tf_wide_num <- tf_wide[,-1]
tf_wide_num <- as.data.frame(tf_wide_num)

setnames <- colnames(tf_wide_num)

library(UpSetR)

upset(tf_wide_num, keep.order = TRUE, sets=setnames, nintersects = NA)


tf_wide$num.comp <- rowSums(tf_wide[2:6])
library(stringr)
names(tf_wide)<-str_replace_all(names(tf_wide), c(" " = "." , "," = "" ))
colnames(tf_wide)

tf_wide_all_meio <- subset(tf_wide, A.LEP.vs.M.LEP == 1 &
                             A.PAC.vs.M.PAC == 0)
tf_wide_all_meio <- subset(tf_wide, A.LEP.vs.M.LEP == 0 &
                             A.PAC.vs.M.PAC == 1)



