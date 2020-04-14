expr.dat <- read.csv("TranscriptomeResultsUpdate_2.csv")
print(levels(expr.dat$ExpressionFilter))
 
pass.expr <- subset(expr.dat, ExpressionFilter != "LowExpressionFilter")
present.anther <- subset(pass.expr, ExpressionFilter != "EmbryoDataset")
 
write.csv(present.anther, "anther_transcripts_pass_expression.csv")
expr.dat <- read.csv("anther_transcripts_pass_expression.csv")

# # Exclude ncRNAs 
# colnames(expr.dat)
# print(levels(expr.dat$coding))
# coding <- subset(expr.dat, coding =="coding")
# 
# # replace 0 and N/A values
# print(levels(coding$COG))
# levels(coding$COG) <- c(levels(coding$COG), "None") 
# coding$COG[coding$COG=="#N/A"]  <- "None"
# coding$COG[coding$COG=="N/A"]  <- "None"
# coding$COG[coding$COG=="0"] <- "None"

# A.LEP.vs.PRE
lep_vs_pre_up <- subset(expr.dat,
                        A.LEP.vs.A.PRE.Pval < 0.01 &
                        A.LEP.vs.A.PRE.log2FC > 1)
lep_vs_pre_up$comp <- "A.LepZyg-A.Pre"
lep_vs_pre_up$change <- "Up-regulated"

lep_vs_pre_down <- subset(expr.dat,
                          A.LEP.vs.A.PRE.Pval < 0.01 &
                          A.LEP.vs.A.PRE.log2FC < -1)
lep_vs_pre_down$comp <- "A.LepZyg-A.Pre"
lep_vs_pre_down$change <- "Down-regulated"


pac_vs_lep_up <- subset(expr.dat,
                        A.PAC.vs.A.LEP.Pval < 0.01 &
                        A.PAC.vs.A.LEP.log2FC > 1)
pac_vs_lep_up$comp <- "A.PacDip-A.LepZyg"
pac_vs_lep_up$change <- "Up-regulated"

pac_vs_lep_down <- subset(expr.dat,
                          A.PAC.vs.A.LEP.Pval < 0.01 &
                          A.PAC.vs.A.LEP.log2FC < -1)
pac_vs_lep_down$comp <- "A.PacDip-A.LepZyg"
pac_vs_lep_down$change <- "Down-regulated"


# 
tet_vs_pac_up <- subset(expr.dat,
                        A.TET.vs.A.PAC.Pval < 0.01 &
                        A.TET.vs.A.PAC.log2FC > 1)
tet_vs_pac_up$comp <- "A.MetTet-A.PacDip"
tet_vs_pac_up$change <- "Up-regulated"

# 
tet_vs_pac_down <- subset(expr.dat,
                          A.TET.vs.A.PAC.Pval < 0.01 &
                          A.TET.vs.A.PAC.log2FC < -1)
tet_vs_pac_down$comp <- "A.MetTet-A.PacDip"
tet_vs_pac_down$change <- "Down-regulated"


all_sig_comb <- rbind(tet_vs_pac_down,
                          tet_vs_pac_up,
                          pac_vs_lep_down,
                          pac_vs_lep_up,
                          lep_vs_pre_down,
                          lep_vs_pre_up)

all_sig <- subset(expr.dat,
                  A.LEP.vs.A.PRE.Pval < 0.01 & A.LEP.vs.A.PRE.log2FC > 1 |
                  A.LEP.vs.A.PRE.Pval < 0.01 & A.LEP.vs.A.PRE.log2FC < -1 |
                  A.PAC.vs.A.LEP.Pval < 0.01 & A.PAC.vs.A.LEP.log2FC > 1 |
                  A.PAC.vs.A.LEP.Pval < 0.01 & A.PAC.vs.A.LEP.log2FC < -1 |
                  A.TET.vs.A.PAC.Pval < 0.01 & A.TET.vs.A.PAC.log2FC > 1 |
                  A.TET.vs.A.PAC.Pval < 0.01 & A.TET.vs.A.PAC.log2FC < -1)


library(ggalluvial)
all_sig_comb$comp <- as.factor(all_sig_comb$comp)
print(levels(all_sig_comb$comp))
all_sig_comb$comp <- factor(all_sig_comb$comp,
                                levels(all_sig_comb$comp)[c(1,3,2)])

colnames(all_sig_comb)
plot.data <- all_sig_comb[,c(3,52,53)]

lep_vs_pre_sig <- subset(expr.dat,
                         A.LEP.vs.A.PRE.Pval < 0.01 &
                         A.LEP.vs.A.PRE.log2FC > 1 |
                         A.LEP.vs.A.PRE.Pval < 0.01 &
                         A.LEP.vs.A.PRE.log2FC < -1)

BAnTr <- setdiff(all_sig$BAnTr, lep_vs_pre_sig$BAnTr)
comp <- rep("A.LepZyg-A.Pre", length(BAnTr))
change <- rep("Not significant", length(BAnTr))

lep_vs_pre_ns <- data.frame(BAnTr, comp, change)


pac_vs_lep_sig <- subset(expr.dat,
                         A.PAC.vs.A.LEP.Pval < 0.01 &
                         A.PAC.vs.A.LEP.log2FC > 1 |
                         A.PAC.vs.A.LEP.Pval < 0.01 &
                         A.PAC.vs.A.LEP.log2FC < -1)

BAnTr <- setdiff(all_sig$BAnTr, pac_vs_lep_sig$BAnTr)
comp <- rep("A.PacDip-A.LepZyg", length(BAnTr))
change <- rep("Not significant", length(BAnTr))

pac_vs_lep_ns <- data.frame(BAnTr, comp, change)


tet_vs_pac_sig <- subset(expr.dat,
                         A.TET.vs.A.PAC.Pval < 0.01 &
                         A.TET.vs.A.PAC.log2FC > 1 |
                         A.TET.vs.A.PAC.Pval < 0.01 &
                         A.TET.vs.A.PAC.log2FC < -1)

BAnTr <- setdiff(all_sig$BAnTr, tet_vs_pac_sig$BAnTr)
comp <- rep("A.MetTet-A.PacDip", length(BAnTr))
change <- rep("Not significant", length(BAnTr))

tet_vs_pac_ns <- data.frame(BAnTr, comp, change)

plot.data <- rbind(plot.data, lep_vs_pre_ns,
                   pac_vs_lep_ns,
                   tet_vs_pac_ns)

cbPalette <- c("#E69F00",
               "#56B4E9",
               "#CC79A7")

# cbPalette <- c("#E69F00",
#                "#56B4E9",
#                "#009E73",
#                "#F0E442",
#                "#0072B2",
#                "#D55E00",
#                "#CC79A7")

class(plot.data$change)
plot.data$change <- as.factor(plot.data$change)
print(levels(plot.data$change))
plot.data$change <- factor(plot.data$change,
                                levels(plot.data$change)[c(3,2,1)])

alluvial_plot <- ggplot(plot.data,
                        aes(x = comp, stratum = change, alluvium = BAnTr,
                            fill = change, label = change, color=change)) +
                        scale_fill_manual(values = cbPalette) +
                        scale_color_manual(values = cbPalette) +
                        geom_flow(stat = "alluvium",
                                  lode.guidance = "frontback",
                                  alpha=0.1) +
                        geom_stratum(alpha=0.75) +
                        theme(legend.position = "bottom") +
                        xlab("Comparison") +
                        ylab("Differentially Expressed Genes") + 
                        theme_grey(base_size = 8) +
                        theme(legend.position = "none")


alluvial_plot

# 
# ggsave("Differential.Expression.Anther.only.tiff",
#        dpi=600, 
#        width=8,
#        height=6.621,
#        units="cm")








