expr.dat <- read.csv("anther_transcripts_pass_expression.csv")

# Exclude ncRNAs 
colnames(expr.dat)
print(levels(expr.dat$coding))
coding <- subset(expr.dat, coding =="coding")

# replace 0 and N/A values
print(levels(coding$COG))
levels(coding$COG) <- c(levels(coding$COG), "None") 
coding$COG[coding$COG=="#N/A"]  <- "None"
coding$COG[coding$COG=="N/A"]  <- "None"
coding$COG[coding$COG=="0"] <- "None"

# A.LepZyg-A.Pre
lep_vs_pre_up <- subset(coding,
                        A.LEP.vs.PRE.Pval < 0.01 &
                          A.LEP.vs.PRE.log2FC > 1)
lep_vs_pre_up$comp <- "A.LepZyg-A.Pre up-regulated"

lep_vs_pre_down <- subset(coding,
                          A.LEP.vs.PRE.Pval < 0.01 &
                            A.LEP.vs.PRE.log2FC < -1)
lep_vs_pre_down$comp <- "A.LepZyg-A.Pre down-regulated"

# A.PacDip-A.LepZyg
pac_vs_lep_up <- subset(coding,
                        A.PAC.vs.A.LEP.Pval < 0.01 &
                          A.PAC.vs.A.LEP.log2FC > 1)
pac_vs_lep_up$comp <- "A.PacDip-A.LepZyg up-regulated"

pac_vs_lep_down <- subset(coding,
                          A.PAC.vs.A.LEP.Pval < 0.01 &
                            A.PAC.vs.A.LEP.log2FC < -1)
pac_vs_lep_down$comp <- "A.PacDip-A.LepZyg down-regulated"

# A.MetTet-A.PacDip
tet_vs_pac_up <- subset(coding,
                        TET.vs.A.PAC.Pval < 0.01 &
                          TET.vs.A.PAC.log2FC > 1)
tet_vs_pac_up$comp <- "A.MetTet-A.PacDip up-regulated"

tet_vs_pac_down <- subset(coding,
                          TET.vs.A.PAC.Pval < 0.01 &
                            TET.vs.A.PAC.log2FC < -1)
tet_vs_pac_down$comp <- "A.MetTet-A.PacDip down-regulated"

# M.LepZyg-A.LepZyg
m.lep_vs_a.lep_up <- subset(coding,
                            M.LEP.vs.A.LEP.Pval < 0.01 &
                              M.LEP.vs.A.LEP.log2FC > 1)
m.lep_vs_a.lep_up$comp <- "M.LepZyg-A.LepZyg up-regulated"

m.lep_vs_a.lep_down <- subset(coding,
                              M.LEP.vs.A.LEP.Pval < 0.01 &
                                M.LEP.vs.A.LEP.log2FC < -1)
m.lep_vs_a.lep_down$comp <- "M.LepZyg-A.LepZyg down-regulated"

# M.PacDip-A.PacDip
m.pac_vs_a.pac_up <- subset(coding,
                            M.PAC.vs.A.PAC.Pval < 0.01 &
                              M.PAC.vs.A.PAC.log2FC > 1)
m.pac_vs_a.pac_up$comp <- "M.PacDip-A.PacDip up-regulated"

m.pac_vs_a.pac_down <- subset(coding,
                              M.PAC.vs.A.PAC.Pval < 0.01 &
                                M.PAC.vs.A.PAC.log2FC < -1)
m.pac_vs_a.pac_down$comp <- "M.PacDip-A.PacDip down-regulated"

# M.PacDip-M.LepZyg
m.pac_vs_m.lep_up <- subset(coding,
                            M.PAC.vs.M.LEP.Pval < 0.01 &
                              M.PAC.vs.M.LEP.log2FC > 1)
m.pac_vs_m.lep_up$comp <- "M.PacDip-M.LepZyg up-regulated"

m.pac_vs_m.lep_down <- subset(coding,
                              M.PAC.vs.M.LEP.Pval < 0.01 &
                                M.PAC.vs.M.LEP.log2FC < -1)
m.pac_vs_m.lep_down$comp <- "M.PacDip-M.LepZyg down-regulated"


all_sig_comb_anther <- rbind(tet_vs_pac_down,
                             tet_vs_pac_up,
                             pac_vs_lep_down,
                             pac_vs_lep_up,
                             lep_vs_pre_down,
                             lep_vs_pre_up)

all_sig_comb_meiocyte <- rbind(m.lep_vs_a.lep_up,
                               m.lep_vs_a.lep_down,
                               m.pac_vs_a.pac_up,
                               m.pac_vs_a.pac_down,
                               m.pac_vs_m.lep_up,
                               m.pac_vs_m.lep_down)


# Generate an intersect plot of genes with significant changes in expression
library(dplyr)
library(tidyr)
gene_count <- all_sig_comb_anther %>% count(BAnTr, comp)
gene_count_wide <- spread(gene_count, comp, n)

gene_count_wide[is.na(gene_count_wide)] <- 0
colnames(gene_count_wide)
gene_count_wide_num <- gene_count_wide[,-1]
gene_count_wide_num <- as.data.frame(gene_count_wide_num)
setnames <- colnames(gene_count_wide_num)

merged <- merge(gene_count_wide, expr.dat, by = "BAnTr")
colnames(merged)

library(UpSetR)

upset(merged,
      keep.order = TRUE,
      sets=setnames,
      nintersects = NA,
      queries = list(list(query = elements, params = list("coding", c("coding", "ncRNA")), color = "#bf812d", active = T),
                     list(query = elements, params = list("coding", "ncRNA"), color = "#01665e", active = T)),
      text.scale = 2
)

library(magick)
library(cowplot)

alluvial_plot <-cowplot::ggdraw() + cowplot::draw_image("Differential.Expression.Anther.only.svg", scale = 0.8)
meio_upset <- cowplot::ggdraw() + cowplot::draw_image("Meiocyte_upset.R.svg", scale = 0.9)
anth_upset <- cowplot::ggdraw() + cowplot::draw_image("Anther_upset.R.svg", scale = 0.9)

plot_grid(alluvial_plot, anth_upset, meio_upset, ncol=2, nrow=2, labels="auto", align = "v")

draw_image(image_read_svg("Meiocyte_upset.R.svg"), scale = 1)

