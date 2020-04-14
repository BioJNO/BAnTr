expr.dat <- read.csv("expression_dat_with_annotation.csv")

# Exclude ncRNAs 
colnames(expr.dat)
print(levels(expr.dat$coding.ncRNA))
coding <- subset(expr.dat, coding.ncRNA =="coding")

# replace 0 and N/A values
print(levels(coding$COG))
levels(coding$COG) <- c(levels(coding$COG), "None") 
coding$COG[coding$COG=="#N/A"]  <- "None"
coding$COG[coding$COG=="N/A"]  <- "None"
coding$COG[coding$COG=="0"] <- "None"


all_sig <- subset(coding,
                  alep.mlep.adj.P.Val < 0.01 & alep.mlep.logFC > 1 |
                  alep.mlep.adj.P.Val < 0.01 & alep.mlep.logFC > 1 |
                  apac.mpac.adj.P.Val <0.01 & apac.mpac.logFC > 1 |
                  apac.mpac.adj.P.Val <0.01 & apac.mpac.logFC < -1 |
                  lep.pre.adj.P.Val <0.01 & lep.pre.logFC > 1 |
                  lep.pre.adj.P.Val <0.01 & lep.pre.logFC < -1 |
                  pac.lep.adj.P.Val < 0.01 & pac.lep.logFC > 1 |
                  pac.lep.adj.P.Val < 0.01 & pac.lep.logFC < -1 |
                  tet.pac.adj.P.Val < 0.01 & tet.pac.logFC > 1 |
                  tet.pac.adj.P.Val < 0.01 & tet.pac.logFC < -1)
                    
COG_count_total <- all_sig %>% count(COG, WGCNA.Module)
COG_count_total
write.csv(COG_count_total, "total_sig_COG_count_per_module.csv")

colnames(all_sig)

log_fc <- all_sig[,c(2,21:30)]
colnames(log_fc)


for (x in 2:ncol(log_fc)) {
    log_fc[,x] <- as.numeric(log_fc[,x])
}

for (x in c(2,4,6,8,10)) {
    for (y in 1:nrow(log_fc)) {
        if(log_fc[y,x] > 1 & log_fc[y,x+1] < 0.01) {
           log_fc[y,x] <- "Up-regulated"
      } else if(log_fc[y,x] < -1 & log_fc[y,x+1] < 0.01) {
           log_fc[y,x] <- "Down-regulated"
      } else {
           log_fc[y,x] <- "Not significant"
      }
    }
}


long_fc <- gather(log_fc, key = GENEID, value = n)

count <- log_fc %>% count(lep.pre.logFC, alep.mlep.logFC,
                          pac.lep.logFC, apac.mpac.logFC,
                          tet.pac.logFC)

count_stage <- log_fc %>% count(lep.pre.logFC, pac.lep.logFC,
                          tet.pac.logFC)

library(RColorBrewer)
cbPalette <- c("#E69F00",
               "#56B4E9",
               "#009E73",
               "#F0E442",
               "#0072B2",
               "#D55E00",
               "#CC79A7")
library(alluvial)

alluvial(count_stage[,1:3], freq=count_stage$n,
         cex = 0.7
)


library(ggalluvial)
data(majors)
majors$curriculum <- as.factor(majors$curriculum)
ggplot(majors,
       aes(x = semester, stratum = curriculum, alluvium = student,
           fill = curriculum, label = curriculum)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom") +
  ggtitle("student curricula across several semesters")
