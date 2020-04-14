sig.tfs <- read.csv("sig.tfs.logfc.csv")

library(dplyr)

# Count total each family
tf_count <- sig.tfs %>% count(tf.family)
tf_count
#write.csv(tf_count, "tf.sig.lfc.counts.csv")

colnames(sig.tfs)
siglfc_count <- sig.tfs %>% count(wgcna.mod)
siglfc_count

sig_lep_vs_pre_up <- subset(sig.tfs, DiffExpr.lepvspre=="TRUE" &
                              logfc.lepvspre > 0)
TF_count_lepvspre_up <- sig_lep_vs_pre_up %>% count(tf.family)
TF_count_lepvspre_up$comp <- "Leptotene vs Premeiosis"
TF_count_lepvspre_up$change <- "up"

sig_lep_vs_pre_down <- subset(sig.tfs, DiffExpr.lepvspre=="TRUE" &
                                logfc.lepvspre < 0)
TF_count_lepvspre_down <- sig_lep_vs_pre_down %>% count(tf.family)
TF_count_lepvspre_down$comp <- "Leptotene vs Premeiosis"
TF_count_lepvspre_down$change <- "down"

sig_meilep_vs_anthlep_up <- subset(sig.tfs, DiffExpr.meilepvsanthlep=="TRUE" &
                                     logfc.meilepvsanthlep > 0)
TF_count_meivsanth_lep_up <- sig_meilep_vs_anthlep_up %>% count(tf.family)
TF_count_meivsanth_lep_up$comp <- "Meiocyte vs anthers leptotene"
TF_count_meivsanth_lep_up$change <- "up"

sig_meilep_vs_anthlep_down <- subset(sig.tfs, DiffExpr.meilepvsanthlep=="TRUE" &
                                       logfc.meilepvsanthlep < 0)
TF_count_meivsanth_lep_down <- sig_meilep_vs_anthlep_down %>% count(tf.family)
TF_count_meivsanth_lep_down$comp <- "Meiocyte vs anthers leptotene"
TF_count_meivsanth_lep_down$change <- "down"

sig_meipac_vs_anthpac_up <- subset(sig.tfs, DiffExpr.meipacvsanthpac=="TRUE" &
                                     logfc.meipacvsanthpac > 0)
TF_count_meivsanth_pac_up <- sig_meipac_vs_anthpac_up %>% count(tf.family)
TF_count_meivsanth_pac_up$comp <- "Meiocyte vs anthers pachytene"
TF_count_meivsanth_pac_up$change <- "up"

sig_meipac_vs_anthpac_down <- subset(sig.tfs, DiffExpr.meipacvsanthpac=="TRUE" &
                                       logfc.meipacvsanthpac < 0)
TF_count_meivsanth_pac_down <- sig_meipac_vs_anthpac_down %>% count(tf.family)
TF_count_meivsanth_pac_down$comp <- "Meiocyte vs anthers pachytene"
TF_count_meivsanth_pac_down$change <- "down"

sig_pac_vs_lep_up <- subset(sig.tfs, DiffExpr.pacvslep=="TRUE" &
                              logfc.pacvslep > 0)
TF_count_pacvslep_up <- sig_pac_vs_lep_up %>% count(tf.family)
TF_count_pacvslep_up$comp <- "Pachytene vs leptotene"
TF_count_pacvslep_up$change <- "up"

sig_pac_vs_lep_down <- subset(sig.tfs, DiffExpr.pacvslep=="TRUE" &
                              logfc.pacvslep < 0)
TF_count_pacvslep_down <- sig_pac_vs_lep_down %>% count(tf.family)
TF_count_pacvslep_down$comp <- "Pachytene vs leptotene"
TF_count_pacvslep_down$change <- "down"

sig_anatet_vs_pac_up <- subset(sig.tfs, DiffExpr.anavspac=="TRUE" &
                              log.fc.anavspac > 0)
TF_count_anavspac_up <- sig_anatet_vs_pac_up %>% count(tf.family)
TF_count_anavspac_up$comp <- "Anaphase/tetrad vs pachytene"
TF_count_anavspac_up$change <- "up"

sig_anatet_vs_pac_down <- subset(sig.tfs, DiffExpr.anavspac=="TRUE" &
                                 log.fc.anavspac < 0)
TF_count_anavspac_down <- sig_anatet_vs_pac_down %>% count(tf.family)
TF_count_anavspac_down$comp <- "Anaphase/tetrad vs pachytene"
TF_count_anavspac_down$change <- "down"

all_TF_counts <- rbind(TF_count_lepvspre_up, TF_count_lepvspre_down,
                       TF_count_anavspac_up, TF_count_anavspac_down,
                       TF_count_meivsanth_lep_up, TF_count_meivsanth_lep_down,
                       TF_count_meivsanth_pac_up, TF_count_meivsanth_pac_down,
                       TF_count_pacvslep_up, TF_count_pacvslep_down)

write.csv(all_TF_counts, "all_TF_counts.csv")

colnames(all_TF_counts)
TF_count_summary <- all_TF_counts %>% count(comp, tf.family)

mybs <- sig.tfs[grepl("MYB", sig.tfs$tf.family),]
write.csv(mybs, "mybs.csv")
