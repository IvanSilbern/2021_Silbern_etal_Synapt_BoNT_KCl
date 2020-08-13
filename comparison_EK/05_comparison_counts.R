
# Do:
# Calculate number of quantified sites
# Check overlap in sequence windows

# INPUT:
# "comparison_EK\\temp\\Phosphosites_prepared_IS.tsv"
# "comparison_EK\\temp\\Phosphosites_prepared_EK.tsv"

# OUTPUT:
# "comparison_EK\\plots\\BarPlot_Sites_Count.pdf"
# "comparison_EK\\plots\\SeqWind_venn.tif"
# "Figures\\Fig_2A\\nr_quant_sites.txt"
# "comparison_EK\\temp\\nr_quant_sites.txt"
# "comparison_EK\\temp\\nr_quant_prot.txt"
# "comparison_EK\\temp\\nr_quant_seqwind.txt"

local({
    
  if(!dir.exists("Figures\\Fig_2A")) dir.create("Figures\\Fig_2A", recursive = TRUE)
  if(!dir.exists("comparison_EK\\plots")) dir.create("comparison_EK\\plots")
  
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(VennDiagram)


  ph_s <- fread("comparison_EK\\temp\\Phosphosites_prepared_IS.tsv")
  ph_r <- fread("comparison_EK\\temp\\Phosphosites_prepared_EK.tsv")
  dim(ph_s)
  dim(ph_s[Localization.prob > 0.75])
  dim(ph_s[Localization.prob > 0.75 & (Intensity.Ca_EGTA_01 != 0 | Intensity.Ca_EGTA_02 != 0)])


  # check that site has at least 3 valid value (reporter intensities) per experiment
  cnames_caegta_01_1 <- grep("Reporter\\.intensity\\.[1-6]\\.Ca_EGTA_01___1", names(ph_s), value = TRUE)
  cnames_caegta_01_2 <- grep("Reporter\\.intensity\\.[1-6]\\.Ca_EGTA_01___2", names(ph_s), value = TRUE)
  cnames_caegta_01_3 <- grep("Reporter\\.intensity\\.[1-6]\\.Ca_EGTA_01___3", names(ph_s), value = TRUE)

  take_caegta_01_1 <- apply(ph_s[, ..cnames_caegta_01_1], 1, function(x) sum(x == 0) < 4)
  take_caegta_01_2 <- apply(ph_s[, ..cnames_caegta_01_2], 1, function(x) sum(x == 0) < 4)
  take_caegta_01_3 <- apply(ph_s[, ..cnames_caegta_01_3], 1, function(x) sum(x == 0) < 4)

  cnames_caegta_02_1 <- grep("Reporter\\.intensity\\.[1-6]\\.Ca_EGTA_02___1", names(ph_s), value = TRUE)
  cnames_caegta_02_2 <- grep("Reporter\\.intensity\\.[1-6]\\.Ca_EGTA_02___2", names(ph_s), value = TRUE)
  cnames_caegta_02_3 <- grep("Reporter\\.intensity\\.[1-6]\\.Ca_EGTA_02___3", names(ph_s), value = TRUE)

  take_caegta_02_1 <- apply(ph_s[, ..cnames_caegta_02_1], 1, function(x) sum(x == 0) < 4)
  take_caegta_02_2 <- apply(ph_s[, ..cnames_caegta_02_2], 1, function(x) sum(x == 0) < 4)
  take_caegta_02_3 <- apply(ph_s[, ..cnames_caegta_02_3], 1, function(x) sum(x == 0) < 4)

  take_caegta_01 <- take_caegta_01_1 | take_caegta_01_2 | take_caegta_01_3
  take_caegta_02 <- take_caegta_02_1 | take_caegta_02_2 | take_caegta_02_3

  cnames_BoNT_01_1 <- grep("Reporter\\.intensity\\.[1-6]\\.BoNT_AD_01___1", names(ph_s), value = TRUE)
  cnames_BoNT_01_2 <- grep("Reporter\\.intensity\\.[1-6]\\.BoNT_AD_01___2", names(ph_s), value = TRUE)
  cnames_BoNT_01_3 <- grep("Reporter\\.intensity\\.[1-6]\\.BoNT_AD_01___3", names(ph_s), value = TRUE)

  take_BoNT_01_1 <- apply(ph_s[, ..cnames_BoNT_01_1], 1, function(x) sum(x == 0) < 4)
  take_BoNT_01_2 <- apply(ph_s[, ..cnames_BoNT_01_2], 1, function(x) sum(x == 0) < 4)
  take_BoNT_01_3 <- apply(ph_s[, ..cnames_BoNT_01_3], 1, function(x) sum(x == 0) < 4)

  cnames_BoNT_02_1 <- grep("Reporter\\.intensity\\.[1-6]\\.BoNT_CB_02___1", names(ph_s), value = TRUE)
  cnames_BoNT_02_2 <- grep("Reporter\\.intensity\\.[1-6]\\.BoNT_CB_02___2", names(ph_s), value = TRUE)
  cnames_BoNT_02_3 <- grep("Reporter\\.intensity\\.[1-6]\\.BoNT_CB_02___3", names(ph_s), value = TRUE)

  take_BoNT_02_1 <- apply(ph_s[, ..cnames_BoNT_02_1], 1, function(x) sum(x == 0) < 4)
  take_BoNT_02_2 <- apply(ph_s[, ..cnames_BoNT_02_2], 1, function(x) sum(x == 0) < 4)
  take_BoNT_02_3 <- apply(ph_s[, ..cnames_BoNT_02_3], 1, function(x) sum(x == 0) < 4)

  take_BoNT_01 <- take_BoNT_01_1 | take_BoNT_01_2 | take_BoNT_01_3
  take_BoNT_02 <- take_BoNT_02_1 | take_BoNT_02_2 | take_BoNT_02_3

  sum(take_caegta_01)
  sum(take_caegta_02)
  sum(take_BoNT_01)
  sum(take_BoNT_02)

  take_caegta <- take_caegta_01 | take_caegta_02
  take_all    <- take_BoNT_01 | take_BoNT_02 | take_caegta_01 | take_caegta_02

  sum(take_caegta)
  sum(take_all)

  df <- data.table(data = c("Engholm_all", "Engholm_1st_class",
                            "Silbern_total", "Silbern_total_1st_class",
                            "Silbern_CaEGTA", "Silbern_CaEGTA_1st_class"),
                   counts = c(nrow(ph_r), nrow(ph_r[Localization.prob > 0.75]),
                              nrow(ph_s[take_all]), nrow(ph_s[take_all & Localization.prob > 0.75]),
                              nrow(ph_s[take_caegta]),
                              nrow(ph_s[take_caegta & 
                                          Localization.prob > 0.75])),
                   Study = c("Engholm-Keller et al", "Engholm-Keller et al",
                             "Silbern et al", "Silbern et al",
                             "Silbern et al", "Silbern et al"),
                   Experiment = c("Engholm-Keller et al:\nKCl vs Mock-Stim.", "Engholm-Keller et al:\nKCl vs Mock-Stim.",
                                  "Silbern et al:\nCa vs EGTA &\nNo-Toxin vs Toxin", "Silbern et al:\nCa vs EGTA &\nNo-Toxin vs Toxin",
                                  "Silbern et al:\nCa vs EGTA", "Silbern et al:\nCa vs EGTA"),
                   Localization_prob = c("Any", "> 0.75",
                                         "Any", "> 0.75",
                                         "Any", "> 0.75")) 
  
  df_prot <- data.table(data = c("Engholm_all", "Engholm_1st_class",
                                       "Silbern_total", "Silbern_total_1st_class",
                                       "Silbern_CaEGTA", "Silbern_CaEGTA_1st_class"),
                              counts = c(length(unique(ph_r$Gene.name)),
                                         length(unique(ph_r$Gene.name[ph_r$Localization.prob > 0.75])),
                                         length(unique(ph_s$Gene.name[take_all])),
                                         length(unique(ph_s$Gene.name[take_all & ph_s$Localization.prob > 0.75])),
                                         length(unique(ph_s$Gene.name[take_caegta])),
                                         length(unique(ph_s$Gene.name[take_caegta & ph_s$Localization.prob > 0.75]))),
                              Study = c("Engholm-Keller et al", "Engholm-Keller et al",
                                        "Silbern et al", "Silbern et al",
                                        "Silbern et al", "Silbern et al"),
                              Experiment = c("Engholm-Keller et al:\nKCl vs Mock-Stim.", "Engholm-Keller et al:\nKCl vs Mock-Stim.",
                                             "Silbern et al:\nCa vs EGTA &\nNo-Toxin vs Toxin", "Silbern et al:\nCa vs EGTA &\nNo-Toxin vs Toxin",
                                             "Silbern et al:\nCa vs EGTA", "Silbern et al:\nCa vs EGTA"),
                              Localization_prob = c("Any", "> 0.75",
                                                    "Any", "> 0.75",
                                                    "Any", "> 0.75")) 
  
  df_sw <- data.table(data = c("Engholm_all", "Engholm_1st_class",
                                 "Silbern_total", "Silbern_total_1st_class",
                                 "Silbern_CaEGTA", "Silbern_CaEGTA_1st_class"),
                        counts = c(length(unique(ph_r$Sequence.window_15)),
                                   length(unique(ph_r$Sequence.window_15[ph_r$Localization.prob > 0.75])),
                                   length(unique(ph_s$Sequence.window_15[take_all])),
                                   length(unique(ph_s$Sequence.window_15[take_all & ph_s$Localization.prob > 0.75])),
                                   length(unique(ph_s$Sequence.window_15[take_caegta])),
                                   length(unique(ph_s$Sequence.window_15[take_caegta & ph_s$Localization.prob > 0.75]))),
                        Study = c("Engholm-Keller et al", "Engholm-Keller et al",
                                  "Silbern et al", "Silbern et al",
                                  "Silbern et al", "Silbern et al"),
                        Experiment = c("Engholm-Keller et al:\nKCl vs Mock-Stim.", "Engholm-Keller et al:\nKCl vs Mock-Stim.",
                                       "Silbern et al:\nCa vs EGTA &\nNo-Toxin vs Toxin", "Silbern et al:\nCa vs EGTA &\nNo-Toxin vs Toxin",
                                       "Silbern et al:\nCa vs EGTA", "Silbern et al:\nCa vs EGTA"),
                        Localization_prob = c("Any", "> 0.75",
                                              "Any", "> 0.75",
                                              "Any", "> 0.75"))
  
  sw <- list(s = unique(ph_s$Sequence.window_7[ph_s$Localization.prob > 0.75 & take_caegta]),
             e = unique(ph_r$Sequence.window_7[ph_r$Localization.prob > 0.75]))

  venn.diagram(sw, filename = "comparison_EK\\plots\\SeqWind_venn.tif",
               category.names = c("Silbern et al", "Engholm-Keller et al"),
               cex = 2, cat.pos = c(0, 21), cat.cex = 1.6)
  
  df[, data := factor(data, levels = data)]
  df[, Experiment := factor(Experiment, levels = c("Silbern et al:\nCa vs EGTA &\nNo-Toxin vs Toxin",
                                                   "Silbern et al:\nCa vs EGTA",
                                                   "Engholm-Keller et al:\nKCl vs Mock-Stim."))]
  df[, Localization_prob := factor(Localization_prob, levels = c("Any", "> 0.75"))]
  
  fwrite(df, "comparison_EK\\temp\\nr_quant_sites.txt")
  fwrite(df, "Figures\\Fig_2A\\nr_quant_sites.txt")
  fwrite(df_prot, "comparison_EK\\temp\\nr_quant_prot.txt")
  fwrite(df_sw, "comparison_EK\\temp\\nr_quant_seqwind.txt")

  pdf("comparison_EK\\plots\\BarPlot_Sites_Count.pdf", width = 10)
  g <- ggplot(df, aes(y = counts, x = Localization_prob, fill = Localization_prob))
  g <- g + facet_grid(~Experiment)
  g <- g + geom_col(alpha = 0.8)
  g <- g + scale_fill_manual(values = c("lightblue", "orange"))
  g <- g + scale_y_continuous(expand = c(0.01, 0.1))
  g <- g + ylab("Number of quantified Phosphorylation sites\n") + xlab ("")
  g <- g + theme_bw()
  g <- g + theme(axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.y = element_text(size = 16),
                 axis.title.y = element_text(size = 24),
                 strip.text = element_text(size = 15),
                 legend.text = element_text(size = 18),
                 legend.title = element_text(size = 18))
  g <- g + guides(fill = guide_legend(title = "Localization Prob.\n(MaxQuant)"))
  g
  dev.off()

})
