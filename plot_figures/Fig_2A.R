local({
  
  library(data.table)
  library(ggplot2)
  library(VennDiagram)
  
  df <- fread("Figures\\Fig_2A\\nr_quant_sites.txt")
  
  df[, data := factor(data, levels = data)]
  df[, Experiment := factor(Experiment, levels = c("Silbern et al:\nCa vs EGTA &\nNo-Toxin vs Toxin",
                                                   "Silbern et al:\nCa vs EGTA",
                                                   "Engholm-Keller et al:\nKCl vs Mock-Stim."))]
  df[, Localization_prob := factor(Localization_prob, levels = c("Any", "> 0.75"))]
  
  pdf("Figures\\Fig_2A\\Fig_2A.pdf", width = 10)
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
  print(g)
  dev.off()
  
  
})