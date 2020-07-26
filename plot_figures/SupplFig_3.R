local({
  
  library(data.table)
  library(ggplot2)
  
  df_long <- fread("Figures\\SupplFig_3\\Proteingroup_Intensity_Suprenatants_Cbotulinum.txt")
  df_long <- df_long[order(-iBAQ.log)]
  selected_genes <- c(unique(df_long$Gene.name)[1:30], c("botA", "botB", "botC", "botD", "CBO3541", "CBO3379", "CBO2497"))
  df_long[, Gene.name := factor(Gene.name, levels = unique(df_long$Gene.name))]
  
  pdf("Figures\\SupplFig_3\\SupplFig_3.pdf", width = 12)
  g <- ggplot(df_long[Gene.name %in% selected_genes], aes(y = iBAQ.log, x = Gene.name, color = strain))
  g <- g + geom_point(size = 2, alpha = 0.8)
  g <- g + theme_bw()
  g <- g + theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
                 axis.text.y = element_text(size = 12),
                 axis.title = element_text(size = 14))
  g <- g + xlab("Gene name") + ylab("Protein intensity (log2 iBAQ)")
  g <- g + guides(color = guide_legend(title = "C.botulinum strain"))
  print(g)
  dev.off()
  
  })