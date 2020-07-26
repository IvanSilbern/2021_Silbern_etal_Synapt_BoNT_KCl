local({
  
  library(data.table)
  library(ggplot2)
  
  temp <- fread("Figures\\SupplFig_4ABC\\ProteinGroups_MockBTX_AD.txt")
  temp_sub <- temp[Gene.name %in% c("Vamp2", "Stx1a", "Snap25", "Actb")]
  
  int_cols_norm <- c("Mock_01.norm", "Mock_02.norm", "Mock_03.norm",
                     "BTX_01.norm", "BTX_02.norm", "BTX_03.norm")
  temp_sub <- melt(temp_sub[, c("Gene.name", ..int_cols_norm)], measure.vars = int_cols_norm)
  temp_sub[, variable := gsub("\\.norm$", "", variable)]
  temp_sub[, variable := factor(variable, levels = c("Mock_01", "Mock_02", "Mock_03", "BTX_01", "BTX_02", "BTX_03"))]
  
  
  pdf("Figures\\SupplFig_4ABC\\SupplFig_4A.pdf")
  boxplot(temp[, c("Mock_01", "Mock_02", "Mock_03", "BTX_01", "BTX_02", "BTX_03")])
  dev.off()
  
  pdf("Figures\\SupplFig_4ABC\\SupplFig_4B.pdf")
  plot(y = -log10(temp$p.mod),
       x = temp$log2FC,
       type = "n",
       xlab = "log2-Ratio",
       ylab = "-log10 modarated p-value",
       xlim = c(-3, 2),
       ylim = c(0, 7))
  points(y = -log10(temp$p.mod[!temp$Candidate]),
         x = temp$log2FC[!temp$Candidate],
         pch = 21,
         bg = "lightgrey")
  points(y = -log10(temp$p.mod[temp$Candidate]),
         x = temp$log2FC[temp$Candidate],
         bg = "orange",
         pch = 21)
  
  abline(v = log2(1/1.5), col = "blue", lty = 3, lwd = 2)
  abline(v = log2(1.5), col = "blue", lty = 3, lwd = 2)
  
  dev.off()
  
  pdf("Figures\\SupplFig_4ABC\\SupplFig_4C.pdf", width = 6, height = 4)
  g <- ggplot(temp_sub, aes(x = variable, y = value, group = Gene.name, color = Gene.name))
  g <- g + geom_line()
  g <- g + geom_point()
  g <- g + scale_y_continuous(limits = c(-0.6, 0.6))
  g <- g + geom_hline(yintercept = c(-0.3, 0.3), color = "blue", linetype = 2)
  g <- g + theme_bw()
  g <- g + xlab("Replicate") + ylab("log2 normalized intensity") 
  g <- g + guides(color = guide_legend(title = "Protein"))
  print(g)
  dev.off()
  
})