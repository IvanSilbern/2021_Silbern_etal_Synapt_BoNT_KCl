
local({
  
  library(data.table)
  library(ggplot2)
  
  temp <- fread("Figures\\SupplFig_17\\ProteinGroups_HeLa_PhSites_Contrasts.txt", check.names = TRUE)
  
  names(temp)
  
  int_cols <- c("BTX01", "BTX02", "BTX03", #  BTX = BoNT
                "Mock01", "Mock02", "Mock03",
                "UT01", "UT02", "UT03")
  
  # Boxplots
  pdf("Figures\\SupplFig_4\\SupplFig_17A.pdf", width = 9)
  boxplot(log2(temp[, ..int_cols]))
  dev.off()
  
  # BoNT vs Mock
  pdf("Figures\\SupplFig_4\\SupplFig_17B.pdf")
  plot(y = -log10(temp$q.val),
       x = temp$log2FC,
       type = "n",
       xlab = "log2-Ratio",
       ylab = "-log10 qValue")
  points(y = -log10(temp$q.val[!temp$Candidate & temp$Contrast == "BoNT_Mock"]),
         x = temp$log2FC[!temp$Candidate & temp$Contrast == "BoNT_Mock"],
         pch = 21,
         bg = "lightgrey")
  points(y = -log10(temp$q.val[temp$Candidate & temp$Contrast == "BoNT_Mock"]),
         x = temp$log2FC[temp$Candidate & temp$Contrast == "BoNT_Mock"],
         bg = "orange",
         pch = 21)
  abline(v = log2(1/1.5), col = "blue", lty = 3, lwd = 2)
  abline(v = log2(1.5), col = "blue", lty = 3, lwd = 2)
  abline(h = -log10(0.01), col = "blue", lty = 3, lwd = 2)
  
  dev.off()
  
  # Mock vs UT
  pdf("Figures\\SupplFig_4\\SupplFig_17C.pdf")
  plot(y = -log10(temp$q.val),
       x = temp$log2FC,
       type = "n",
       xlab = "log2-Ratio",
       ylab = "-log10 qValue")
  points(y = -log10(temp$q.val[!temp$Candidate & temp$Contrast == "Mock_UT"]),
         x = temp$log2FC[!temp$Candidate & temp$Contrast == "Mock_UT"],
         pch = 21,
         bg = "lightgrey")
  points(y = -log10(temp$q.val[temp$Candidate & temp$Contrast == "Mock_UT"]),
         x = temp$log2FC[temp$Candidate & temp$Contrast == "Mock_UT"],
         bg = "orange",
         pch = 21)
  abline(v = log2(1/1.5), col = "blue", lty = 3, lwd = 2)
  abline(v = log2(1.5), col = "blue", lty = 3, lwd = 2)
  abline(h = -log10(0.01), col = "blue", lty = 3, lwd = 2)
  
  dev.off()
  
  # BoNT vs UT
  pdf("Figures\\SupplFig_4\\SupplFig_17D.pdf")
  plot(y = -log10(temp$q.val),
       x = temp$log2FC,
       type = "n",
       xlab = "log2-Ratio",
       ylab = "-log10 qValue")
  points(y = -log10(temp$q.val[!temp$Candidate & temp$Contrast == "BoNT_UT"]),
         x = temp$log2FC[!temp$Candidate & temp$Contrast == "BoNT_UT"],
         pch = 21,
         bg = "lightgrey")
  points(y = -log10(temp$q.val[temp$Candidate & temp$Contrast == "BoNT_UT"]),
         x = temp$log2FC[temp$Candidate & temp$Contrast == "BoNT_UT"],
         bg = "orange",
         pch = 21)
  abline(v = log2(1/1.5), col = "blue", lty = 3, lwd = 2)
  abline(v = log2(1.5), col = "blue", lty = 3, lwd = 2)
  abline(h = -log10(0.01), col = "blue", lty = 3, lwd = 2)
  dev.off()
  
  })