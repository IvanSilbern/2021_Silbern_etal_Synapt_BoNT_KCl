local({
  
  library(data.table)
  library(ggplot2)
  
  wnt <- fread("Figures\\SupplFig_16\\log2FC_Wnt.txt")
  kch <- fread("Figures\\SupplFig_16\\log2FC_K_channels.txt")
  spectr <- fread("Figures\\SupplFig_16\\log2FC_Spectr_Add.txt")
    
  g <- ggplot(wnt, aes(y = Site_id, x = Experiment, fill = log2FC))
  g <- g + coord_equal()
  g <- g + scale_fill_gradient2(high = "#eeaaff", low = "#aad400", na.value = "lightgrey", limits = c(-1.5, 1.5))
  g <- g + geom_tile()
  g <- g + xlab("") + ylab("")
  g <- g + guides(fill = guide_colorbar(title = expression(paste("log"[2], " FC"))))
  g <- g + theme_minimal()
  g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pdf("Figures\\SupplFig_16\\SupplFig_16_Wnt.pdf", width = 4, height = 10)
  print(g)
  dev.off()
  
  g <- ggplot(kch, aes(y = Site_id, x = Experiment, fill = log2FC))
  g <- g + coord_equal()
  g <- g + scale_fill_gradient2(high = "#eeaaff", low = "#aad400", na.value = "lightgrey", limits = c(-1.5, 1.5))
  g <- g + geom_tile()
  g <- g + xlab("") + ylab("")
  g <- g + guides(fill = guide_colorbar(title = expression(paste("log"[2], " FC"))))
  g <- g + theme_minimal()
  g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pdf("Figures\\SupplFig_16\\SupplFig_16_K_Channel.pdf", width = 4, height = 6)
  print(g)
  dev.off()
  
  g <- ggplot(spectr, aes(y = Site_id, x = Experiment, fill = log2FC))
  g <- g + coord_equal()
  g <- g + scale_fill_gradient2(high = "#eeaaff", low = "#aad400", na.value = "lightgrey", limits = c(-1.5, 1.5))
  g <- g + geom_tile()
  g <- g + xlab("") + ylab("")
  g <- g + guides(fill = guide_colorbar(title = expression(paste("log"[2], " FC"))))
  g <- g + theme_minimal()
  g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pdf("Figures\\SupplFig_16\\SupplFig_16_Spectr.pdf", width = 4, height = 4.8)
  print(g)
  dev.off()
  
  
  })