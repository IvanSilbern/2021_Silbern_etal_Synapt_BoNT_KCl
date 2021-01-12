local({
  
  library(data.table)
  library(ggplot2)
  
  df <- fread("Figures\\SupplFig_7\\regulated_phosphatases.txt")
  
  df[, Experiment := factor(Experiment, levels = c("CaEGTA", "BoNT"))]
  df[, Site_id := paste0(Gene.name, "-", Position, "(", Multiplicity, ")")]
  df <- df[, Gene.name := factor(Gene.name, levels = c("Synj1",
                                                       "Ssh3",
                                                       "Ptprd",
                                                       "Ptpn4",
                                                       "Phactr3",
                                                       "Phactr2",
                                                       "Phactr1",
                                                       "Pcp2",
                                                       "Ppp6r2",
                                                       "Ppp3ca",
                                                       "Ppp2r5c",
                                                       "Ppp1r21",
                                                       "Ppp1r16b",
                                                       "Ppp1r13b",
                                                       "Ppp1r12c",
                                                       "Ppp1r12a",
                                                       "Ppp1r9b",
                                                       "Ppp1r9a",
                                                       "Ppp1r7",
                                                       "Ppp1r3e",
                                                       "Ppp1r3d",
                                                       "Ppp1r2",
                                                       "Ppp1r1b",
                                                       "Ppp1r1a"))]
  df <- df[order(Gene.name, -Position, - Multiplicity)]
  df[, Site_id := factor(Site_id, levels = unique(Site_id))]
  
  
  g <- ggplot(df, aes(y = Site_id, x = Experiment, fill = log2FC))
  g <- g + coord_equal()
  g <- g + scale_fill_gradient2(high = "#eeaaff", low = "#aad400", na.value = "lightgrey", limits = c(-1.5, 1.5))
  g <- g + geom_tile()
  g <- g + xlab("") + ylab("")
  g <- g + guides(fill = guide_colorbar(title = expression(paste("log"[2], " FC"))))
  g <- g + theme_minimal()
  g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  pdf("Figures\\SupplFig_7\\SupplFig_7A_regulated_phosphatases.pdf", width = 4)
  print(g)
  dev.off()
  
  })