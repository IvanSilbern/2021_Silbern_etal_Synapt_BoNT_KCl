local({
  
  library(ggplot2)
  library(data.table)
  library(stringr)
  
  if(!dir.exists("Figures\\SupplFig_7")) dir.create("Figures\\SupplFig_7", recursive = TRUE)
  
  df <- fread("temp\\PhPeptIntensities6.tsv")
  
  # subset phosphatase-related sites
  
  df <- df[grepl("Phosphatase", Keyword_manual)]  
  
  # keep regulated sites
  
  df <- df[Regulation_group_resolved != "not-affected"]
  
  
 
  id_cols <- c("Gene.name", "Position", "Multiplicity", "Regulation_group_resolved")
  measure_cols <- c("log2FC.CaEGTA", "log2FC.BoNT")
  
  df <- melt(df[, .SD, .SDcols = c(id_cols, measure_cols)],
             mearsure.vars = measure_cols,
             id.vars = id_cols,
             variable.name = "Experiment", value.name = "log2FC")
  
  df[, Experiment := gsub("log2FC\\.", "", Experiment)]
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
  
  pdf("plots\\phosphatases.pdf", width = 4)
  print(g)
  dev.off()
  
  # export figure data
  fwrite(df, "Figures\\SupplFig_7\\regulated_phosphatases.txt", sep = "\t")

})
  