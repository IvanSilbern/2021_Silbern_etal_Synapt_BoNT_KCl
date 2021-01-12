# Do: plot log2FC of Wnt-related proteins

# INPUT: 
# "temp\\PhPeptIntensities7.tsv"
# "temp\\Protein_classification.txt"

# OUTPUT:
# "plots\\K-Channel_log2FC.pdf"
# "Figures\\SupplFig_16\\log2FC_K_channels.txt"


local({
  
  library(ggplot2)
  library(data.table)
  library(stringr)
  
  if(!dir.exists("Figures\\SupplFig_16")) dir.create("Figures\\SupplFig_16")
  
  # phosphosites data
  df <- fread("temp\\PhPeptIntensities6.tsv")
  df <- df[Regulation_group_resolved != "not-affected"]
  
  # protein classification
  prot_class <- fread("temp\\Protein_classification.tsv")
  
  selected <- prot_class[grepl("K-Channel", Keyword_manual)]$Gene.name
  
  df <- df[Gene.name %in% selected]
  df[, nk_group := unlist(lapply(str_split(df$netphorest_group, ";"), function(x){
    
    if(length(x) > 1){ 
      
      return(x[length(x) - 1])
      
    } else {
      
      return(NA)
      
    }
    
  }))]
  df[, nk_group := gsub("_group", "", nk_group)]
  
  id_cols <- c("Gene.name", "Position", "Multiplicity", "Regulation_group_resolved", "nk_group")
  measure_cols <- c("log2FC.CaEGTA", "log2FC.BoNT")
  
  df <- melt(df[, .SD, .SDcols = c(id_cols, measure_cols)],
             mearsure.vars = measure_cols,
             id.vars = id_cols,
             variable.name = "Experiment", value.name = "log2FC")
  
  df[, Experiment := gsub("log2FC\\.", "", Experiment)]
  df[Experiment == "CaEGTA", Experiment := "Ca/EGTA"]
  df[Experiment == "BoNT", Experiment := "Mock/BoNT"]
  
  df[, Experiment := factor(Experiment, levels = c("Ca/EGTA", "Mock/BoNT"))]
  df[, Site_id := paste0(Gene.name, "-", Position, "(", Multiplicity, ")")]
  
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
  
  pdf("plots\\K-Channel_log2FC.pdf", width = 4, height = 6)
  print(g)
  dev.off()
  
  # save figure data
  fwrite(df, "Figures\\SupplFig_16\\log2FC_K_channels.txt", sep = "\t")
  
})