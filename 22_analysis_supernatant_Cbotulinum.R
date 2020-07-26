
# Do:
# Analyse protein intensities
# C.botulinum cell culture supernatants

# INPUT:
# "search_results//Supernatants_Cbotulinum//combined//txt//proteinGroups.txt

# OUTPUT:
# "plots\\supernatants_Cbotulinum_proteins_ranked_first30.pdf"
# "Figures\\SupplFig_3\\Proteingroup_Intensity_Suprenatants_Cbotulinum.txt"


local({

  if(!dir.exists("Figures\\SupplFig_3\\")) dir.create("Figures\\SupplFig_3\\", recursive = TRUE)
  if(!dir.exists("plots")) dir.create("plots")
  
  library(ggplot2)
  library(stringr)
  library(data.table)
  
  df <- fread("search_results//Supernatants_Cbotulinum//combined//txt//proteinGroups.txt", check.names = TRUE)
  dim(df)
  head(df)
  
  df <- df[Peptide.counts..all. > 1]
  dim(df)
  
  df <- df[Potential.contaminant != "+"]
  df <- df[Reverse != "+"]
  dim(df)
  
  df[, Gene.name := str_match(Fasta.headers, "GN=([^\\s]+)")[, 2]]
  df[, Protein.description := str_match(Fasta.headers, "\\s(.+) OS=")[, 2]]
  df[Protein.IDs == "sp|P18640|BXC_CBCP", Gene.name := "botC"]
  
  int_cols <- paste0("iBAQ.", c("A", "B", "C", "D"))
                   
  df_long <- melt(df[, c("id", "Gene.name", ..int_cols)], measure.vars = int_cols)
  df_long[, iBAQ.log := log2(value)]
  
  df_long <- df_long[order(-iBAQ.log)]
  selected_genes <- c(unique(df_long$Gene.name)[1:30], c("botA", "botB", "botC", "botD", "CBO3541", "CBO3379", "CBO2497"))
  
  df_long[, Gene.name := factor(Gene.name, levels = unique(df_long$Gene.name))]
  df_long[, variable := gsub("iBAQ\\.", "", variable)]
  
  names(df_long)[names(df_long) == "value"] <- "iBAQ"
  names(df_long)[names(df_long) == "variable"] <- "strain"
  
  pdf("plots\\supernatants_Cbotulinum_proteins_ranked_first30.pdf", width = 12)
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
  
  # provide source data
  fwrite(df_long, "Figures\\SupplFig_3\\Proteingroup_Intensity_Suprenatants_Cbotulinum.txt", sep = "\t")
  
})
