
# Do:
# prepare Rdata files for shiny app

# INPUT:
# "temp\\PhPeptInt_long.tsv"
# "temp\\PhPeptInt_BoNT.tsv"
# "temp\\PhPeptInt_CaEGTA.tsv"
# "temp\\Domains.tsv"

# OUTPUT:
# "ShinyApp\\dat_int.Rdata"
# "ShinyApp\\df_bont.Rdata"
# "ShinyApp\\df_caegta.Rdata"
# "ShinyApp\\dom.Rdata"
# "ShinyApp\\gn_acc.Rdata"
# "ShinyApp\\genes_all.Rdata"
# "ShinyApp\\genes_reg.Rdata"
# "ShinyApp\\positions.Rdata"



local({
  
  if(!dir.exists("ShinyApp")) dir.create("ShinyApp")
  library(data.table)
  
  # prepare data 
  
  # normalized intensities in long format
  dat_int <- fread("temp\\PhPeptInt_long.tsv")
  
  save(dat_int, file = "ShinyApp\\dat_int.Rdata")
  
  # log2FC
  df_bont <- fread("temp\\PhPeptInt_BoNT.tsv")
  save(df_bont, file = "ShinyApp\\df_bont.Rdata")
  
  df_caegta <- fread("temp\\PhPeptInt_CaEGTA.tsv")
  save(df_caegta, file = "ShinyApp\\df_caegta.Rdata")
  
  # annotated domains
  dom <- fread("temp\\Domains.tsv")
  save(dom, file = "ShinyApp\\dom.Rdata")
  
  # gene names and accessions
  gn_acc <- rbind(df_bont[, c("Gene.name", "Accession", "Regulation_group", "Protein.description", "Function")],
                  df_caegta[, c("Gene.name", "Accession", "Regulation_group", "Protein.description", "Function")])
  
  # extract gene names
  genes_all <- sort(unique(gn_acc$Gene.name))
  genes_reg <- sort(unique(gn_acc$Gene.name[gn_acc$Regulation_group != "not-regulated"]))
  
  # sort accessions
  gn_acc <- gn_acc[, list(Gene.name = Gene.name[1],
                          Protein.description = Protein.description[1],
                          Function = Function[1],
                          Nr_reg_sites = sum(Regulation_group != "not-regulated")),
                   by = "Accession"]
  gn_acc <- gn_acc[order(-Nr_reg_sites, Gene.name)]
  
  save(gn_acc, file = "ShinyApp\\gn_acc.Rdata")
  save(genes_all, file = "ShinyApp\\genes_all.Rdata")
  save(genes_reg, file = "ShinyApp\\genes_reg.Rdata")
  
  
  # modified positions
  positions <- rbind(df_bont[, c("Accession", "Position")],
                     df_caegta[, c("Accession", "Position")])
  positions <- positions[!duplicated(positions)]
  positions <- positions[order(Position)]
  positions <- positions[, list(list(Position)), by = "Accession"]
  names(positions)[names(positions) == "V1"] <- "Position"
  
  save(positions, file = "ShinyApp\\positions.Rdata")

})
