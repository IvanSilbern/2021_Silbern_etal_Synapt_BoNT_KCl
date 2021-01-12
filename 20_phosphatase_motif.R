# Do: check presence of PP1 and Calcineurin docking motifs within protein sequences

# INPUT:
# "external\\data_Uniprot\\201902_uniprot-rat.fasta"
# "temp\\Stringids_Rat_SynaptosomalProteins.txt"
# "external\\data_Uniprot\\MappingTable_rat_Stringid.txt"
# "temp\\PhPeptIntensities5.tsv"

# OUTPUT:
# "plots\\PP1_substrates_pie_chart.pdf"
# "plots\\CN_substrates_pie_chart.pdf"
# "Figures\\SupplFig_7\\sites_pp1_per.txt"
# "Figures\\SupplFig_7\\sites_pp1_count.txt"
# "Figures\\SupplFig_7\\sites_cn_per.txt"
# "Figures\\SupplFig_7\\sites_cn_count.txt"
# "temp\\PhPeptIntensities6.tsv"

local({
  
  library(ggplot2)
  library(data.table)
  library(stringr)
  
  colors <- c("#fcb533",
              "#51baf4")
  
  if(!dir.exists("Figures\\SupplFig_7")) dir.create("Figures\\SupplFig_7", recursive = TRUE)
  
  # read in fasta file
  fasta.file    <- fread("external\\data_Uniprot\\201902_uniprot-rat.fasta", sep = NULL, header = FALSE)[[1]]
  
  # extract header positions
  ff_header_pos <- which(grepl("^>", fasta.file)) # Lines containing headers 
  ff_headers    <- str_split_fixed(fasta.file[ff_header_pos], pattern = " ", n = 2)[, 1] # extract protein id without description
  ff_headers    <- gsub("^>", "", ff_headers)
  
  ff_acc <- str_match(ff_headers, "\\|([^|]+)\\|")[, 2]
  
  # Load stringids of proteins found in Synaptosomes
  prot_ids <- readLines("temp\\Stringids_Rat_SynaptosomalProteins.txt")
  
  # UP mapping table
  up_map <- fread("external\\data_Uniprot\\MappingTable_rat_Stringid.txt", check.names = TRUE)
  
  # genes names of proteins found in synaptosomes
  prot_genes <- up_map[Stringid %in% prot_ids]$Gene.name.main.up
  prot_genes <- unique(prot_genes)
  prot_genes <- prot_genes[!is.na(prot_genes)]
  prot_genes <- prot_genes[prot_genes != ""]
  
  # read phosphosite data
  df <- fread("temp\\PhPeptIntensities5.tsv")
  
  # add protein sequences
  pseq <- character()
  for(i in 1:nrow(df)){
    
    sstart <- ff_header_pos[which(ff_acc == df$Accession[i])] + 1
    send   <- ff_header_pos[which(ff_acc == df$Accession[i]) + 1] - 1
    if(is.na(send)) send <- length(fasta.file)
    
    pseq[i] <- paste(fasta.file[sstart:send], collapse = "")
    
  }
  
  df[, pseq := ..pseq]
  
  # check presence of PP1 or Calcineurin docking motif
  
  df[, PP1_RVXF_1 := grepl("[KR].{0,1}[VI][^P][FW]", pseq)]
  df[, PP1_RVXF_2 := grepl("[HKR][ACHKMNQRSTV][V][CHKNQRST][FW]", pseq)]
  df[, PP1_RVXF_3 := grepl("KSQKW", pseq)]
  df[, PP1_SILK := grepl("[GS]IL[RK]", pseq)]
  
  sum(df$PP1_RVXF_1)
  sum(df$PP1_RVXF_2)
  sum(df$PP1_RVXF_3)
  sum(df$PP1_SILK)
  
  df[, CN_PXIXIT := grepl("[PG].[ILAV].[ILAV][THDI]", pseq)]
  df[, CN_LXVP := grepl("L.VP", pseq)]
  
  sum(df$CN_PXIXIT)
  sum(df$CN_LXVP)
  
  df[, PP1_substr := PP1_RVXF_1 | PP1_RVXF_2 | PP1_RVXF_3 | PP1_SILK]
  df[, CN_substr  := CN_PXIXIT | CN_LXVP]
  
  sum(df$PP1_substr)
  sum(df$CN_substr)
  
  # subset regulated sites
  not_regulated <- df[Regulation_group_resolved == "not-affected"]
  df <- df[Regulation_group_resolved != "not-affected"]
   
  count_pp1 <- df[PP1_substr == TRUE, .N, by = "Regulation_group_resolved"]
  count_cn  <- df[CN_substr == TRUE, .N,  by = "Regulation_group_resolved"]
  
  per_pp1 <- df[PP1_substr == TRUE, .N/df[PP1_substr == TRUE, .N], by = "Regulation_group_resolved"]
  per_cn  <- df[CN_substr == TRUE, .N/df[CN_substr == TRUE, .N], by = "Regulation_group_resolved"]
  
  g <- ggplot(per_pp1, aes(x = "", y = V1, fill = Regulation_group_resolved))
  g <- g + geom_bar(stat = "identity")
  g <- g + annotate("segment", x = 0.55, xend = 1.45,
                    y = per_pp1[Regulation_group_resolved == "SV-cycling-dependent"]$V1,
                    yend = per_pp1[Regulation_group_resolved == "SV-cycling-dependent"]$V1,
                    color = "white", size = 1.25)
  g <- g + annotate("segment", x = 0.55, xend = 1.45,
                    y = 0, yend = 0,
                    color = "white", size = 1.25)
  g <- g + annotate("text", label = paste0(count_pp1[Regulation_group_resolved == "SV-cycling-dependent"]$N,
                                           "\n(",
                                           round(100*per_pp1[Regulation_group_resolved == "SV-cycling-dependent"]$V1, 1),
                                           "%)"),
                    x = 1, y = 0.25, color = "white", size = 6)
  
  g <- g + annotate("text", label = paste0(count_pp1[Regulation_group_resolved == "primary Ca-dependent"]$N,
                                           "\n(",
                                           round(100*per_pp1[Regulation_group_resolved == "primary Ca-dependent"]$V1, 1),
                                           "%)"),
                    x = 1, y = 0.75, color = "white", size = 6)
  
  g <- g + coord_polar("y")
  g <- g + theme_void()
  g <- g + scale_fill_manual(values = scales::alpha(colors, 0.8),
                             labels = c(expression(paste("primary Ca"^"2+", "-dependent")),
                                        "SV-cycling-dependent"))
  g <- g + guides(fill = guide_legend(title = "Regulation group", label.hjust = 0))
  print(g)
  
  pdf("plots\\PP1_substrates_pie_chart.pdf", width = 6, height = 4)
  print(g)
  dev.off()
  
  g <- ggplot(per_cn, aes(x = "", y = V1, fill = Regulation_group_resolved))
  g <- g + geom_bar(stat = "identity")
  g <- g + annotate("segment", x = 0.55, xend = 1.45,
                    y = per_cn[Regulation_group_resolved == "SV-cycling-dependent"]$V1,
                    yend = per_cn[Regulation_group_resolved == "SV-cycling-dependent"]$V1,
                    color = "white", size = 1.25)
  g <- g + annotate("segment", x = 0.55, xend = 1.45,
                    y = 0, yend = 0,
                    color = "white", size = 1.25)
  g <- g + annotate("text", label = paste0(count_cn[Regulation_group_resolved == "SV-cycling-dependent"]$N,
                                           "\n(",
                                           round(100*per_cn[Regulation_group_resolved == "SV-cycling-dependent"]$V1, 1),
                                           "%)"),
                    x = 1, y = 0.25, color = "white", size = 6)
  
  g <- g + annotate("text", label = paste0(count_cn[Regulation_group_resolved == "primary Ca-dependent"]$N,
                                           "\n(",
                                           round(100*per_cn[Regulation_group_resolved == "primary Ca-dependent"]$V1, 1),
                                           "%)"),
                    x = 1, y = 0.75, color = "white", size = 6)
  
  g <- g + coord_polar("y")
  g <- g + theme_void()
  g <- g + scale_fill_manual(values = scales::alpha(colors, 0.8),
                             labels = c(expression(paste("primary Ca"^"2+", "-dependent")),
                                        "SV-cycling-dependent"))
  g <- g + guides(fill = guide_legend(title = "Regulation group", label.hjust = 0))
  
  print(g)
  
  pdf("plots\\CN_substrates_pie_chart.pdf", width = 6, height = 4)
  print(g)
  dev.off()
  
  # write figure data
  fwrite(per_pp1,   "Figures\\SupplFig_7\\sites_pp1_per.txt",   sep = "\t")
  fwrite(count_pp1, "Figures\\SupplFig_7\\sites_pp1_count.txt", sep = "\t")
  fwrite(per_cn,    "Figures\\SupplFig_7\\sites_cn_per.txt",    sep = "\t")
  fwrite(count_cn,  "Figures\\SupplFig_7\\sites_cn_count.txt",  sep = "\t")
  
  # update sites table
  fwrite(rbind(df, not_regulated, fill = TRUE), "temp\\PhPeptIntensities6.tsv", sep = "\t")
  
})