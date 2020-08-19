# DO:
# count number of regulated sites per protein and regulateion group
# assign protein to one of the regulation groups
# based on the phosphosites the protein carry
# if > 60% of the sites are Ca-dependent, Ca-compensating, or Cycling-dependent,
# then the protein is also considered as being Ca-dependent, Ca-compensating, or Cycling-dependent, respectively
# if > 60% are Ca-compensating and Cycling-dependent together, then the protein is "Ca-comp." or "Cycling-dep."
# if > 60% are Ca-dependent and Ca-enhancing toghether, then the protein is classified as "Ca-dep. or Ca-enh."
# in other cases the protein is classified as "mixed"

# INPUT:
# "temp\\PhPeptIntensities4.tsv"
#

# OUTPUT:
# "plots//circle_prot_grouped.pdf"
# "Figures\\Fig_4B\\Gene_RegulationGroups.txt"
# "Gene_RegulationGroups.tsv"

local({
  
  if(!dir.exists("Figures\\Fig_4B")) dir.create("Figures\\Fig_4B", recursive = TRUE)
  
  library(data.table)
  library(stringr)
  library(ggplot2)
  
  # minimal percentage of the sites belonging to the same regulation group
  # to consider the protein being mostly regulated the same way
  min_per <- 75
  
  df <- fread("temp\\PhPeptIntensities4.tsv")
  dim(df)
  
  df <- df[df$Significance]
  dim(df)
  
  # remove duplicates due to several multiplicity-states
  df <- df[!duplicated(df$Site_id4), c("Site_id4", "Accession", "Accession.noIso", "Stringid", "REVIEWED",
                                       "Gene.name", "Protein.description", "Regulation_group_resolved", "Function")]
  
  prot_groups_count <- df[, list(Ca_dependent = sum(Regulation_group_resolved == "primary Ca-dependent"),
                                 Cycling_dependent = sum(Regulation_group_resolved == "SV-cycling-dependent")), by = Gene.name]
  prot_groups_per <- prot_groups_count[, lapply(.SD, function(x) 100*x/sum(apply(.SD, 1, sum))),
                                       by = Gene.name, .SDcols = c("Ca_dependent",
                                                                   "Cycling_dependent")]
  prot_groups_count <- merge(prot_groups_count, prot_groups_per, by = "Gene.name", suffixes = c("", "_per"))
  
  prot_groups_count <- merge(prot_groups_count, df[, c("Gene.name", "Accession",
                                                        "Accession.noIso", "Protein.description", "Stringid", "REVIEWED", "Function")],
                             by = "Gene.name", all.x = TRUE)
  prot_groups_count <- prot_groups_count[order(REVIEWED)]
  prot_groups_count <- prot_groups_count[!duplicated(Gene.name)]
  
  # Number of regulated sites in total
  prot_groups_count[, N_reg_sites := Ca_dependent + Cycling_dependent]
  
  prot_groups_count[, Regulation := "Mixed"]
  prot_groups_count[Ca_dependent_per >= min_per, Regulation := "mostly primary Ca-dependent"]
  prot_groups_count[Cycling_dependent_per >= min_per, Regulation := "mostly SV-cycling-dependent"]
 
  n_ca_dep  <- prot_groups_count[Ca_dependent_per > min_per, .N]
  n_cycling <- prot_groups_count[Cycling_dependent_per > min_per, .N]
  other <- prot_groups_count[, .N] - n_ca_dep - n_cycling
  
  pt <- data.table(type = c("mostly primary Ca-dependent",
                            "mostly SV-cycling-dependent",
                            "Mixed"),
                   N = c(n_ca_dep, 
                         n_cycling,
                         other))
  pt[, Percent := 100*N/nrow(prot_groups_count)]
  pt[, type := factor(type, pt$type)]
  
  pt[, type := factor(type, levels = pt$type)]
    
  colors <- c("#fcb533",
              "#51baf4",
              "#c79a8d")
  
  pdf("plots//circle_prot_grouped.pdf", width = 8)
  
  g <- ggplot(pt, aes(x = "", y = Percent, fill = type))
  g <- g + geom_bar(stat = "identity", width = 0.5)
  g <- g + annotate("segment", x = 0.5, xend = 1.25, y = 100-sum(pt$Percent[1:2]), yend = 100-sum(pt$Percent[1:2]), color = "white", size = 2)
  g <- g + annotate("segment", x = 0.5, xend = 1.25, y = 0, yend = 0, color = "white", size = 2)
  g <- g + annotate("segment", x = 0.5, xend = 1.25, y = sum(pt$Percent[2:3]), yend = sum(pt$Percent[2:3]), color = "white", size = 2)
  g <- g + scale_fill_manual(values = scales::alpha(colors, 0.6))
  g <- g + coord_polar("y", start=0)
  g <- g + guides(fill = guide_legend(""))
  g <- g + theme_void()
  g <- g + theme(legend.text = element_text(size = 20),
                 legend.title = element_text(size = 20),
                 legend.key.size = unit(1.5, "cm"))
  
  print(g)
  
  dev.off()
  
  # write tables
  fwrite(prot_groups_count, "SupplData04_Protein_classification.tsv", sep = "\t")
  fwrite(prot_groups_count, "Figures\\Fig_4B\\Gene_RegulationGroups.txt", sep = "\t")
  
})
