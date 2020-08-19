
# DO:
# test enrichment of kinase-substrate relationships
# in each regulation group
# use pairwise test to assess significance

# INPUT:
# "temp\\PhPeptIntensities4.tsv"
#
# OUTPUT:
# "Figures\\Fig_4C\\np_reg_sites_count.txt"
# "Figures\\Fig_4C\\np_reg_sites_pval.txt"
# "plots\\Fig_NPh_difference_Reg//kinase_reg_group_count.pdf"
# "plots\\Fig_NPh_difference_Reg//kinase_reg_group_pval.pdf"



local({
  
  if(!dir.exists("Figures\\Fig_4C")) dir.create("Figures\\Fig_4C", recursive = TRUE)
  if(!dir.exists("plots")) dir.create("plots")
  
  library(data.table)
  library(stringr)
  library(ggplot2)
  
  df <- fread("temp\\PhPeptIntensities4.tsv")
  dim(df)
  
  df <- df[df$Significance]
  dim(df)
  df <- df[!duplicated(df$Site_id4), c("Site_id4", "Position", "Amino.acid", "Multiplicity", "Accession",
                                       "Gene.name", "Protein.description", "Regulation_group_resolved",
                                       "Sequence.window_7", "Kin_gen", "Kin_mapping", "netphorest_group")]
  df[, netphorest_group := str_split(netphorest_group, ";")]
  
  # extract the second outmost kinase group
  df[, netphorest_group2 := lapply(netphorest_group, function(x) x[max(c(length(x) - 1, 2))])]
  df[, netphorest_group2 := unlist(netphorest_group2)]
  df[, netphorest_group2 := gsub("_group", "", netphorest_group2)]
  
  # extract the highest level kinase group
  df[, netphorest_group_high := unlist(lapply(netphorest_group, function(x) x[length(x)]))]
  
  
  # Approx. the same % of sites (~43-48%) from each regulation group has not kinase assigned
  # Remove sites with no kinase assigned from further analysis
  
  df <- df[!is.na(netphorest_group2)]
  
  # sites per netphorest group and regulation group
  sites_by_np_reg <- df[,.N, by = c("Regulation_group_resolved", "netphorest_group_high", "netphorest_group2")]
  sites_by_np_reg <- sites_by_np_reg[order(netphorest_group2, Regulation_group_resolved)]
  
  # sites per regulation group
  sites_by_reg    <- df[,list(Total = .N), by = c("Regulation_group_resolved")]  
  
  # % sites per regulation group
  sites_by_np_reg <- merge(sites_by_np_reg, sites_by_reg, by = c("Regulation_group_resolved"))
  sites_by_np_reg[, N_per := 100*N / Total]
  
  # max number of sites per netphorest_group2
  max_N <- sites_by_np_reg[, list(max_N = unlist(lapply(.SD, max, na.rm = TRUE))), by = "netphorest_group2", .SDcols = "N"]
  sites_by_np_reg <- merge(sites_by_np_reg, max_N, by = "netphorest_group2")
  
  # remove sites with max_N < 10
  sites_by_np_reg <- sites_by_np_reg[max_N >= 10]
  
  # order Regulation groups
  sites_by_np_reg[, Regulation_group_resolved := factor(Regulation_group_resolved, levels = c("primary Ca-dependent", "SV-cycling-dependent"))]
  
  # order netphorest groups at the highest level
  np_groups_high <- sites_by_np_reg[, sum(max_N), by = netphorest_group_high]
  np_groups_high <- np_groups_high[order(-V1)]
  
  sites_by_np_reg[, netphorest_group_high := factor(netphorest_group_high, levels = np_groups_high$netphorest_group_high)]
  
  # order Netphorest groups
  np_groups_ordered <- unique(sites_by_np_reg$netphorest_group2[order(sites_by_np_reg$netphorest_group_high, -sites_by_np_reg$max_N)])
  sites_by_np_reg[, netphorest_group2 := factor(netphorest_group2, levels = np_groups_ordered)]
  
  # number of sites in the Regulation group that are not regulated by a given kinase
  sites_by_np_reg[, Total_ := Total - N]
  
  # prepare matrix for pairwise tests
  sites_split <- split(sites_by_np_reg, f = sites_by_np_reg$netphorest_group2)
  
  list_matr <- vector("list", length(sites_split))
  for(i in seq_along(sites_split)){
    
    m <- as.matrix(sites_split[[i]][, c("N", "Total_")])
    rownames(m) <- sites_split[[i]]$Regulation_group_resolved
    list_matr[[i]] <- m
    names(list_matr)[i] <- as.character(sites_split[[i]]$netphorest_group2)[1]
    
  }
  
  # testing
  test <- vector("list", length(list_matr))
  names(test) <- names(list_matr)
  for(i in seq_along(list_matr)){
    
    test[[i]] <- tryCatch({fisher.test(list_matr[[i]])}, error = function(e) NA)
    
  }
  
  # extract p-values
  df_tests <- data.table()
  for(i in seq_along(test)){
    
    df_tests <- rbind(df_tests,
                      data.table(netphorest_group2 = names(test)[i],
                                 comparison = c("primary Ca-dependent vs\nSV-cycling-dependent"),
                                 p.val = as.numeric(test[[i]]$p.value))
    )
    
  }
  df_tests <- df_tests[!is.na(p.val)]
  
  # BH correction
  df_tests[, p.adj := p.adjust(p.val, "BH")]
  
  # order np groups
  df_tests[, netphorest_group2 := factor(netphorest_group2, levels = np_groups_ordered)]
  
  df_tests[df_tests$p.adj < 0.10]
  
  # provide figure source
  fwrite(sites_by_np_reg, "Figures\\Fig_4C\\np_reg_sites_count.txt", sep = "\t")
  fwrite(df_tests,        "Figures\\Fig_4C\\np_reg_sites_pval.txt", sep = "\t")
  
  ##### plots #####  
  
  pdf("plots\\kinase_reg_group_count.pdf", width = 8, height = 10)
  
  g <- ggplot(sites_by_np_reg, aes(Regulation_group_resolved, netphorest_group2))
  
  g <- g + geom_tile(aes(fill = N_per), color = "white", size = 1)
  g <- g + coord_equal()
  g <- g + geom_text(aes(label = paste0(round(N_per, 1), "%\n(", N, ")")), color = "grey40")
  g <- g + annotate("text", x = c(1, 2), y = c(0.5, 0.5), label = c("", ""))
  g <- g + annotate("text", x = c(1, 2), y = c(0.95, 0.95),
                    label = c(sites_by_np_reg$Total[sites_by_np_reg$Regulation_group_resolved == "primary Ca-dependent"][1],
                              sites_by_np_reg$Total[sites_by_np_reg$Regulation_group_resolved == "SV-cycling-dependent"][1]),
                    fontface = "bold")
  g <- g + annotate("segment", x = 0.5, xend = 2.5, y = 1.4, yend = 1.4, color = "black", size = 1.5)
  g <- g + scale_fill_gradient(high = "#fcb533", low = "#ebf7e9", na.value = "white", limits = c(0, 30), oob = scales::squish)
  g <- g + scale_x_discrete(expand = c(0, 0), position = "top", breaks = levels(sites_by_np_reg$Regulation_group_resolved))
  g <- g + scale_y_discrete(expand = c(0, 0), limits = rev(c(levels(sites_by_np_reg$netphorest_group2), "Total")))
  g <- g + theme(plot.margin = unit(c(3,3,3,3), "cm"),
                 axis.text.x = element_text(angle = 30, hjust = 0, face = "bold", size = 11),
                 axis.text.y = element_text(size = 16),
                 legend.key.size = unit(1, "cm"),
                 legend.text = element_text(size = 14),
                 legend.title = element_text(size = 18),
                 panel.grid = element_blank(),
                 panel.background = element_rect(fill = "white"),
                 axis.ticks.y = element_blank()
  )
  g <- g + xlab("") + ylab("")
  g <- g + guides(fill = guide_colorbar(title = "% Sites"))
  print(g)
  
  dev.off()
  
  pdf("plots\\kinase_reg_group_pval.pdf", width = 8, height = 10)
  
  g <- ggplot(df_tests, aes(comparison, netphorest_group2))
  g <- g + geom_tile(aes(fill = p.adj), color = "white", size = 1)
  g <- g + coord_equal()
  g <- g + geom_text(aes(label = round(p.adj, 3)), color = "grey40")
  g <- g + scale_fill_gradient(low = "#51baf4", high = "#ebf7e9", na.value = "white", limits = c(0, 0.2), oob = scales::squish)
  g <- g + scale_x_discrete(expand = c(0, 0), position = "top", breaks = levels(df_tests$comparison))
  g <- g + scale_y_discrete(expand = c(0, 0), limits = rev(levels(df_tests$netphorest_group2)))
  g <- g + theme(plot.margin = unit(c(3,3,3,3), "cm"),
                 axis.text.x = element_text(angle = 30, hjust = 0, face = "bold", size = 10),
                 axis.text.y = element_text(size = 16),
                 legend.key.size = unit(1, "cm"),
                 legend.text = element_text(size = 14),
                 legend.title = element_text(size = 18),
                 panel.grid = element_blank(),
                 panel.background = element_rect(fill = "white"),
                 axis.ticks.y = element_blank()
  )
  g <- g + xlab("") + ylab("")
  g <- g + guides(fill = guide_colorbar(title = "adj. p-Value"))
  print(g)
  
  dev.off()
})