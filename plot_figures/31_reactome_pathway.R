local({
  
  library(data.table)
  library(rjson)
  library(stringr)
  library(ggplot2)
  
  sites <- fread("temp\\PhPeptIntensities4.tsv")
  prot  <- fread("temp\\Protein_classification.tsv")
  
  rdb <- fread("external\\data_Reactome\\reactome_result_all.csv", check.names = TRUE)
  rdb_map <- fread("external\\data_Reactome\\reactome_mapping_all.csv", check.names = TRUE)
  
  ns_map <- fromJSON(file = "external\\data_Reactome\\containedEvents_R-HSA-112316_NeuronalSystem.json")
  stids <- unlist(lapply(ns_map, "[[", "stId"))
  
  rdb <- rdb[rdb$Pathway.identifier %in% stids & rdb$Entities.FDR < 0.1]
  
  gene_pathway <- rdb[, list(Gene.name = unlist(str_split(Submitted.entities.found, ";"))), by = Pathway.name]
  
  ca_dep <- sites[Regulation_group == "primary Ca-dependent"]
  ca_dep <- ca_dep[!duplicated(ca_dep$Site_id4)]
  ca_dep <- ca_dep[, .N, by = "Gene.name"]
  names(ca_dep) <- c("Gene.name", "N_ca")
  total_ca <- sum(ca_dep$N_ca)
  
  cycl_dep <- sites[Regulation_group == "SV-cycling-dependent"]
  cycl_dep <- cycl_dep[!duplicated(cycl_dep$Site_id4)]
  cycl_dep <- cycl_dep[, .N, by = "Gene.name"]
  names(cycl_dep) <- c("Gene.name", "N_cycl")
  total_cycl <- sum(cycl_dep$N_cycl)
  
  gene_pathway <- merge(gene_pathway, ca_dep, by = "Gene.name", all.x = TRUE)
  gene_pathway <- merge(gene_pathway, cycl_dep, by = "Gene.name", all.x = TRUE)
  
  path_count <- gene_pathway[, lapply(.SD, sum, na.rm = TRUE), by = "Pathway.name", .SDcols = c("N_ca", "N_cycl")]
  path_count[, Total_ca := ..total_ca - N_ca]
  path_count[, Total_cycl := ..total_cycl - N_ca]
  
  list_mtx <- vector("list", length = path_count[, .N])
  for(i in seq_along(path_count$Pathway.name)){
    
    list_mtx[[i]] <- matrix(as.integer(path_count[i, -c("Pathway.name")]), ncol = 2)
    colnames(list_mtx[[i]]) <- c("N", "Total_")
    rownames(list_mtx[[i]]) <- c("primary Ca-dependent", "SV-cycling-dependent")
    
  }
  names(list_mtx) <- path_count$Pathway.name
  
  # testing
  test <- vector("list", length(list_mtx))
  names(test) <- names(list_mtx)
  for(i in seq_along(list_mtx)){
    
    test[[i]] <- tryCatch({fisher.test(list_mtx[[i]])}, error = function(e) NA)
    
  }
  
  # extract p-values
  df_tests <- data.table()
  for(i in seq_along(test)){
    
    df_tests <- rbind(df_tests,
                      data.table(Pathway.name = names(test)[i],
                                 comparison = c("primary Ca-dependent vs\nSV-cycling-dependent"),
                                 p.val = as.numeric(test[[i]]$p.value))
    )
    
  }
  df_tests <- df_tests[!is.na(p.val)]
  
  # BH correction
  df_tests[, p.adj := p.adjust(p.val, "BH")]
  
  # calculate enrichments
  path_count[, Ca_per := 100*N_ca / (N_ca + total_ca)]
  path_count[, Cycl_per := 100*N_cycl / (N_ca + total_ca)]
  
  path_count <- path_count[order(-Ca_per, -Cycl_per)]
  
  # long format
  path_long_count <- melt(path_count[, c("Pathway.name", "N_ca", "N_cycl")], measure.vars = c("N_ca", "N_cycl"), value.name = "count", variable.name = "Regulation.group")
  path_long_per   <- melt(path_count[, c("Pathway.name", "Ca_per", "Cycl_per")], measure.vars = c("Ca_per", "Cycl_per"), value.name = "percentage", variable.name = "Regulation.group")
  path_long <- cbind(path_long_count, path_long_per[, "percentage"])
  path_long[Regulation.group == "N_ca", Regulation.group := "primary Ca-dependent"]
  path_long[Regulation.group == "N_cycl", Regulation.group := "SV-cycling-dependent"]
  
  # order np groups
  df_tests[, Pathway.name := factor(Pathway.name, levels = path_count$Pathway.name)]
  path_long[, Pathway.name := factor(Pathway.name, levels = path_count$Pathway.name)]
  
  pdf("plots\\Reactome_count.pdf", width = 12, height = 9)
  g <- ggplot(path_long, aes(y = Pathway.name, x = Regulation.group))
  g <- g + geom_tile(aes(fill = percentage), color = "white", size = 1)
  g <- g + coord_equal()
  g <- g + scale_fill_gradient(high = "#fcb533", low = "#ebf7e9", na.value = "white", limits = c(0, 10), oob = scales::squish)
  g <- g + theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 11),
                 axis.text.y = element_text(size = 16),
                 )
  g <- g + xlab("") + ylab("")
  g <- g + guides(fill = guide_colorbar(title = "% Sites"))
  print(g)
  dev.off()
  
  pdf("plots\\Reactome_pvalue.pdf", width = 12, height = 9)
  g <- ggplot(df_tests, aes(comparison, Pathway.name))
  g <- g + geom_tile(aes(fill = p.adj), color = "white", size = 1)
  g <- g + coord_equal()
  g <- g + scale_fill_gradient(low = "#51baf4", high = "#ebf7e9", na.value = "white", limits = c(0, 0.2), oob = scales::squish)
  g <- g + scale_x_discrete(expand = c(0, 0), position = "top", breaks = levels(df_tests$comparison))
  g <- g + scale_y_discrete(expand = c(0, 0), limits = rev(levels(df_tests$netphorest_group2)))
  g <- g + theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 11),
                 axis.text.y = element_text(size = 16),
  )
  g <- g + xlab("") + ylab("")
  g <- g + guides(fill = guide_colorbar(title = "adj. p-Value"))
  print(g)
  dev.off()
  
 
  })
  