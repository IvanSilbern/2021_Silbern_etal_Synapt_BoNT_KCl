# Do: make a graph of phosphatase interactions

# INPUT:
# "temp\\PhPeptIntensities6.tsv"
# "temp\\Protein_classification.txt"
# "temp\\Stringids_Rat_SynaptosomalProteins.txt"
# "external\\data_Uniprot\\MappingTable_rat_Stringid.txt"
# "temp\\stringDB_graph.Rdata"

# OUTPUT:
# "plots\\Graph_PP1_substrates_endo_exocyt.pdf"
# "Figures\\SupplFig_7\\PP1_graph.Rdata"
# "Figures\\SupplFig_7\\PP1_slices.Rdata"
# "plots\\Graph_CaN_substrates_endo_exocyt.pdf"
# "Figures\\SupplFig_7\\CN_graph.Rdata"
# "Figures\\SupplFig_7\\CN_slices.Rdata"


local({
  
  library(ggplot2)
  library(igraph)
  library(data.table)
  library(stringr)
  
  if(!dir.exists("Figures\\SupplFig_7")) dir.create("Figures\\SupplFig_7", recursive = TRUE)
  
  # phosphosites data
  df <- fread("temp\\PhPeptIntensities6.tsv")
  df <- df[Regulation_group_resolved != "not-affected"]
  
  # protein classification
  prot_class <- fread("temp\\Protein_classification.tsv")
  
  # Load stringids of proteins found in Synaptosomes
  prot_ids <- readLines("temp\\Stringids_Rat_SynaptosomalProteins.txt")
  
  # UP mapping table
  up_map <- fread("external\\data_Uniprot\\MappingTable_rat_Stringid.txt", check.names = TRUE)
  
  # load string db graph
  load("temp\\stringDB_graph.Rdata")
  
  # keep only vertices that were found in rat synaptosomes
  graph <- delete_vertices(graph, v = names(V(graph))[!names(V(graph)) %in% prot_ids])
  
  # putative pp1 phosphatase substrates
  
  # select ids of regulated proteins
  cand_genes <- unique(df[PP1_substr == TRUE & Gene.name %in% prot_class[grepl("Exocytosis", Keyword_manual) | grepl("Endocytosis", Keyword_manual)]$Gene.name]$Gene.name)
  cand_id <- unique(df[PP1_substr == TRUE & Gene.name %in% prot_class[grepl("Exocytosis", Keyword_manual) | grepl("Endocytosis", Keyword_manual)]$Gene.name]$Stringid)
  
  ph_genes <- c("PPP1CA", "PPP1CB", "PPP1CC")
  ph_id <- c("10116.ENSRNOP00000025282", "10116.ENSRNOP00000006190", "10116.ENSRNOP00000047912")
  
  sel_id <- unique(c(cand_id, ph_id))
  
  sp <- shortest_paths(graph, from = ph_id, to = c(cand_id[cand_id %in% as_ids(V(graph))], ph_id), output = "both")
  sg <- do.call(igraph::union, sp$epath)
  sg <- matrix(unlist(str_split(as_ids(sg), "\\|")), ncol = 2, byrow = TRUE)
  sg <- graph_from_edgelist(sg, directed = FALSE)
  sg <- simplify(sg)
  
  V(sg)$Gene.name <- df$Gene.name[match(as_ids(V(sg)), df$Stringid)]
  V(sg)$Gene.name[is.na(V(sg)$Gene.name)] <- up_map$Gene.name.main[match(as_ids(V(sg))[is.na(V(sg)$Gene.name)], up_map$Stringid)]
  
  V(sg)$col <- "grey"
  V(sg)$col[V(sg)$Gene.name %in% prot_class$Gene.name] <- "cyan"
  V(sg)$frame_color <- "black"
  V(sg)$frame_color[V(sg)$Gene.name %in% cand_genes] <- "red"
  V(sg)$col[as_ids(V(sg)) %in% ph_id]  <- "red"
  V(sg)$shape <- "pie"
  V(sg)$shape[V(sg)$col == "grey"] <- "circle"
  V(sg)$shape[as_ids(V(sg)) %in% ph_id] <- "square"
  V(sg)$size <- 8
  V(sg)$size[V(sg)$col == "grey"] <- 4
  
  V(sg)$Label <- ""
  V(sg)$Label[V(sg)$col != "grey"] <- V(sg)$Gene.name[V(sg)$col != "grey"]
  V(sg)$Label[degree(sg) > 2] <- V(sg)$Gene.name[degree(sg) > 2]
  
  slices <- vector("list", length = length(V(sg)))
  for(i in seq_along(V(sg)$Gene.name)){
    
    gene <- V(sg)$Gene.name[i]
    temp <- unlist(prot_class[Gene.name == gene, .SD, .SDcols = c("Ca_dependent", "Cycling_dependent")])
    
    # replace NA with 0 
    temp[is.na(temp)] <- 0
    
    # if there are no entry -> use NA's
    if(length(temp) == 0) temp <- c(NA, NA)
    
    slices[[i]] <- temp
    next
    
  }
  
  pdf("plots\\Graph_PP1_substrates_endo_exocyt.pdf", width = 5, height = 5)
  plot(sg,
       vertex.label = V(sg)$Label,
       vertex.label.cex = 0.6,
       vertex.color = V(sg)$col,
       vertex.frame.color = V(sg)$frame_color,
       vertex.size = V(sg)$size,
       vertex.label = "",
       vertex.shape = V(sg)$shape,
       vertex.pie = slices,
       vertex.pie.color = list(c(scales::alpha("#fcb533", 0.8), scales::alpha("#51baf4", 0.8))))
  dev.off()
  
  # save figure data
  save(sg,     file = "Figures\\SupplFig_7\\PP1_graph.Rdata")
  save(slices, file = "Figures\\SupplFig_7\\PP1_slices.Rdata")
  
  # putative calcineurin phosphatase substrates
  
  # select ids of regulated proteins
  cand_genes <- unique(df[CN_substr == TRUE & Gene.name %in% prot_class[grepl("Exocytosis", Keyword_manual) | grepl("Endocytosis", Keyword_manual)]$Gene.name]$Gene.name)
  cand_id <- unique(df[CN_substr == TRUE & Gene.name %in% prot_class[grepl("Exocytosis", Keyword_manual) | grepl("Endocytosis", Keyword_manual)]$Gene.name]$Stringid)
  
  ph_genes <- c("PPP3CA", "PPP3CB", "PPP3CC")
  ph_id <- c("10116.ENSRNOP00000049674", "10116.ENSRNOP00000010476", "10116.ENSRNOP00000012993")
  ph_id %in% prot_ids
  
  sel_id <- unique(c(cand_id, ph_id))
  
  sp <- shortest_paths(graph, from = ph_id[ph_id %in% as_ids(V(graph))], to = c(cand_id[cand_id %in% as_ids(V(graph))], ph_id), output = "both")
  sg <- do.call(igraph::union, sp$epath)
  sg <- matrix(unlist(str_split(as_ids(sg), "\\|")), ncol = 2, byrow = TRUE)
  sg <- graph_from_edgelist(sg, directed = FALSE)
  #sg <- simplify(sg)
  
  V(sg)$Gene.name <- df$Gene.name[match(as_ids(V(sg)), df$Stringid)]
  V(sg)$Gene.name[is.na(V(sg)$Gene.name)] <- up_map$Gene.name.main[match(as_ids(V(sg))[is.na(V(sg)$Gene.name)], up_map$Stringid)]
  
  V(sg)$col <- "grey"
  V(sg)$col[V(sg)$Gene.name %in% prot_class$Gene.name] <- "cyan"
  V(sg)$frame_color <- "black"
  V(sg)$frame_color[V(sg)$Gene.name %in% cand_genes] <- "red"
  V(sg)$col[as_ids(V(sg)) %in% ph_id]  <- "red"
  V(sg)$shape <- "pie"
  V(sg)$shape[V(sg)$col == "grey"] <- "circle"
  V(sg)$shape[as_ids(V(sg)) %in% ph_id] <- "square"
  V(sg)$size <- 8
  V(sg)$size[V(sg)$col == "grey"] <- 4
  
  V(sg)$Label <- ""
  V(sg)$Label[V(sg)$col != "grey"] <- V(sg)$Gene.name[V(sg)$col != "grey"]
  V(sg)$Label[degree(sg) > 2] <- V(sg)$Gene.name[degree(sg) > 2]
  
  slices <- vector("list", length = length(V(sg)))
  for(i in seq_along(V(sg)$Gene.name)){
    
    gene <- V(sg)$Gene.name[i]
    temp <- unlist(prot_class[Gene.name == gene, .SD, .SDcols = c("Ca_dependent", "Cycling_dependent")])
    
    # replace NA with 0 
    temp[is.na(temp)] <- 0
    
    # if there are no entry -> use NA's
    if(length(temp) == 0) temp <- c(NA, NA)
    
    slices[[i]] <- temp
    next
    
  }
  
  pdf("plots\\Graph_CaN_substrates_endo_exocyt.pdf", width = 5, height = 5)
  plot(sg,
       vertex.label = V(sg)$Label,
       vertex.label.cex = 0.6,
       vertex.color = V(sg)$col,
       vertex.frame.color = V(sg)$frame_color,
       vertex.size = V(sg)$size,
       vertex.label = "",
       vertex.shape = V(sg)$shape,
       vertex.pie = slices,
       vertex.pie.color = list(c(scales::alpha("#fcb533", 0.8), scales::alpha("#51baf4", 0.8))))
  dev.off()
  
  # save figure data
  save(sg,     file = "Figures\\SupplFig_7\\CN_graph.Rdata")
  save(slices, file = "Figures\\SupplFig_7\\CN_slices.Rdata")
  
})