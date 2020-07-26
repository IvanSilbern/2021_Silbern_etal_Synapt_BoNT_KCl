# DO:

# Add Networkin predictions to phosphosites data table
# Decide on the consensus kinase based on networkin predictions and 
# Phosphosite plus data base
# Add netphorest groups information

# INPUT:
# "external\\data_NetworKIN\\networking_predictions.tsv.bz2"
# "temp\\Stringids_Rat_SynaptosomalProteins.txt"
# "external\\data_Uniprot\\MappingTable_rat_Stringid.txt"
# "temp\\Human_Stringid_PhSites.tsv"
# "temp\\Human_vs_Rat_Networking_subset.txt"
# "temp\\netphorest_groups.tsv"
# "temp\\Phosphosites_prepared.tsv"
# "temp\\stringDB_graph.Rdata"

# OUTPUT:
# temp\\Phosphosites_prepared2.tsv

local({
    
  if(!dir.exists("temp")) dir.create("temp")
  
  library(data.table)
  library(stringr)
  library(igraph)
  
  # Load Networkin results
  nk <- fread("external\\data_NetworKIN\\networkin_predictions.tsv.bz2", check.names = TRUE)
  
  # Load stringids of proteins found in Synaptosomes
  prot_ids <- readLines("temp\\Stringids_Rat_SynaptosomalProteins.txt")
  
  # UP mapping table
  up_map <- fread("external\\data_Uniprot\\MappingTable_rat_Stringid.txt", check.names = TRUE)
  
  # load human - rat string id mapping table
  Rat_Human_match <- fread("temp\\Human_Stringid_PhSites.tsv", check.names = TRUE)
  blast <- fread("temp\\Human_vs_Rat_Networkin_subset.txt", check.names = TRUE)
  setkey(blast, "Human_Stringid")
  
  # netphorest groups file
  np <- fread("temp\\netphorest_groups.tsv", sep = "\t")
  
  # phosphosites data
  ph <- fread("temp\\Phosphosites_prepared.tsv", check.names = TRUE)
  
  # add information about human Stringid to the main phosphosite data table
  ph <- merge(ph, Rat_Human_match, by.x = c("Protein", "Position", "Amino.acid"), by.y = c("Rat_Protein", "Rat_Position", "Rat_Residue"), all.x = TRUE)
  
  # load string db graph
  load("temp\\stringDB_graph.Rdata")
  
  # keep only vertices that were found in rat synaptosomes
  graph <- delete_vertices(graph, v = names(V(graph))[!names(V(graph)) %in% prot_ids])
  
  # extract kinase and other protein stringids from networkin results
  nk <- nk[nk$tree == "KIN"]
  nk[, string_scores := lapply(nk$string_path, function(x) str_match_all(x, pattern = ",\\s(\\d\\.\\d+)\\s")[[1]][, 2])] # string scores in as a list
  nk[, string_path   := lapply(nk$string_path, function(x) str_match_all(x, pattern = "(ENSP[0-9]+)")[[1]][, 2])]    # stringids as a list without scores
  
  nk[, stringid_kin  := unlist(lapply(string_path, "[", 1))] # exctract kinase ids: the first id in each path
  stringid_all  <- unlist(nk$string_path)                    # all ids in a simple character vector
  
  # convert human string kinase id to rat string kinase id
  nk <- merge(nk, blast[, c("Human_Stringid", "Rat_Stringid")], by.x = "stringid_kin", by.y = "Human_Stringid")
  names(nk)[names(nk) == "Rat_Stringid"] <- "Rat_Stringid_kin"

  # keep nk entries that have kinases found among synaptosomal proteins
  nk <- nk[nk$Rat_Stringid_kin %in% prot_ids]

  # convert Human Stringids in the NetworKIN path to
  nk[, id_ := 1:.N]
  nk[, Rat_string_path := list(list(blast[string_path, Rat_Stringid])), by = "id_"]

  # check how many links are missing in synaptosomal proteome
  nk[, nlinks_missing := unlist(lapply(Rat_string_path, function(x) sum(!x[-c(1, length(x))] %in% prot_ids)))]

  # order nk predictions based on the number of missing link proteins, networkin, netphorest and string scores
  # remove dulplicate entries based on substrat id and position
  nk <- nk[order(-round(networkin_score, 1), nlinks_missing, -netphorest_score, -string_score)]
  nk2 <- nk[!duplicated(nk[, c("X.substrate", "position")])]
  dim(nk2)
  
  # add networkin information to phosphosites data table
  
  ph <- merge(ph, nk2[, c("X.substrate", "position",
                          "Rat_Stringid_kin", "Rat_string_path",
                          "networkin_score", "nlinks_missing", "netphorest_group")], by.x = c("Human_Protein", "Human_Position"), by.y = c("X.substrate", "position"), all.x = TRUE, sort = FALSE)
  
  # convert Rat_Stringid_kin to Gene Symbol using uniprot mapping
  ph[, Rat_Gene_kin := up_map$Gene.name.main[!is.na(up_map$Stringid)][match(Rat_Stringid_kin, up_map$Stringid[!is.na(up_map$Stringid)])]]
  ph[, Rat_Gene_kin := toupper(Rat_Gene_kin)]
  
  # add new column with kinase gene names. Prefer Phosphosites Plus data base over Networkin
  ph[, Kin_mapping  := "NetworKIN"]
  ph[, Kin_gen      := toupper(Rat_Gene_kin)]
  ph[, Kin_stringid := Rat_Stringid_kin]
  
  
  ph[ph$GENE != "", Kin_mapping  := "PhosphositePlus"]
  ph[ph$GENE != "", Kin_gen      := toupper(GENE)]
  ph[ph$GENE != "", Kin_stringid := up_map$Stringid[match(toupper(GENE), toupper(up_map$Gene.name.main))]]
  
  
  # kinases from phosphosite plus db without stringid
  dim(ph)
  ph[is.na(ph$Kin_stringid) & ph$GENE != "", 380:400]
  
  # use networkin annotation for these sites
  ph[is.na(ph$Kin_stringid) & ph$GENE != "" & ph$Rat_Gene_kin != "", Kin_mapping  := "NetworKIN"]
  ph[is.na(ph$Kin_stringid) & ph$GENE != "" & ph$Rat_Gene_kin != "", Kin_gen      := toupper(Rat_Gene_kin)]
  ph[is.na(ph$Kin_stringid) & ph$GENE != "" & ph$Rat_Gene_kin != "", Kin_stringid := Rat_Stringid_kin]
  
  # update netphorest group based on consensus kinase
  # check if there kinases in the data set that match to ambiguous gene names (some overlapping gene synonyms)
  kinases <- unique(ph$Kin_gen)
  kinases[kinases %in% np$Gene.names[duplicated(np$Gene.names)]]
  
  ph <- merge(ph[, -c("netphorest_group")], np[, c("Gene.names", "Group_all")], by.x = "Kin_gen", by.y = "Gene.names", all.x = TRUE, all.y = FALSE, allow.cartesian = FALSE)
  names(ph)[names(ph) == "Group_all"] <- "netphorest_group"
  
  ph[unlist(lapply(Rat_string_path, function(x) length(x) == 0)), Rat_string_path := ""]
  ph[, Rat_string_path := unlist(lapply(Rat_string_path, paste, collapse = ";"))]
  
  save(ph, file = "temp\\ph_prepared.RData")
  
  # update string path for kinase-substrate interactions mapped from Phosphosite-plus data base
  temp <- ph[ph$Kin_mapping == "PhosphositePlus"]
  temp2 <- Map(function(from, to, graph, output){
    
    tryCatch({
      
      paste(names(shortest_paths(graph = graph, from = from, to = to, output = output)$vpath[[1]]), collapse = ";")
      
    }, error = function(e) return(NA))
    
    
  }, from = temp$Kin_stringid, to = temp$Stringid, MoreArgs = list(graph = graph, output = "vpath"))
  
  temp2 <- unlist(temp2)
  writeLines(temp2, "temp\\Rat_path_stringid_from_graph.txt")
  
  ph[ph$Kin_mapping == "PhosphositePlus", Rat_string_path := temp2]
  ph[ph$Rat_string_path == "", Rat_string_path := NA]
  
  fwrite(ph, "temp\\Phosphosites_prepared2.tsv", sep = "\t", sep2 = c("",";",""))

})