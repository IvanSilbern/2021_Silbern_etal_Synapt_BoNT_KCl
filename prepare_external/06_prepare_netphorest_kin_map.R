
# DO:
# use netphorest human protein name map
# convert kinase gene names into a suitable format
# save as "temp\\np_genes.txt"
# use Uniprot conversion tool (externally) to get UniProtKB for Homo sapience [9606]


# INPUT:
# "external\\data_NetworKIN\\group_human_protein_name_map.txt" (can be found in "NetworKIN_release3.0\data")
# "external\\data_Uniprot\\np_up_gene_conversion.txt" (output of the Uniprot id conversion tool)
#

# OUTPUT:
# "temp\\np_genes.txt" (list of gene names, used for Uniprot id conversion)
# "temp\\netphorest_groups_all.tsv" (inculdes duplicated gene names; long format for netphorest groups)
# "temp\\netphorest_groups.tsv" (duplicated genes removed, netphorest kinase tree is collapsed) 

local({

  if(!dir.exists("temp")) dir.create("temp")
  
  library(data.table)
  library(stringr)

  np <- fread("external\\data_NetworKIN\\group_human_protein_name_map.tsv", header = FALSE, fill = TRUE)
  
  names(np) <- c("Tree", "Group", "Gene_id")
  np <- np[np$Tree == "KIN"]

  # add identifier that can be used as a rank
  # since the table is sorted
  # more general groups will have lower id
  np[, id := 1:.N]

  # # example
  # np[np$Gene_id == "MAPK4"]

  # remove "leaves":
  np <- np[grepl("_group", np$Group)]

  # convert Gene_id to a common Gene symbol
  np[, Gene := toupper(Gene_id)]

  np[, Gene := gsub("^PKA", "PRKAC", Gene)]
  np[, Gene := gsub("PKBALPHA$", "AKT1", Gene)]
  np[, Gene := gsub("PKBBETA$",  "AKT2", Gene)]
  np[, Gene := gsub("PKBGAMMA$", "AKT3", Gene)]
  np[, Gene := gsub("^PKC",   "PRKC", Gene)]
  np[, Gene := gsub("^PKD",   "PRKD", Gene)]
  np[, Gene := gsub("^PKG",   "PRKG", Gene)]
  np[, Gene := gsub("CGK[I]+",   "", Gene)]
  
  
  np[, Gene := gsub("(.)ALPHA",   "\\1A", Gene)]
  np[, Gene := gsub("(.)BETA",    "\\1B", Gene)]
  np[, Gene := gsub("(.)GAMMA",   "\\1G", Gene)]
  np[, Gene := gsub("(.)DELTA",   "\\1D", Gene)]
  np[, Gene := gsub("(.)EPSILON", "\\1E", Gene)]
  np[, Gene := gsub("(.)ZETA",    "\\1Z", Gene)]
  np[, Gene := gsub("(.)THETA",   "\\1Q", Gene)]
  np[, Gene := gsub("(.)IOTA",    "\\1I", Gene)]
  np[, Gene := gsub("(.)ETA",     "\\1H", Gene)]
  
  np[, Gene := gsub("ALPHAK([1-3])", "ALPK\\1", Gene)]
  np[, Gene := gsub("AMPKA([1-9])", "PRKAA\\1", Gene)]
  np[, Gene := gsub("ANKRD3", "RIPK4", Gene)]
  np[, Gene := gsub("ANPA", "NPR1", Gene)]
  np[, Gene := gsub("ANPB", "NPR2", Gene)]
  np[, Gene := gsub("ADCK3", "COQ8A", Gene)]
  np[, Gene := gsub("ADCK4", "COQ8B", Gene)]
  np[, Gene := gsub("ACTR2B", "ACVR2B", Gene)]
  np[, Gene := gsub("ALK2", "ACVR1", Gene)]
  np[, Gene := gsub("^AURORA", "AURK", Gene)]
  np[, Gene := gsub("^AURORA", "AURK", Gene)]
  np[, Gene := gsub("CAMKIV", "CAMK4", Gene)]
  np[, Gene := gsub("CAMKII", "CAMK2", Gene)]
  np[, Gene := gsub("CAMKIA", "CAMK1", Gene)]
  np[, Gene := gsub("CAMKIB", "PNCK", Gene)]
  np[, Gene := gsub("CAMKIG", "CAMK1G", Gene)]
  
  np[, Gene := gsub("CAMLCK", "MYLK3", Gene)]
  np[, Gene := gsub("CK1A$", "CSNK1A1", Gene)]
  np[, Gene := gsub("CK1A2$", "CSNK1A1L", Gene)]
  np[, Gene := gsub("CK1D$", "CSNK1D", Gene)]
  np[, Gene := gsub("CK1G1$", "CSNK1G1", Gene)]
  np[, Gene := gsub("CK1G3$", "CSNK1G3", Gene)]
  np[, Gene := gsub("CK2A$", "CSNK2A1", Gene)]
  np[, Gene := gsub("CYGD$", "GUCY2D", Gene)]
  np[, Gene := gsub("CYGF$", "GUCY2F", Gene)]
  np[, Gene := gsub("DMPK1$", "DMPK", Gene)]
  np[, Gene := gsub("DNAPK$", "PRKDC", Gene)]
  np[, Gene := gsub("DOMAIN2_GCN2$", "EIF2AK4", Gene)]
  np[, Gene := gsub("DOMAIN2MSK1$", "RPS6KA5", Gene)]
  np[, Gene := gsub("DOMAIN2MSK2$", "RPS6KA4", Gene)]
  np[, Gene := gsub("DOMAIN2OBSCN$", "OBSCN", Gene)]
  np[, Gene := gsub("DOMAIN2RSK1$", "RPS6KA1", Gene)]
  np[, Gene := gsub("DOMAIN2RSK2$", "RPS6KA3", Gene)]
  np[, Gene := gsub("DOMAIN2RSK3$", "RPS6KA2", Gene)]
  np[, Gene := gsub("DOMAIN2RSK4$", "RPS6KA6", Gene)]
  np[, Gene := gsub("DOMAIN2SPEG$", "PAEP", Gene)]
  np[, Gene := gsub("FUSED$", "STK36", Gene)]
  np[, Gene := gsub("HH498$", "TNNI3K", Gene)]
  np[, Gene := gsub("HSER$", "GUCY2C", Gene)]
  np[, Gene := gsub("JAKB1$", "JAK1", Gene)]
  np[, Gene := gsub("JAKB2$", "JAK2", Gene)]
  np[, Gene := gsub("JAKB3$", "JAK3", Gene)]
  np[, Gene := gsub("KHS2$", "MAP4K3", Gene)]
  np[, Gene := gsub("LMR3$", "LMTK3", Gene)]
  np[, Gene := gsub("MRCKA$", "CDC42BPA", Gene)]
  np[, Gene := gsub("MRCKB$", "CDC42BPB", Gene)]
  np[, Gene := gsub("P70S6K$", "RPS6KB1", Gene)]
  np[, Gene := gsub("P70S6KB$", "RPS6KB2", Gene)]
  np[, Gene := gsub("PFTAIRE1$", "CDK14", Gene)]
  np[, Gene := gsub("PFTAIRE2$", "CDK15", Gene)]
  np[, Gene := gsub("PITSLRE$", "CDK11B", Gene)]
  np[, Gene := gsub("RSKL1$", "RPS6KC1", Gene)]
  np[, Gene := gsub("RSKL2$", "RPS6KL1", Gene)]
  np[, Gene := gsub("SBK$", "SBK1", Gene)]
  np[, Gene := gsub("SGK269$", "PEAK1", Gene)]
  np[, Gene := gsub("SKMLCK$", "MYLK2", Gene)]
  np[, Gene := gsub("SMMLCK$", "MYLK", Gene)]
  np[, Gene := gsub("STLK3$", "STK39", Gene)]
  np[, Gene := gsub("STLK5$", "STRADA", Gene)]
  np[, Gene := gsub("STLK6$", "STRADB", Gene)]
  np[, Gene := gsub("SURTK106$", "STYK1", Gene)]
  np[, Gene := gsub("TAO([1-3])$", "TAOK\\1", Gene)]
  np[, Gene := gsub("TYK2B$", "TYK2", Gene)]
  np[, Gene := gsub("VACAMKL$", "CAMKV", Gene)]
  
  writeLines(unique(np$Gene), "temp\\np_genes.txt")

  conv <- fread("external\\data_Uniprot\\np_up_gene_conversion.txt", check.names = TRUE)
  names(conv)[1] <- "Gene"
  
  conv[, Gene.name.main := unlist(lapply(str_split(Gene.names, " "), "[", 1))]
  conv <- conv[order(Status)]
  conv[, id := 1:.N]
  
  # several gene names are possible
  temp <- conv[, list(Gene.names = unlist(str_split(Gene.names, " "))), by = id]
  conv <- merge(conv[, -c("Gene.names")], temp, by = "id", all = TRUE)
  conv <- conv[order(id)]
  conv <- conv[!duplicated(conv[, c("Gene", "Gene.names")])]
  conv[, id2 := 1:.N]
  
  temp <- conv[, list(Gene = unlist(str_split(Gene, ","))), by = id2]
  conv <- merge(conv[, -c("Gene")], temp, by = "id2", all = TRUE)
  conv <- conv[order(id)]
  
  # add to netphorest kinase mapping
  np <- merge(np, conv[, -c("id", "id2")], by = "Gene", all = TRUE, allow.cartesian = TRUE)
  dim(np)
  
  # custom changes
  add_row <- np[Gene == "DLK"][1]
  add_row[, Gene.names := "MAP3K12"]
  np <- rbind(np, add_row)
  
  add_row <- np[Gene == "CSNK2A1"][1]
  add_row[, Gene.names := "CSNK2B"]
  np <- rbind(np, add_row)
  
  add_row <- np[Gene == "PHKG1"][1]
  add_row[, Gene.names := "PHKA1"]
  np <- rbind(np, add_row)
  
  add_row <- np[Gene == "PHKG1"][1]
  add_row[, Gene.names := "PHKB"]
  np <- rbind(np, add_row)
  
  add_row <- np[Gene == "CK1E"][1]
  add_row[, Gene.names := "CSNK1E"]
  np <- rbind(np, add_row)
  
  # sort np from higher id (more specific group) to lower id (more genral group)
  np <- np[order(-id)]
  
  # Remove ambiguos annotations due to overlapping gene synonyms
  np <- np[!(np$Gene == "NEK8" & np$Gene.names == "NEK9")]
  np <- np[!(np$Gene == "NEK9" & np$Gene.names == "NEK8")]
  np <- np[!(np$Gene == "PAK5" & np$Gene.names == "PAK6")]
  np <- np[!(np$Gene == "PAK6" & np$Gene.names == "PAK5")]
  np <- np[!(np$Gene == "PKN1" & np$Gene.names == "PAK1")]
  np <- np[!(np$Gene == "PAK1" & np$Gene.names == "PKN1")]
  np <- np[!(np$Gene == "RPS6KA3" & np$Gene.names == "RPS6KA3")]
  np <- np[!(np$Gene == "RPS6KA1" & np$Gene.names == "RPS6KA1")]
  np <- np[!(np$Gene == "PRKCA" & np$Gene.names == "PRKACA")]
  np <- np[!(np$Gene == "PRKACA" & np$Gene.names == "PRKCA")]
  np <- np[!(np$Gene == "PDK1" & np$Gene.names == "PDK1")]
  np <- np[!(np$Gene == "PDK1" & np$Gene.names == "PDHK1")]
  np <- np[!(np$Gene == "MST3" & np$Gene.names == "STK3")]
  np <- np[!(np$Gene == "MAK" & np$Gene.names == "ALPK3")]
  np <- np[!(np$Gene %in% c("CDK7", "DDR1") & np$Gene.names == "CAK")]
  np <- np[!(np$Gene %in% c("CAMK4", "CAMK2G") & np$Gene.names == "CAMK")]
  
  fwrite(np, "temp\\netphorest_groups_all.tsv", sep = "\t")
  
  # check duplicate pair Gene - Gene.names
  dim(np)
  sum(duplicated(np[, c("Gene", "Group", "Gene.names")]))
  np <- np[!duplicated(np[, c("Gene", "Group", "Gene.names")])]
  np[is.na(Gene.name.main), Gene.names := "RSK5"]
  np[is.na(Gene.name.main), Gene.name.main := "RSK5"]
  
  np_red <- np[!duplicated(np[, c("Gene.name.main", "Group")])]
  tree <- np_red[, list(id = list(id), Group = list(c(Gene.name.main, Group))), by = "Gene.name.main"]
  np <- merge(np, tree, by = "Gene.name.main", suffixes = c("", "_all"))
  np <- np[!duplicated(np[, c("Gene", "Gene.names")])]
  
  fwrite(np, "temp\\netphorest_groups.tsv", sep = "\t", sep2 = c("", ";", ""))
  
})