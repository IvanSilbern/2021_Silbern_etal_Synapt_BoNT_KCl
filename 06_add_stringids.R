# DO:
# Add Stringids to phosphosites data set
# a. using uniprot mapping
# b. using stringid alias file (for stringids missing after a.)
# check if a stringid belongs to a putative kinase
# combine unique string id from phosphosites data set and 
# in-house data base of synaptosomal proteins

# INPUT:
# "temp\\Phosphosites_prepared.tsv (Phosphosites data set)
# "external\\data_Uniprot\\MappingTable_rat_Stringid.txt (Uniprot Stringid mapping table)
# "kinase_Stringid.tsv" (extracted kinases from kinase-substrate data set (phosphosite plus db))
# "external\\Synaptosomes_PG_identified.tsv" (synaptosomal proteins, in-house dataset)

# OUTPUT:
# updated "temp\\Phosphosites_prepared.tsv"
# "temp\\Stringids_Rat_SynaptosomalProteins.txt" (list of protein Stringids in Synaptosomes)

#### Load data

local({
  
  if(!dir.exists("temp")) dir.create("temp")
  
  library(data.table)
  library(stringr)

  # phospho data
  ph <- fread("temp\\Phosphosites_prepared.tsv", check.names = TRUE) 

  # UP mapping table
  up_map <- fread("external\\data_Uniprot\\MappingTable_rat_Stringid.txt", check.names = TRUE)

  # kinases based on Phosphosite plus data base
  phsp_kin <- fread("temp\\kinase_Stringid.tsv", check.names = TRUE)

  # full proteome data
  pg_full <- fread("external\\Synaptosomes_PG_identified.tsv", check.names = TRUE)

  # string id aliases
  alias <- fread("external\\data_Stringdb\\10116.protein.aliases.v10.5.txt.bz2", check.names = TRUE, skip = 1, header = FALSE)
  names(alias) <- c("string_protein_id", "alias", "source")
  alias[, alias := toupper(alias)]
  setkey(alias, "alias")

  # map Stringids to the phospho data
  ph[, Stringid := up_map$Stringid[match(toupper(ph$Accession.noIso), up_map$Entry)]]
  sum(is.na(ph$Stringid))
  
  ph[, Gene.name.up := toupper(Gene.name)]

  unmapped <- unique(ph$Gene.name.up[is.na(ph$Stringid) | ph$Stringid == ""])
  unmapped <- unmapped[!(is.na(unmapped) | unmapped == "")]
  unmapped <- toupper(unmapped)
  
  alias_sub <- alias[unmapped]
  alias_sub <- alias_sub[!is.na(string_protein_id)]
  alias_sub <- alias_sub[!duplicated(alias_sub$alias)]
  
  temp <- ph[is.na(ph$Stringid) | ph$Stringid == ""]
  temp <- merge(temp, alias_sub[, c("string_protein_id", "alias")], by.x = "Gene.name.up", by.y = "alias", all.x = TRUE)
  temp[, Stringid := string_protein_id]
  ph <- rbind(ph[!(is.na(ph$Stringid) | ph$Stringid == "")], temp[, -c("string_protein_id")])
  
  # mark putative kinases  
  ph[, Putative.kinase := ph$Stringid %in% phsp_kin$kin_ids]
  length(unique(ph$Stringid[ph$Putative.kinase]))

  fwrite(ph, "temp\\Phosphosites_prepared.tsv", sep = "\t")

  # extract protein ids (stringids) from 
  # in-house database of synaptosomal proteins with
  # identified phosphoproteins
  prot_ids <- unique(c(pg_full$Stringid, ph$Stringid))
  prot_ids <- prot_ids[!is.na(prot_ids)]
  prot_ids <- prot_ids[prot_ids != ""]
  length(prot_ids)
  writeLines(prot_ids, "temp\\Stringids_Rat_SynaptosomalProteins.txt")

})
