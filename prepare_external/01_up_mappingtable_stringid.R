
# DO:
# prepare mapping table
# rat proteins - string database id

# INPUT:
# "external\\data_Uniprot\\MappingTable_rat.tab.bz2"
# "external\\data_Uniprot\\MappingTable_Accession_Stringid_rat.tab.bz2"
# "external\\data_Stringdb\\10116.protein.aliases.v10.5.txt.bz2"

# OUTPUT:
# "external\\data_Uniprot\\MappingTable_rat_Stringid.txt"

local({

  library(data.table)
  library(stringr)
  
  # Load Data
  up_acc_gene   <- fread("external\\data_Uniprot\\MappingTable_rat.tab.bz2", check.names = TRUE)
  up_acc_string <- fread("external\\data_Uniprot\\MappingTable_Accession_Stringid_rat.tab.bz2", check.names = TRUE)
  
  # string id aliases
  alias <- fread("external\\data_Stringdb\\10116.protein.aliases.v10.5.txt.bz2", check.names = TRUE)
  names(alias) <- c("string_protein_id", "alias", "source")
  setkey(alias, "alias")
  
  # Mapping Table for Genes, UP-Accessions and Stringids
  up_map <- merge(up_acc_gene, up_acc_string, by.x = "Entry", by.y = "From", all.x = TRUE)
  names(up_map)[names(up_map) == "To"] <- "Stringid"
  
  # recommended gene name
  up_map[, Gene.name.main := unlist(lapply(str_split(up_map$Gene.names, " "), "[", 1))]
  up_map[!is.na(up_map$Stringid), Mapped_by := "UP"]
  up_map[, Gene.name.main.up := toupper(Gene.name.main)]
  
  unmapped_genes <- unique(up_map$Gene.name.main.up[is.na(up_map$Stringid)])
  unmapped_genes <- unmapped_genes[!is.na(unmapped_genes)]
  unmapped_genes <- unmapped_genes[unmapped_genes != ""]
  alias_sub      <- alias[unmapped_genes]
  alias_sub      <- alias_sub[!duplicated(alias)]
  
  
  temp <- merge(up_map[is.na(up_map$Stringid)],
                alias[, c("alias", "string_protein_id")],
                by.x = "Gene.name.main.up", by.y = "alias",
                all.x = TRUE)
  temp[, Mapped_by := "StringDB"]
  temp[, Stringid := string_protein_id]
  temp <- temp[, -c("string_protein_id")]
  
  up_map <- rbind(up_map[!is.na(up_map$Stringid)], temp)
  dim(up_map)
  
  sum(is.na(up_map$Stringid))
  tail(up_map)
  
  fwrite(up_map, "external\\data_Uniprot\\MappingTable_rat_Stringid.txt", sep = "\t")

})
