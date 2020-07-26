# DO:
# Prepare Kinase-Substrate data set from Phosphopeptide plus database
# Add Stringids based on uniprot mapping and aliases from stringid database
# save Kinase and substrate Stringids in separate files: "temp\\kinase_Stringid.tsv" and "temp\\substrate_Stringid.tsv"

#INPUT:
# "external\\data_Phosphositeplus\\Kinase_Substrate_Dataset.txt
# "external\\data_Uniprot\\MappingTable_rat_Stringid.txt"
# "external\\data_Uniprot\\UP_AnnotationScore.tsv"
# "external\\data_Stringdb\\10116.protein.aliases.v10.5.txt.bz2"

#OUTPUT:
# "temp\\kinase_Stringid.tsv"
# "temp\\substrate_Stringid.tsv"

local({
  
  if(!dir.exists("temp")) dir.create("temp")
 
  library(data.table)
  library(stringr)
  
  # PhosphositePlus data base
  kin_sub <- fread("external\\data_Phosphositeplus\\Kinase_Substrate_Dataset.txt", skip = 2, check.names = TRUE)
  
  # UP mapping table
  up_map <- fread("external\\data_Uniprot\\MappingTable_rat_Stringid.txt", check.names = TRUE)

  # up annotation scores
  up_score <- fread("external\\data_Uniprot\\UP_AnnotationScore.tsv", check.names = TRUE)
  up_score <- up_score[, SCORE2 := str_match(SCORE, "^([0-5])\\s")[, 2]]

  # add up score information
  up_map <- merge(up_map, up_score[, c("UNIPROTKB", "SCORE2")], by.x = "Entry", by.y = "UNIPROTKB", all.x = TRUE)
  up_map <- up_map[order(-(up_map$Status == "reviewed"), -up_map$SCORE2, -Entry)]

  # string id aliases
  alias <- fread("external\\data_Stringdb\\10116.protein.aliases.v10.5.txt.bz2", check.names = TRUE, skip = 1, header = FALSE)
  names(alias) <- c("string_protein_id", "alias", "source")
  setkey(alias, "alias")

  # Phosphosite Plus Kinase-Substrate Data Set
  kin <- unique(toupper(kin_sub$GENE))
  kin <- kin[kin != ""]
  length(kin)

  sub <- unique(toupper(kin_sub$SUB_GENE))
  length(sub)

  # convert kin genes to stringid first based on uniprot-mapping
  temp <- up_map[Gene.name.main.up %in% kin]
  temp <- temp[!is.na(Stringid)]
  temp <- temp[!duplicated(temp$Gene.name.main.up)]
  kin_ids <- temp$Stringid[match(kin, temp$Gene.name.main.up)]
  sum(is.na(kin_ids))

  # convert remaining kinase genes to string ids based on the alias mapping provided by string database
  temp <- alias[kin[is.na(kin_ids)]]
  temp <- temp[!duplicated(temp$alias)]
  kin_ids[is.na(kin_ids)] <- temp[match(kin[is.na(kin_ids)], temp$alias)]$string_protein_id 
  #kin[is.na(kin_ids)]
  kin_ids[kin == "CDK8"] <- alias[grepl("CDK8", toupper(alias)), "string_protein_id"]

  length(kin_ids)
  sum(is.na(kin_ids))

  kin     <- kin[!is.na(kin_ids)]
  kin_ids <- kin_ids[!is.na(kin_ids)]
  sum(duplicated(kin_ids))

  
  # find stringids for substrates
  temp <- up_map[Gene.name.main.up %in% sub]
  temp <- temp[!is.na(Stringid)]
  temp <- temp[!Stringid == ""]
  temp <- temp[!duplicated(temp$Gene.name.main.up)]
  sub_ids <- temp$Stringid[match(sub, temp$Gene.name.main.up)]
  sum(is.na(sub_ids))
  
  temp <- alias[sub[is.na(sub_ids)]]
  temp <- temp[!duplicated(temp$alias)]
  sub_ids[is.na(sub_ids)] <- temp[match(sub[is.na(sub_ids)], temp$alias)]$string_protein_id 
  #sub[is.na(sub_ids)]
  
  sub_ids[sub == "NHE-3"]  <- "10116.ENSRNOP00000020711"

  # phosphorylated amino acid
  aa <- unlist(
    
    lapply(
      
      strsplit(kin_sub$SITE_...7_AA, split = ""),
      
      "[", 8)
    
  )

  kin_sub[, Amino.acid := aa]
  
  # assess kinase specificity
  kin_specificity <- data.table()
  for(i in seq_along(kin)){
    
    aa <- kin_sub$Amino.acid[grepl(kin[i], toupper(kin_sub$GENE))]
    
    kin_specificity <- rbind(kin_specificity,
                            data.table(kin = kin[i],
                                             s = sum(aa == "s")*100/length(aa),
                                             t = sum(aa == "t")*100/length(aa),
                                             y = sum(aa == "y")*100/length(aa),
                                             stringsAsFactors = FALSE
                            ))
    
    
  }

  # View(kin_specificity)
  
  kin_specificity$st <- unlist(Map(sum, kin_specificity$s, kin_specificity$t, na.rm = T))
  kin_specificity$st_specific  <- kin_specificity$st > 60
  kin_specificity$y_specific   <- kin_specificity$y  > 60 
  kin_specificity$y_ambiguous  <- kin_specificity$y > 0 & kin_specificity$y < 60 & kin_specificity$st < 60
  
  sum(kin_specificity$st_specific, na.rm = T)
  sum(kin_specificity$y_specific, na.rm = T)
  sum(kin_specificity$y_ambiguous, na.rm = T)


  #kinase gene names + stringids + kinase specificity
  phsp_kin <- data.table(kin = kin, kin_ids = kin_ids, stringsAsFactors = FALSE)
  phsp_kin <- merge(phsp_kin, kin_specificity, by = "kin")
  fwrite(phsp_kin, "temp\\kinase_Stringid.tsv", sep = "\t")

  # substrate gene names + stringids
  phsp_sub <- data.table(sub = sub, sub_ids = sub_ids, stringsAsFactors = FALSE)
  fwrite(phsp_sub, "temp\\substrate_Stringid.tsv", sep = "\t")

})