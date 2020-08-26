# DO:
# add uniprot information about protein length and function
# (use Uniprot id mapping service to retrieve the information)
# Retrieve functional information for rat, mouse and human proteins
# use information for rat protein in the first line, supplement missing
# annotation by annotation of respective mouse and human proteins (if available)
# add manual annotation of key-words associated with protein function/localization
#
# INPUT:
# "temp\\PhPeptIntensities2.tsv"
# "external\\data_Uniprot\\PhPeptIntensities2_annotation_upID_rat.tsv"
# "external\\data_Uniprot\\UP_Genes_Function_Rat.tsv"
# "external\\data_Uniprot\\UP_Genes_Function_Mouse.tsv"
# "external\\data_Uniprot\\UP_Genes_Function_Human.tsv"
#
# OUTPUT:
# "temp\\PhPeptIntensities3.tsv"


local({
  
  library(data.table)
  library(stringr)
    
  df <- fread("temp\\PhPeptIntensities2.tsv")
  annot_rat <- fread("external\\data_Uniprot\\PhPeptIntensities2_annotation_upID_rat.tsv", check.names = TRUE)
  annot_man <- fread("external\\manual_annotation_proteins.txt")
  names(annot_man)[names(annot_man) == "Groups"] <- "Keyword_manual"
  
  # add information about protein length
  # based on uniprot id
  
  df <- merge(df, annot_rat[, c("Entry", "Length")], by.x = "Accession", by.y = "Entry", all.x = TRUE)
  sum(is.na(df$Length))
  df[is.na(df$Length), c("Accession", "Gene.name", "Length")]
  df[Accession == "A0A0G2K2D6", Length := 218]
  df[Accession == "A0A0G2JWA9", Length := 2055]
  
  writeLines(unique(df$Gene.name), "temp\\Gene_names_1stClass.txt") # use for uniprot id mapping

  # Add uniprot annotation
  annot_rat <- fread("external\\data_Uniprot\\UP_Genes_Function_Rat.tsv", check.names = TRUE)  
  annot_mouse <- fread("external\\data_Uniprot\\UP_Genes_Function_Mouse.tsv", check.names = TRUE)  
  annot_human <- fread("external\\data_Uniprot\\UP_Genes_Function_Human.tsv", check.names = TRUE)  
  
  names(annot_rat) <- gsub("\\.\\.CC\\.", "", names(annot_rat))
  names(annot_mouse) <- gsub("\\.\\.CC\\.", "", names(annot_mouse))
  names(annot_human) <- gsub("\\.\\.CC\\.", "", names(annot_human))
  
  # extract protein scores
  annot_rat[, SCORE := as.integer(str_match(Annotation, "^([0-5])")[, 2])]
  annot_mouse[, SCORE := as.integer(str_match(Annotation, "^([0-5])")[, 2])]
  annot_human[, SCORE := as.integer(str_match(Annotation, "^([0-5])")[, 2])]
  
  # Remove entries without function 
  annot_rat <- annot_rat[Function != ""]
  annot_rat <- annot_rat[Function != "FUNCTION: Its physiological role is not yet clear."]
  
  annot_mouse <- annot_mouse[Function != ""]
  annot_mouse <- annot_mouse[Function != "FUNCTION: Its physiological role is not yet clear."]
  
  annot_human <- annot_human[Function != ""]
  annot_human <- annot_human[Function != "FUNCTION: Its physiological role is not yet clear."]
  
  # Remove "Function"
  annot_rat[, Function := gsub("FUNCTION:", "", Function)]
  annot_mouse[, Function := gsub("FUNCTION:", "", Function)]
  annot_human[, Function := gsub("FUNCTION:", "", Function)]
  
  
  # Expand Gene names
  annot_rat[, Gene.name := str_split(Gene.names, " ")]
  temp <- annot_rat[, list(Gene.name = unlist(Gene.name)), by = Entry]
  annot_rat <- merge(temp, annot_rat[, -c("Gene.name")], by = c("Entry"), all.x = TRUE)
  dim(annot_rat)
  
  annot_mouse[, Gene.name := str_split(Gene.names, " ")]
  temp <- annot_mouse[, list(Gene.name = unlist(Gene.name)), by = Entry]
  annot_mouse <- merge(temp, annot_mouse[, -c("Gene.name")], by = c("Entry"), all.x = TRUE)
  dim(annot_mouse)
  
  annot_human[, Gene.name := str_split(Gene.names, " ")]
  temp <- annot_human[, list(Gene.name = unlist(Gene.name)), by = Entry]
  annot_human <- merge(temp, annot_human[, -c("Gene.name")], by = c("Entry"), all.x = TRUE)
  dim(annot_human)
  
  # sort based on the score
  annot_rat   <- annot_rat[order(-SCORE)]
  annot_mouse <- annot_mouse[order(-SCORE)]
  annot_human <- annot_human[order(-SCORE)]

  # remove duplicated entries (gene name)
  annot_rat   <- annot_rat[!duplicated(annot_rat$Gene.name)]
  annot_mouse <- annot_mouse[!duplicated(Gene.name)]
  annot_human <- annot_human[!duplicated(Gene.name)]
    
  # add functional annotation to the cand tables
  df[, Function := annot_rat$Function[match(df$Gene.name, annot_rat$Gene.name)]]
  dim(df)
  sum(is.na(df$Function))
  
  df[is.na(df$Function), Function := annot_mouse$Function[match(df$Gene.name[is.na(df$Function)], annot_mouse$Gene.name)]]
  df[is.na(df$Function), Function := annot_human$Function[match(df$Gene.name[is.na(df$Function)], annot_human$Gene.name)]]
  
  # add manual keyword annotation
  df <- merge(df, annot_man, by = "Gene.name", all.x = TRUE)
  
  # correct annotation
  df[df$Gene.name == "Pcp2", Function := c("May function as a cell-type specific modulator for G protein-mediated cell signaling")]
  fwrite(df, "temp\\PhPeptIntensities3.tsv", sep = "\t")

})