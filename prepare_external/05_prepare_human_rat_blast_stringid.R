# DO:
# order BLAST results after matching Human Networking Stringids to Rat protein Data base
# Sorting is based on bitscore -> pident -> gapopen
# Best match for each query sequence is preserved
# Corresponding Rat Stringids are added from Uniprot Stringid mapping table

# INPUT:
# "external\\data_BLAST\\Human_vs_Rat_Networkin.txt.bz2" (BLAST results)
# "external\\data_Uniprot\\MappingTable_rat_Stringid.txt" (Uniprot mapping table)

# OUTPUT:
# "temp\\Human_vs_Rat_Networking_subset.txt"

local({
  
  if(!dir.exists("temp")) dir.create("temp")
  
  library(data.table)
  library(stringr)

  blast <- fread("external\\data_BLAST\\Human_vs_Rat_Networkin.txt.bz2", header = FALSE)
  names(blast) <- c("qacc", "sacc", "evalue", "bitscore", "nident", "pident",
                    "mismatch", "gapopen", "gaps", "qstart", "qend", "sstart",
                    "send", "qlen", "slen", "length", "qseq", "sseq", "qframe", "sframe")
  dim(blast)

  rat_up_map <- fread("external\\data_Uniprot\\MappingTable_rat_Stringid.txt", check.names = TRUE)

  #calculate % coverage of subject sequence (rat protein)
  blast[, sseq_coverage := 100*(send - sstart + 1) / slen]

  # order & keep best matches
  blast  <- blast[order(-sseq_coverage, -pident, -bitscore, gapopen)]

  blast2 <- blast[!duplicated(blast$qacc)]
  blast2 <- blast2[, -c("qseq" ,"sseq")]
  dim(blast2)

  #blast$Human_UP_id    <- extractCapturing("^[a-z][a-z]\\|(.+)\\|", blast$qseqid)
  blast2[, Human_Stringid := str_match(qacc, "^9606\\.(.+)")[, 2]]
  blast2[, Rat_UP_id      := str_match(sacc, "\\|([^|]+)\\|")[, 2]]
  blast2[, Rat_UP_id      := gsub("-\\d+", "", Rat_UP_id)]
  blast2[, Rat_Stringid   := rat_up_map$Stringid[match(toupper(blast2$Rat_UP_id), toupper(rat_up_map$Entry))]]
  sum(is.na(blast2$Rat_Stringid))

  fwrite(blast2, "temp\\Human_vs_Rat_Networkin_subset.txt", sep = "\t")

})