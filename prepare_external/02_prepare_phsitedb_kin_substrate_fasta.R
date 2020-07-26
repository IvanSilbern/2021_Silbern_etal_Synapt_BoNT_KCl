# DO:
# Prepare Kinase-Substrate data set from Phosphopeptide plus database
# Save Sequence windows as separate fasta file "temp\\Phsitedb_kin_sub.fasta"

# INPUT:
# "external\\data_Phosphositeplus\\Kinase_Substrate_Dataset.txt"

# OUTPUT:
# "temp\\Phsitedb_kin_sub.fasta"


local({
  
  library(data.table)
  
  # PhosphositePlus data base
  kin_sub <- fread("external\\data_Phosphositeplus\\Kinase_Substrate_Dataset.txt", skip = 2, check.names = TRUE)
  
  # create kin_substrate id
  kin_sub[, id := 1:.N]
  
  # sequence windows without space
  kin_sub[, SEQ_WIND_NOSPACE := toupper(gsub("_", "", SITE_...7_AA))]
  
  # substrates per sequence window and kinase-substrate interactions ids
  kin_sub_sw <- kin_sub[, lapply(.SD, paste, collapse = "/"), by = SEQ_WIND_NOSPACE, .SDcols = c("id", "SUB_GENE")]
  
  # prepare fasta file from sequence windows
  sw_fasta <- vector("character", length = nrow(kin_sub_sw)*2)
  sw_fasta[seq(1, length(sw_fasta), by = 2)] <- paste0(">", "db|", gsub(" ", "", kin_sub_sw$SUB_GENE), "|", kin_sub_sw$id)
  sw_fasta[seq(2, length(sw_fasta), by = 2)] <- kin_sub_sw$SEQ_WIND_NOSPACE
  
  if(!dir.exists("temp")) dir.create("temp")
  writeLines(sw_fasta, "temp\\Phsitedb_kin_sub.fasta")
  
})