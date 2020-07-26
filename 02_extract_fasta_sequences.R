# DO:
# Prepare fasta file
# Extract fasta protein sequences of phosphoproteins
# Prepare sequence windows in the form of the fasta file

# INPUT:
# "temp\\Phosphosites_prepared.tsv"
# "external\\data_Uniprot\\201902_uniprot-rat.fasta" fasta file

# OUTPUT:
# "Phosphoprotein_sequences.fasta"
# "Sequence_window_7.fasta"
# "Sequence_window_15.fasta"

local({
  
  if(!dir.exists("temp")) dir.create("temp")
  
  library(data.table)
  library(stringr)
  
  # read in fasta file
  fasta.file    <- fread("external\\data_Uniprot\\201902_uniprot-rat.fasta", sep = NULL, header = FALSE)[[1]]
  
  # extract header positions
  ff_header_pos <- which(grepl("^>", fasta.file)) # Lines containing headers 
  ff_headers    <- str_split_fixed(fasta.file[ff_header_pos], pattern = " ", n = 2)[, 1] # extract protein id without description
  ff_headers    <- gsub("^>", "", ff_headers)
  #ff_headers    <- str_match(ff_headers, "\\|([^\\s]+)\\|")[, 2] 
  
  # phosphosites data
  ph <- fread("temp\\Phosphosites_prepared.tsv", check.names = TRUE)
  
  # full protein sequence fasta
  phprot_acc <- unique(ph$Protein)
  acc.matched <- which(ff_headers %in% phprot_acc)
  length(acc.matched)
  
  start  <- ff_header_pos[acc.matched]
  end    <- ff_header_pos[acc.matched + 1] - 1 # line where the next header starts  - 1
  end[is.na(end)] <- length(fasta.file)        # in the case of the last entry
  
  writeLines(unlist(Map(function(start, end) fasta.file[start:end], start = start, end = end)), "temp\\Phosphoprotein_sequences.fasta")
  
  # prepare sequence window in fasta format
  ph[, Sequence.window_15_noSpace := gsub("_", "", Sequence.window_15)]
  ph_sw7  <- ph[, lapply(.SD, paste, collapse = "/"), by = "Sequence.window_7_noSpace",  .SDcols = c("Protein", "Amino.acid", "Position")]
  ph_sw15 <- ph[, lapply(.SD, paste, collapse = "/"), by = "Sequence.window_15_noSpace", .SDcols = c("Protein", "Amino.acid", "Position")]
  
  sw_fasta <- vector("character", length = nrow(ph_sw7)*2)
  sw_fasta[seq(1, length(sw_fasta), by = 2)] <- paste0(">", ph_sw7$Protein, "___", ph_sw7$Position)
  sw_fasta[seq(2, length(sw_fasta), by = 2)] <- gsub("_", "", ph_sw7$Sequence.window_7_noSpace)
  writeLines(sw_fasta, "temp\\Sequence_window_7.fasta")
  
  sw_fasta <- vector("character", length = nrow(ph_sw15)*2)
  sw_fasta[seq(1, length(sw_fasta), by = 2)] <- paste0(">", ph_sw15$Protein, "___", ph_sw15$Position)
  sw_fasta[seq(2, length(sw_fasta), by = 2)] <- gsub("_", "", ph_sw15$Sequence.window_15_noSpace)
  writeLines(sw_fasta, "temp\\Sequence_window_15.fasta")

})
