
# DO:
# remove potential contaminants and reversed sequences
# Add Uniprot annotation scores to proteins
# Extract sequence windows for each protein
# Rank proteins based on number of (unique) identified phosphorylation sites >
# > reviewed protein sequence > protein score > is protein from Swiss Prot db >
# > is protein an isoform > does it have a gene name
# sort protein ids and corresponding sequence positions accordingly
# save modified table under "Phosphosites_prepared.tsv"

# input:
# "search_results\\Synapt_Ph\\combined\\txt\\Phospho (STY)Sites.txt" table loaded as 'ph'
# "external\\data_Uniprot\\UP_AnotationScore.tsv" containing uniprot scores
# "external\\data_Uniprot\\201902_uniprot-rat.fasta" fasta file used for database search

# output:
# "comparison_EK\\temp\\Phosphosites_prepared_IS.tsv"       accordingly modified "Phosphosites (STY).txt" table

local({
  
  if(!dir.exists("temp")) dir.create("temp")
  if(!dir.exists("comparison_EK\\temp")) dir.create("comparison_EK\\temp")
  
  library(data.table)
  library(stringr)
  
  # read in fasta file
  fasta.file    <- fread("external\\data_Uniprot\\201902_uniprot-rat.fasta", sep = NULL, header = FALSE)[[1]]
  
  # extract header positions
  ff_header_pos <- which(grepl("^>", fasta.file)) # Lines containing headers 
  ff_headers    <- str_split_fixed(fasta.file[ff_header_pos], pattern = " ", n = 2)[, 1] # extract protein id without description
  ff_headers    <- gsub("^>", "", ff_headers) # remove '>'
  #ff_headers    <- str_match(ff_headers, "\\|([^\\s]+)\\|")[, 2] # useful if only Uniprot id is reported in Phospho (STY).txt
  
  # load scoring information from UP
  scores <- fread("external\\data_Uniprot\\UP_AnnotationScore.tsv")
  
  # integer scores
  scores[, SCORE2 := as.integer(str_match(scores$SCORE, "^([0-5])\\s")[, 2])]
  
  # re-order the score  table and remove duplicate entries based on UNIPROTKB
  scores <- scores[order(REVIEWED, -SCORE2, UNIPROTKB), ]
  scores <- scores[!duplicated(scores$UNIPROTKB), ]
  
  #import Phospho STY MQ output table
  df <- fread( "search_results\\Synapt_Ph\\combined\\txt\\Phospho (STY)Sites.txt.gz", check.names = TRUE)
  
  # remove phosphosites belonging to contaminants or reverse sequences
  message("Remove ", sum(df$Potential.contaminant == "+"),
          " Potential contaminants and ", sum(df$Reverse == "+"), " Reversed sequences")
  
  df <- df[Potential.contaminant != "+"]
  df <- df[Reverse != "+"]
  
  # Site is unique if no other proteins are reported for this site
  df[, Unique.site := !grepl(";", df$Proteins)]
  
  # one protein - site per row
  df_prot <- df[, list(Protein_single = unlist(strsplit(Proteins, ";")),
                       ProSite_single = unlist(strsplit(as.character(Positions.within.proteins),  ";"))),
                by = id]
  
  # remove potential contaminants and reverse protein sequences if still present
  df_prot <- df_prot[!grepl("CON_", Protein_single)]
  df_prot <- df_prot[!grepl("REV_", Protein_single)]
  
  # find header position in the fasta file
  df_prot[, ff_header_pos :=  ..ff_header_pos[match(df_prot$Protein_single, ..ff_headers)]]
  df_prot[, ff_seq_stop   := (..ff_header_pos[match(df_prot$Protein_single, ..ff_headers) + 1] - 1)]
  df_prot[is.na(ff_seq_stop), ff_seq_stop := length(fasta.file)]
  
  # extract sequence windows
  
  # +/- 7 aa
  seq_wind_size = 7
  df_prot[, Sequence.window_7 := unlist(lapply(seq_along(df_prot$ff_header_pos), function(i){
    
    temp_header_pos <- df_prot$ff_header_pos[i]
    temp_seq_stop   <- df_prot$ff_seq_stop[i]
    temp_seq_pos    <- as.integer(df_prot$ProSite_single[i])
    
    temp_sequence <- paste0(fasta.file[(temp_header_pos + 1) : temp_seq_stop], collapse = "")
    temp_seq_wind <- substr(temp_sequence, start = temp_seq_pos - seq_wind_size, stop = temp_seq_pos + seq_wind_size)
    
    # if the site is located at the end of the sequence: append "___"
    
    missing_start <- temp_seq_pos - (seq_wind_size + 1)
    missing_end   <- nchar(temp_sequence) - (temp_seq_pos + seq_wind_size)
    
    if(missing_start < 0) temp_seq_wind <- paste0(paste0(rep("_", times = abs(missing_start)), collapse = ""), temp_seq_wind)
    if(missing_end < 0)   temp_seq_wind <- paste0(temp_seq_wind, paste0(rep("_", times = abs(missing_end)),   collapse = ""))
    
    temp_seq_wind
    
  }))]
  
  # +/- 15 aa
  seq_wind_size = 15
  df_prot[, Sequence.window_15 := unlist(lapply(seq_along(df_prot$ff_header_pos), function(i){
    
    temp_header_pos <- df_prot$ff_header_pos[i]
    temp_seq_stop   <- df_prot$ff_seq_stop[i]
    temp_seq_pos    <- as.integer(df_prot$ProSite_single[i])
    
    temp_sequence <- paste0(fasta.file[(temp_header_pos + 1) : temp_seq_stop], collapse = "")
    temp_seq_wind <- substr(temp_sequence, start = temp_seq_pos - seq_wind_size, stop = temp_seq_pos + seq_wind_size)
    
    # if the site is located at the end of the sequence: append "___"
    
    missing_start <- temp_seq_pos - (seq_wind_size + 1)
    missing_end   <- nchar(temp_sequence) - (temp_seq_pos + seq_wind_size)
    
    if(missing_start < 0) temp_seq_wind <- paste0(paste0(rep("_", times = abs(missing_start)), collapse = ""), temp_seq_wind)
    if(missing_end < 0)   temp_seq_wind <- paste0(temp_seq_wind, paste0(rep("_", times = abs(missing_end)),   collapse = ""))
    
    temp_seq_wind
    
  }))]
  
  # Protein / Site annotation
  
  # Amino.acid
  df_prot[, Amino.acid := substr(Sequence.window_7, start = 8, stop = 8)]
  
  # sequence windows without space holders
  df_prot[, Sequence.window_7_noSpace := gsub("_", "", Sequence.window_7)]
  
  # Fasta.header
  df_prot[, fasta.header := gsub(">", "", ..fasta.file[ff_header_pos])]
  
  # Gene.name
  df_prot[, Gene.name := str_match(fasta.header, "GN=([A-Za-z0-9]+)")[, 2]]
  
  # Swiss.prot?
  df_prot[, Swiss.prot := grepl("^sp", fasta.header)]
  
  # Accession
  df_prot[, Accession := Protein_single]
  
  # Accession w/o Isoform number
  df_prot[, Accession.noIso := gsub("-[0-9]+$", "", Accession)]
  
  # is protein an isoform?
  df_prot[, Protein.isoform := grepl("-[0-9]+$", Accession)]
  
  # protein description
  df_prot[, Protein.description := str_match(fasta.header, "\\s(.+)\\sOS=")[, 2]]
  
  # add Unique.site column from df
  df_prot <- merge(df_prot, df[, c("id", "Unique.site")], by = "id", all.x = TRUE)
  
  # add UP scoring information
  df_prot <- merge(df_prot, scores[, c("UNIPROTKB", "REVIEWED", "SCORE2")], by.x = "Accession",  by.y = "UNIPROTKB", all.x = TRUE)
  
  # Extract information about proteins in order to rank them
  # one protein per row
  df_prot_single <- df_prot[, list(Nr_Unique     = sum(Unique.site[!duplicated(ProSite_single)]),
                                   Nr_Sites      = length(unique(ProSite_single)),
                                   Sites         = paste(unique(ProSite_single), collapse = ";"),
                                   Ids           = paste(unique(id), collapse = ";"),
                                   Sequence.window_15 = paste(Sequence.window_15, collapse = ";"),
                                   Sequence.window_7  = paste(Sequence.window_7, collapse = ";"),
                                   Sequence.window_7_noSpace = paste(Sequence.window_7_noSpace, collapse = ";"),
                                   Amino.acid         = paste(Amino.acid, collapse = ";")),
                            by = Protein_single]
  
  df_prot_single <- merge(df_prot_single, df_prot[!duplicated(Protein_single), -c("id",
                                                                                  "ProSite_single",
                                                                                  "Unique.site",
                                                                                  "Sequence.window_15",
                                                                                  "Sequence.window_7",
                                                                                  "Sequence.window_7_noSpace",
                                                                                  "Amino.acid")], by = "Protein_single")
  
  # sort proteins = rank
  df_prot_single <- df_prot_single[order(-Nr_Unique, -Nr_Sites, REVIEWED, -SCORE2, -Swiss.prot, Protein.isoform, is.na(Gene.name), Protein_single)]
  
  # assign Protein_id that can be also used to rank proteins
  df_prot_single[, Protein_id := 1:df_prot_single[, .N]]
  
  # add protein rank to df_prot table
  df_prot <- merge(df_prot, df_prot_single[, c("Protein_single", "Protein_id")], by = "Protein_single", all.x = TRUE)
  
  # lists of protein accessions and corresponding modified sequence positions
  prots <- str_split(df$Proteins, ";")
  posit <- str_split(df$Positions.within.proteins, ";")
  
  temp <- data.table()
  pb <- txtProgressBar(min = 1, max = nrow(df), char = "*", style = 3)
  for (i in 1:nrow(df)){
    
    setTxtProgressBar(pb, i)
    
    rank <- unlist(lapply(prots[[i]], function(x) df_prot_single$Protein_id[df_prot_single$Accession == x]))
    new.order <- order(rank)
    rank <- rank[new.order]
    
    temp <- rbind(temp, data.table(Proteins_sorted     = paste(df_prot_single$Protein_single[rank],
                                                               collapse = ";"),
                                   Accessions_sorted   = paste(df_prot_single$Accession[rank],
                                                               collapse = ";"),
                                   Positions_sorted  = paste(posit[[i]][new.order],
                                                             collapse = ";"),
                                   Gene.names_sorted = paste(df_prot_single$Gene.name[rank],
                                                             collapse = ";"),
                                   Fasta.headers_sorted = paste(df_prot_single$fasta.header[rank],
                                                                collapse = ";"),
                                   stringsAsFactors  = FALSE
    ))
    
  }
  
  temp[, Protein      := unlist(lapply(str_split(Proteins_sorted, ";"), "[", 1))]
  temp[, Accession    := unlist(lapply(str_split(Accessions_sorted, ";"), "[", 1))]
  temp[, Accession    := str_match(Accession, "\\|([^\\s]+)\\|")[, 2]]
  temp[, Position     := unlist(lapply(str_split(Positions_sorted, ";"), "[", 1))]
  temp[, Gene.name    := unlist(lapply(str_split(Gene.names_sorted, ";"), "[", 1))]
  temp[, Fasta.header := unlist(lapply(str_split(Fasta.headers_sorted, ";"), "[", 1))]
  temp[, Accession.noIso := gsub("-[0-9]+$", "", Accession)]
  
  df <- cbind(temp, df[, -c("Proteins", "Positions.within.proteins",
                            "Leading.proteins", "Protein",
                            "Fasta.headers", "Sequence.window", "Amino.acid", "Position")])
  df <- merge(df, df_prot[, c("Accession", "id", "Protein.description",
                              "Amino.acid", "Sequence.window_7", "Sequence.window_15",
                              "Sequence.window_7_noSpace", "REVIEWED")], by.x = c("Protein", "id"), by.y = c("Accession", "id"))
  
  fwrite(df, "comparison_EK\\temp\\Phosphosites_prepared_IS.tsv", sep = "\t")
  })