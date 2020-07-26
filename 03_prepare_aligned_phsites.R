# DO:
# use BLAST alignment of phosphopeptide sequence windows and
# sequence windows in phosphosite plus kinase-substrate data base
# Check that the phosphosite is on the same place within the aligned sequences
# and contains same amino acid
# add information about kinase and substrates from Kinase-Substrate data set

# INPUT:
# "temp\\Phosphosites_prepared.tsv"
# "external\\data_Phosphositeplus\\Kinase_Substrate_Dataset.txt"
# "external\\data_BLASTPhSeq_PhsitePlus_align.txt" (blast alignment results)

# OUTPUT:
# "temp\\Candidate_Phsp_alignment.txt"
# "plots\\Bitscore_density.pdf"

local({

  if(!dir.exists("plots")) dir.create("plots")
  if(!dir.exists("temp")) dir.create("temp")
  
  library(data.table)
  library(stringr)
  library(ggplot2)
  
  # custom functions
  determineSitePosition    <- function(seq_wind, standard.position, wind.length){
    
    # function returns the phosphosite position in the sequence
    # after omitting empty spaces (marked as "_")
    # input:
    # seq_wind = character vector of sequence windows
    # standard.position = int, position of ph site in the sequence window
    # wind.length = int, sequence window length
    
    seq <- gsub("_", "", seq_wind)
    site.position <- integer()
    for(i in seq_along(seq)){
      
      if(nchar(seq[i]) == wind.length) {
        
        site.position[i] <- standard.position
        
      } else {
        
        seq_wind.split   <- unlist(strsplit(seq_wind[i], split = ""))
        empty.positions  <- which(seq_wind.split == "_")
        site.position[i] <- standard.position - sum(empty.positions < standard.position)
      }
      
    }
    
    return(site.position)
    
  }
  alignedSequencePositions <- function(df, seq.col, start.col, end.col){
    # function matches position numbers of aligned part of the sequence
    # gaps are assigned with position 0
    
    sequences <- list()
    seq_wind  <- df[[seq.col]]
    seq_wind  <- str_split(seq_wind, "")
    
    for(i in seq_along(seq_wind)){
      
      start.position <- df[[start.col]][i]
      end.position   <- df[[end.col]][i]
      seq.split <- seq_wind[[i]]
      
      seq.positions <- integer()
      k <- start.position
      for(j in seq_along(seq.split)){
        
        if(seq.split[j] != "-"){ 
          
          seq.positions[j] <- k
          k <- k + 1
          
        } else {
          
          seq.positions[j] <- 0
          
        }
        
      }
      
      if((k-1) != end.position) message(paste0(i, "calculated and reported end positions do not match"))
      
      sequences[[i]] <- seq.positions
      
    }
    
    return(sequences)
    
  }
  
  # load data
  # phosphosite table
  ph <- fread("temp\\Phosphosites_prepared.tsv", check.names = TRUE)
  
  # PhosphositePlus data base
  kin_sub <- fread("external\\data_Phosphositeplus\\Kinase_Substrate_Dataset.txt", skip = 2, check.names = TRUE)
  
  # add kin_sub ids
  kin_sub[, id := 1:.N]
  
  # alignment
  blast.align <- fread("external\\data_BLAST\\PhSeq_PhsitePlus_align.txt", check.names = TRUE)
  names(blast.align) <- c("qacc", "sacc", "evalue", "bitscore", "nident", "pident",
                          "mismatch", "gapopen", "gaps", "qstart", "qend", "sstart",
                          "send", "qseq", "sseq", "qframe", "sframe")
  
  # fix protein name
  blast.align[grepl("db\\|tao2", blast.align$sacc), sacc := paste0("db|tao2_beta|", kin_sub[grepl("tao2", kin_sub$SUB_GENE), id][1])]
  
  # extract sequence ids
  blast.align[, qacc_pos  := lapply(str_split(str_match(qacc, "___(.+)$")[, 2], "/"), as.integer)]
  blast.align[, qacc_prot := str_split(str_match(qacc, "(.+)___.+$")[, 2], "/")]
  
  blast.align[, sacc_id  := lapply(str_split(str_match(sacc, "db\\|[^|]*\\|(.+)$")[, 2], "/"), as.integer)]
  
  # extract 1st id
  blast.align[, qacc_pos1  := unlist(lapply(qacc_pos, "[", 1))]
  blast.align[, qacc_prot1 := unlist(lapply(qacc_prot, "[", 1))]
  blast.align[, sacc_id1   := unlist(lapply(sacc_id, "[", 1))]
  
  blast.align <- merge(blast.align, ph[, c("Protein", "Position", "Sequence.window_15")],
                       by.x = c("qacc_prot1", "qacc_pos1"),
                       by.y = c("Protein", "Position"))
  
  blast.align <- merge(blast.align, kin_sub[, c("id", "SITE_...7_AA")],
                       by.x = c("sacc_id1"),
                       by.y = c("id"))
  
  names(blast.align)[names(blast.align) == "Sequence.window_15"] <- "qacc_sw"
  names(blast.align)[names(blast.align) == "SITE_...7_AA"] <- "sacc_sw"
  
  blast.align[, sacc_sw := toupper(sacc_sw)]
  
  # add query sequence windows (phosphosites table)
  
  # find the ph site position in query and subject sequences
  # standard positions are 16 and 8th positions respectively
  # unless sequence window starts with "_"
  blast.align$ph_position_query <- determineSitePosition(seq_wind = blast.align$qacc_sw,
                                                         standard.position = 16,
                                                         wind.length = 31)
  #table(blast.align$ph_position_query)
  
  blast.align$ph_position_subject <- determineSitePosition(seq_wind = blast.align$sacc_sw,
                                                           standard.position = 8,
                                                           wind.length = 15)
  #table(blast.align$ph_position_subject)
  
  # phosphosites should be present within aligned sequences
  dim(blast.align)
  blast.align <- blast.align[blast.align$qstart <= blast.align$ph_position_query   & blast.align$qend >= blast.align$ph_position_query, ]
  blast.align <- blast.align[blast.align$sstart <= blast.align$ph_position_subject & blast.align$send >= blast.align$ph_position_subject, ]
  dim(blast.align)
  
  
  # find phosphosite position within aligned sequences
  # query sequences
  temp_pos <- alignedSequencePositions(df = blast.align, seq.col = "qseq", start.col = "qstart", end.col = "qend")
  head(temp_pos)
  
  # positions of the aligned aa in the full sequence window
  blast.align[, qpos := unlist(lapply(temp_pos, paste, collapse = ";"))]
  
  # position of the phospho-site in the aligned sequence
  blast.align[, ph_position_qalign := unlist(Map(function(x, y) which(x == y), x = temp_pos, y = blast.align$ph_position_query))]
  
  # extract the central phosphosite based on its position in aligned sequence
  blast.align$qaa <- unlist(Map("substr", x = blast.align$qseq,             
                                start = blast.align$ph_position_qalign,     
                                stop = blast.align$ph_position_qalign))
  #table(blast.align$qaa)
  
  # subject sequences
  temp_pos <- alignedSequencePositions(df = blast.align, seq.col = "sseq", start.col = "sstart", end.col = "send")
  head(temp_pos)
  
  # position of the phospho-site in the aligned sequence
  blast.align[, spos := unlist(lapply(temp_pos, paste, collapse = ";"))]
  
  # positions
  blast.align[, ph_position_salign := unlist(Map(function(x, y) which(x == y), x = temp_pos, y = blast.align$ph_position_subject))]
  
  # extract the central phosphosite based on its position in aligned sequence
  blast.align[, saa := unlist(Map("substr", x = blast.align$sseq,
                                start = blast.align$ph_position_salign,
                                stop = blast.align$ph_position_salign))]
  #table(blast.align$saa)
  
  # check two criteria:
  
  # FIRST: the positions of the phospho-sites in the query and subject sequence should be the same!
  blast.align <- blast.align[blast.align$ph_position_qalign == blast.align$ph_position_salign, ]
  
  # SECOND: 
  # Serins/Threonins should match to Serin/Threonins and Tyrosins to Tyrosins
  # It simplifies to : central tyrosine in query sequences should match to central tyrosines in subject sequences and vice versa
  to.keep <- rep(TRUE, nrow(blast.align))
  
  to.keep[grepl("Y", blast.align$qaa)] <- blast.align$qaa[grepl("Y", blast.align$qaa)] == blast.align$saa[grepl("Y", blast.align$qaa)]
  to.keep[grepl("Y", blast.align$saa)] <- blast.align$qaa[grepl("Y", blast.align$saa)] == blast.align$saa[grepl("Y", blast.align$saa)]
  
  blast.align <- blast.align[to.keep]
  
  # add blast align id
  blast.align[, id:= 1:.N]
  
  #expand query protein sequence & position
  temp <- blast.align[, list(qAccession = unlist(qacc_prot), qPosition = unlist(qacc_pos)), by = "id"]
  blast.align <- merge(temp, blast.align, by = "id", all.x = TRUE)
  
  # expand on subject id
  temp <- blast.align[, list(sId = unlist(sacc_id)), by = "id"]
  blast.align <- merge(temp, blast.align, by = "id", all.x = TRUE)
  
  # add kinase information
  blast.align <- merge(blast.align, kin_sub, by.x = "sId", by.y = "id", all.x = TRUE)
  
  # remove unnecessary columns
  blast.align <- blast.align[, -c("sacc_id1", "qacc_prot1",
                                  "qacc_pos1", "qacc_pos",
                                  "qacc_prot", "sacc_id",
                                  "ph_position_query",
                                  "qpos", "spos")]
  
  fwrite(blast.align, "temp\\Candidate_Phsp_alignment.txt", sep = "\t")
  
  # bitscore density plot
  # highlight nident > 12
  pdf("plots\\Bitscore_density.pdf")
  g <- ggplot(blast.align, aes(x = bitscore, color = blast.align$nident > 12, fill = blast.align$nident > 12))
  g <- g + geom_density(alpha = 0.2)
  g <- g + scale_fill_discrete(label = c("nident < 12", "nident > 12"))
  g <- g + guides(fill = guide_legend(title = "# identical AA"), color = "none")
  g <- g + ggtitle("Sequence window alignment Kinase-Substrate data set: Bitscores")
  g
  dev.off()
  
})
