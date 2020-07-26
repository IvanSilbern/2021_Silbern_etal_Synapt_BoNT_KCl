
# DO:

# Convert identified phospho-sites on rat proteins
# to potential phospho-sites on human proteins
# based on BLAST results of rat protein to Human String v9.0 (used by NetworKIN)
# output is used as an input for NetworKIN algorithm

# INPUT:
# "temp\\Phosphosites_prepared.tsv" (phosphosite data table)
# "external\\data_BLAST\\Human_vs_Rat_Networkin.txt.bz2" (bz2-compressed BLAST results)
# ** "external\\data_NetworKIN\\Networkin_unmapped_sites.txt" **(additional input, sites not recognized by NetworKIN)

# OUTPUT:
# "temp\\Human_Stringid_PhSites.tsv" (table of matched Rat and (potential) Human phosphosites)

local({
  
  if(!dir.exists("temp")) dir.create("temp")
    
  library(data.table)
  library(stringr)
  
  # custom functions
  findIndex <- function(site, seq, sstart, aa){
    
    # Provide site position before alignment
    # and sequence after alignment (with gaps)
    ngaps    <- 0
    seq_len  <- nchar(seq)
    site_new <- 0
    site     <- site - sstart + 1 # account for start of the alignment
    
    while(site + ngaps != site_new && site_new <= seq_len){
      
      site_new <- site + ngaps
      aa_new   <- str_sub(seq, start = site_new, end = site_new)
      
      while(aa_new == "-" && site_new <= seq_len){
        
        site_new <- site_new + 1
        aa_new   <- str_sub(seq, start = site_new, end = site_new)
        
      }
      
      ngaps <- str_count(str_sub(seq, 1, site_new), "-")
      
    }
    
    if(aa == aa_new) site_new
    else NA
    
  }
  findPosition <- function(site_index, qseq, qstart, saa){
    
    qaa <- str_sub(qseq, start = site_index, end = site_index)
    if(!saa %in% c("S", "T", "Y") || !qaa %in% c("S", "T", "Y")) return(NA)
    if(saa %in% c("S", "T") && !qaa %in% c("S", "T")) return(NA)
    if(saa == "Y" && ! qaa != "Y") return (NA)
    
    ngaps <- str_count(str_sub(qseq, start = 1, end = site_index), "-")
    
    site_pos <- site_index - ngaps
    site_pos <- site_pos + qstart - 1
    
    return(site_pos)
    
  }  
  findSite <- function(site_index, qseq, qstart, saa){
    
    qaa <- str_sub(qseq, start = site_index, end = site_index)
    if(!saa %in% c("S", "T", "Y") || !qaa %in% c("S", "T", "Y")) return(NA)
    if(saa %in% c("S", "T") && !qaa %in% c("S", "T")) return(NA)
    if(saa == "Y" && ! qaa != "Y") return (NA)
    
    return(qaa)
    
  }
  
  ph <- fread("temp\\Phosphosites_prepared.tsv", check.names = TRUE)
  
  blast <- fread("external\\data_BLAST\\Human_vs_Rat_Networkin.txt.bz2", header = FALSE)
  names(blast) <- c("qacc", "sacc", "evalue", "bitscore", "nident", "pident",
                    "mismatch", "gapopen", "gaps", "qstart", "qend", "sstart",
                    "send", "qlen", "slen", "length", "qseq", "sseq", "qframe", "sframe")
  
  # calculate % coverage of subject sequence (rat protein)
  blast[, sseq_coverage := 100*(send - sstart + 1) / slen]
  
  # order & keep best matches
  blast  <- blast[order(-sseq_coverage, -pident, -bitscore, gapopen)]
  blast2 <- blast[!duplicated(blast$sacc)]
  
  # keep only phosphoproteins
  blast2 <- blast2[blast2$sacc %in% unique(ph$Protein)]
  
  # combine protien / site with blast results
  temp <- unique(ph[, c("Protein", "Position" ,"Amino.acid")])
  temp <- merge(temp, blast2, by.x = "Protein", by.y = "sacc", all.x = TRUE)
  dim(temp)
  
  # modified site should lie within the alignment
  temp <- temp[temp$Position >= temp$sstart & temp$Position <= temp$send]
  dim(temp)
  
  # add site position in the aligned sequence
  temp[, site_index := unlist(Map(findIndex, site = Position, seq = sseq, sstart = sstart, aa = Amino.acid))]
  
  # add site and its position in the query sequence (HUMAN protein sequence)
  temp[, qseq_pos  := unlist(Map(findPosition, site_index = site_index, qseq = qseq, qstart = qstart, saa = Amino.acid))]
  temp <- temp[!is.na(qseq_pos)]
  
  temp[, qseq_site := unlist(Map(findSite,     site_index = site_index, qseq = qseq, qstart = qstart, saa = Amino.acid))]
  
  # subset columns
  temp2 <- temp[, c("Protein", "Position", "Amino.acid", "qacc", "qseq_pos", "qseq_site")]
  names(temp2) <- c("Rat_Protein", "Rat_Position", "Rat_Residue", "Human_Protein", "Human_Position", "Human_Residue")
  
  temp2[, Human_Protein := str_replace(Human_Protein, "^9606\\.", "")]
  temp2 <- temp2[!duplicated(temp2)]
  
  
  fwrite(temp2, "temp\\Human_Stringid_PhSites.tsv", sep = "\t")
  
  
  # for some reason NetworKIN cannot recognise some sites
  # not recognized sites are removed from the data set
  
  unmapped <- fread("external\\data_NetworKIN\\Networkin_unmapped_sites.txt")
  names(unmapped) <- c("Human_Protein", "Human_Position", "Human_Residue")
  dim(temp2)
  
  temp2[, Siteid := paste0(Human_Protein, "_", Human_Position, "_", Human_Residue)]
  unmapped[, Siteid := paste0(Human_Protein, "_", Human_Position, "_", Human_Residue)]
  
  dim(temp2)
  temp2 <- temp2[!temp2$Siteid %in% unmapped$Siteid]
  dim(temp2)
  
  fwrite(temp2[, -c("Siteid")], "temp\\Human_Stringid_PhSites.tsv", sep = "\t")

})
