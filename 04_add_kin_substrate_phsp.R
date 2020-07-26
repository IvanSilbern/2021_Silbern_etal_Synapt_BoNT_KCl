# DO:
# Add kinase-substrate information to phosphosites table
# based on BLAST results (sequence windows in phosphosites data table
# vs sequence windows in phosphosite plus kinase-substrate data table)

# INPUT:
# "temp\\Candidate_Phsp_alignment.txt"
# "temp\\Phosphosites_prepared.tsv" (phosphosites table)

# OUTPUT:
# updated "temp\\Phosphosites_prepared.tsv"

local({
  
  if(!dir.exists("temp")) dir.create("temp")
  
  library(data.table)

  # phosphosites data
  ph <- fread("temp\\Phosphosites_prepared.tsv")
  
  # load alignment table of candidate siequence windows vs phosphosite plus data base
  cand_phsp_align <- fread("temp\\Candidate_Phsp_alignment.txt", check.names = TRUE)
  
  # use only substrates and kinases from rat, mouse, human, and rabbit
  allow_species <- c("rat", "mouse", "human", "rabbit")
  cand_phsp_align <- cand_phsp_align[cand_phsp_align$KIN_ORGANISM %in% allow_species & 
                                     cand_phsp_align$SUB_ORGANISM %in% allow_species]
  cand_phsp_align[, KIN_ORGANISM := factor(KIN_ORGANISM, levels = allow_species)]
  cand_phsp_align[, SUB_ORGANISM := factor(SUB_ORGANISM, levels = allow_species)]
  
  # set a bitscore cutoff and order matches based on bitscore, nident and gapopen, KIN_ORGANISM, and SUB_ORGANISM
  dim(cand_phsp_align)
  cand_phsp_align <- cand_phsp_align[order(-bitscore, -nident, gapopen, KIN_ORGANISM, SUB_ORGANISM)]
  
  # remove duplicated matches
  cand_phsp_align <- cand_phsp_align[!duplicated(cand_phsp_align[, c("qAccession", "qPosition")])]
  dim(cand_phsp_align)
  
  # bitscore cutoff
  cand_phsp_align <- cand_phsp_align[cand_phsp_align$bitscore > 20]
  
  # merge with blast results
  take.columns <- c("qAccession", "qPosition", "GENE",
                    "KIN_ACC_ID", "KIN_ORGANISM", "SUB_GENE",
                    "SUB_ACC_ID", "SUB_ORGANISM", "SUB_MOD_RSD", "DOMAIN")
  
  ph <- merge(ph, cand_phsp_align[, ..take.columns], by.x = c("Protein", "Position"),
              by.y = c("qAccession", "qPosition"), all.x = TRUE, all.y = FALSE)
  
  fwrite(ph, "temp\\Phosphosites_prepared.tsv", sep = "\t")

})