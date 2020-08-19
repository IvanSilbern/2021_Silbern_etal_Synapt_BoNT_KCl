# DO:
#
# restructure output tables
#
# INPUT:
# "temp\\PhPeptIntensities.tsv"
# "temp\\PhPeptIntensities_slim.tsv"
# "temp\\PhIntensities_cand.tsv"
#
# OUTPUT:
# "temp\\PhPeptIntensities2.tsv"

local({
  
  library(data.table)
    
  # Re-structure output
  
  df <- fread("temp\\PhPeptIntensities.tsv")
  
  names(df)
  dim(df)
  
  # compute CaEGTA - BoNT
  df[,      log2FC.CaEGTA_BoNT := log2FC.CaEGTA - log2FC.BoNT]
  
  #Add Site id
  df[, Site_id := paste0(Gene.name, "-", Position, Amino.acid)]
  df[, Site_id2 := paste0(Gene.name, "-", Position, Amino.acid, "_", Multiplicity)]
  df[, Site_id3 := paste0(Accession.noIso, "-", Position, Amino.acid, "_", Multiplicity)]
  
  sum(duplicated(df$Site_id2))
  sum(duplicated(df$Site_id3))
  
  # distribute sites in different regulation groups
  fc_cut <- log2(1.2)
  
  # not regulated
  df[, Regulation_group := "not-affected"]
  
  # cycling-dependent
  df[abs(log2FC.BoNT) > fc_cut, Regulation_group := "SV-cycling-dependent"]

  # Ca-dependent
  df[abs(log2FC.CaEGTA) >  fc_cut & abs(log2FC.BoNT) <  fc_cut, Regulation_group := "primary Ca-dependent"]
  
  # primary Ca-dependent if regulated in CaEGTA but not found in BoNT experiment
  df[which(abs(log2FC.CaEGTA) > fc_cut & is.na(Candidate.BoNT)), Regulation_group := "primary Ca-dependent"]
  
  # SV-cycling-dependent if regulated in BoNT experiment but not found in CaEGTA experiment
  df[which(abs(log2FC.BoNT) > fc_cut & is.na(Candidate.CaEGTA)), Regulation_group := "SV-cycling-dependent"]
  
  # primary Ca-dependent group should show q.val.CaEGTA < 0.01
  # SV-cycling-dependent group should show q.val.BoNT < 0.01 or be a Candidate in CaEGTA (abs(log2FC.CaEGTA) > log2(1.2) & q.val.CaEGTA < 0.01)
  
  df[, Significance := FALSE]
  
  df[Regulation_group == "primary Ca-dependent" & q.val.CaEGTA < 0.01, Significance := TRUE]
  df[Regulation_group == "SV-cycling-dependent" & (q.val.BoNT < 0.01 | Candidate.CaEGTA), Significance := TRUE]
  
  df[is.na(df$Regulation_group), Regulation_group := "not-affected"]
  
  sum(df$Significance)
  length(unique(df$Gene.name[df$Significance]))
  
  fwrite(df, "temp\\PhPeptIntensities2.tsv", sep = "\t")

})
