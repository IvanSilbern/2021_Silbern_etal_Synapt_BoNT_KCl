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
  
  # compute CaEGTA - BTX
  df[,      log2FC.CaEGTA_BTX := log2FC.CaEGTA - log2FC.BTX]
  
  #Add Site id
  df[, Site_id := paste0(Gene.name, "-", Position, Amino.acid)]
  df[, Site_id2 := paste0(Gene.name, "-", Position, Amino.acid, "_", Multiplicity)]
  df[, Site_id3 := paste0(Accession.noIso, "-", Position, Amino.acid, "_", Multiplicity)]
  
  sum(duplicated(df$Site_id2))
  sum(duplicated(df$Site_id3))
  
  # distribute sites in different regulation groups
  fc_cut <- log2(1.2)
  # not regulated
  df[(log2FC.CaEGTA_BTX <  fc_cut & log2FC.CaEGTA_BTX > -fc_cut) & (log2FC.BTX <  fc_cut & log2FC.BTX > -fc_cut), Regulation_group := "not-regulated"]
  # cycling-dependent
  df[(log2FC.CaEGTA_BTX <  fc_cut & log2FC.CaEGTA_BTX > -fc_cut) & (log2FC.BTX >  fc_cut), Regulation_group := "Cycling-dependent"]
  df[(log2FC.CaEGTA_BTX <  fc_cut & log2FC.CaEGTA_BTX > -fc_cut) & (log2FC.BTX < -fc_cut), Regulation_group := "Cycling-dependent"]
  # Ca-compensating
  df[(log2FC.CaEGTA_BTX < -fc_cut) & (log2FC.BTX >  fc_cut), Regulation_group := "Ca-compensating"]
  df[(log2FC.CaEGTA_BTX >  fc_cut) & (log2FC.BTX < -fc_cut), Regulation_group := "Ca-compensating"]
  # Ca-enhancing
  df[(log2FC.CaEGTA_BTX >  fc_cut) & (log2FC.BTX >  fc_cut), Regulation_group := "Ca-enhancing"]
  df[(log2FC.CaEGTA_BTX < -fc_cut) & (log2FC.BTX < -fc_cut), Regulation_group := "Ca-enhancing"]
  # Ca-dependent
  df[(log2FC.CaEGTA_BTX >  fc_cut) & (log2FC.BTX <  fc_cut & log2FC.BTX > -fc_cut), Regulation_group := "Ca-dependent"]
  df[(log2FC.CaEGTA_BTX < -fc_cut) & (log2FC.BTX <  fc_cut & log2FC.BTX > -fc_cut), Regulation_group := "Ca-dependent"]
  
  # direction
  df[log2FC.CaEGTA_BTX >  fc_cut, Regulation_direction := "pos"]
  df[log2FC.CaEGTA_BTX < -fc_cut, Regulation_direction := "neg"]
  df[Regulation_group == "Cycling-dependent" & log2FC.BTX >  fc_cut, Regulation_direction := "pos"]
  df[Regulation_group == "Cycling-dependent" & log2FC.BTX < -fc_cut, Regulation_direction := "neg"]
  
  # Ca-dependent if regulated in CaEGTA but not found in BTX experiment
  df[which(abs(log2FC.CaEGTA) > fc_cut & is.na(Candidate.BTX)), Regulation_group := "Ca-dependent"]
  df[log2FC.CaEGTA >  fc_cut & is.na(Candidate.BTX), Regulation_direction := "pos"]
  df[log2FC.CaEGTA < -fc_cut & is.na(Candidate.BTX), Regulation_direction := "neg"]
  
  # Cyclin-dependent if regulated in BTX experiment but not found in CaEGTA experiment
  df[which(abs(log2FC.BTX) > fc_cut & is.na(Candidate.CaEGTA)), Regulation_group := "Cycling-dependent"]
  df[log2FC.BTX >  fc_cut & is.na(Candidate.CaEGTA), Regulation_direction := "pos"]
  df[log2FC.BTX < -fc_cut & is.na(Candidate.CaEGTA), Regulation_direction := "neg"]
  
  
  # Cycling-dependent group should show q.val.BTX < 0.01
  # Ca-dependent group should show q.val.CaEGTA < 0.01
  # Ca-compensating or Ca-enhancing group should show q.val < 0.01 in either CaEGTA or BTX
  
  dim(df)
  df[, Significance := FALSE]
  
  df[Regulation_group == "Cycling-dependent" & q.val.BTX < 0.01, Significance := TRUE]
  df[Regulation_group == "Ca-dependent" & q.val.CaEGTA < 0.01, Significance := TRUE]
  df[Regulation_group == "Ca-compensating" & (q.val.CaEGTA < 0.01 | q.val.BTX < 0.01), Significance := TRUE]
  df[Regulation_group == "Ca-enhancing" & (q.val.CaEGTA < 0.01 | q.val.BTX < 0.01), Significance := TRUE]
  
  df[is.na(df$Regulation_group), Regulation_group := "not-regulated"]
  
  sum(df$Significance)
  length(unique(df$Gene.name[df$Significance]))
  
  fwrite(df, "temp\\PhPeptIntensities2.tsv", sep = "\t")

})
