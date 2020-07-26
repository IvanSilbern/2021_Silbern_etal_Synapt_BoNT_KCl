
# DO:
# select regulation group for sites
# that were quantified with several multiplicities
# and two or more "multiplicity" versions of the site are significantly regulated.
# Prefer phosphorylation events for which both log2FC are defined, in CaEGTA and MockBTX.
# If two multiplisity state have only one log2FC defined (only possible if Ca-dep. or Cycling-dep.),
# and their regulation groups are different (one is Ca-dep., another is Cycling-dep.), consider the site
# as Ca-compensating
# If multiplicities of the same site belong to different regulation group and have log2FC defined in both experiments,
# use the regulation group of the phosphorylation event with the highest calcium-effect (Ca-dependent or Ca-compensating
# regulation groups will be preferred in this case)

# INPUT:
# "temp\\PhPeptIntensities3.tsv"

# OUTPUT:
# "temp\\PhPeptIntensities4.tsv"


local({

  library(data.table)
    
  ph <- fread("temp\\PhPeptIntensities3.tsv")
  dim(ph)
  
  # Site_id4 does not take multiplicity into account
  # Site_id3 considers different multiplicities
  
  ph[, Site_id4 := paste0(Accession, "-", Position)]
  ph$Site_id4[1:10]
  
  df <- ph[ph$Significance]
  dim(df)
  
  # select sites that are significantly regulated at several multiplicities
  dupl <- df[df$Significance & (df$Site_id4 %in% df$Site_id4[duplicated(df$Site_id4)]), c("Site_id4", "Site_id3", "Multiplicity", "Regulation_group", "log2FC.CaEGTA", "log2FC.BTX")]
  dim(dupl)
  
  # prefer phospho events that have valid log2FC for CaEGTA and BTX experiments
  dupl[, valid := FALSE]
  dupl[!is.na(log2FC.CaEGTA) & !is.na(log2FC.BTX), valid := TRUE]
  
  # take those events that have valid log2FC in both experiments
  dupl_valid <- dupl[dupl$valid]
  
  # subset events that have a valid log2FC in only one experiment
  # do not take it into account if there is the same site with a different multiplicity
  dupl <- dupl[!dupl$valid & !(dupl$Site_id4 %in% dupl_valid$Site_id4)]
  dim(dupl)
  
  # decide on those sites that have one invalid value
  # in this case, only "Ca-dependent" or "Cycling-dependent" entries are possible
  # if both events have same regulation group, keep it
  # if events have different regulation groups, consider them as ca-compensating
  dupl <- dupl[, list(Regulation_group = unlist(lapply(.SD, function(x){
    
    if(length(unique(x)) == 1) return(unique(x))
    else return("Ca-compensating")

  }))), by = "Site_id4", .SDcols = "Regulation_group"]
  
  # add sites with both valid values and no duplication
  dupl <- rbind(dupl, dupl_valid[!dupl_valid$Site_id4 %in% dupl_valid$Site_id4[duplicated(dupl_valid$Site_id4)],
                                 c("Site_id4", "Regulation_group")])
  
  # several multiplicities with both valid values
  dupl_valid <- dupl_valid[dupl_valid$Site_id4 %in% dupl_valid$Site_id4[duplicated(dupl_valid$Site_id4)]]
  dim(dupl_valid)
  
  # Sites where all multiplicity states belong to the same Regulation_group
  temp <- dupl_valid[, list(Regulation_group = unlist(lapply(.SD, function(x){
    
    if(length(unique(x)) == 1) return(unique(x))
    else return("ambiguous")
    
  }))), by = "Site_id4", .SDcols = "Regulation_group"]
  
  # add to previous resolved cases
  dupl <- rbind(dupl, temp[temp$Regulation_group != "ambiguous"])
  
  # keep unresolved cases
  dupl_valid <- dupl_valid[!dupl_valid$Site_id4 %in% temp$Site_id4[temp$Regulation_group != "ambiguous"]]
  
  # select Regulation_group based on stronges Ca-effect (Cycling-dependent will loose against Ca-dependent or Ca-compensating groups)
  dupl_valid[, log2FC.CaEGTA_BTX := log2FC.CaEGTA - log2FC.BTX]
  dupl_valid <- dupl_valid[, list(Regulation_group = Regulation_group[which.max(abs(log2FC.CaEGTA_BTX))]), by = "Site_id4"]
  
  # add to previous resolved cases
  dupl <- rbind(dupl, dupl_valid)
  names(dupl)[names(dupl) == "Regulation_group"] <- "Regulation_group_resolved"
  
  # combine with initial data table
  df <- merge(ph, dupl, by = "Site_id4", all.x = TRUE)
  df[is.na(Regulation_group_resolved), Regulation_group_resolved := Regulation_group]
  df[!df$Significance, Regulation_group_resolved := "not-regulated"]
  
  fwrite(df, "temp\\PhPeptIntensities4.tsv", sep = "\t")
  
})
