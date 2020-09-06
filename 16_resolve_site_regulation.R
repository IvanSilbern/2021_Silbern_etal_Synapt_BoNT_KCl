
# DO:
# select regulation group for sites
# that were quantified with several multiplicities
# and two or more "multiplicity" versions of the site are significantly regulated.
# If one of the multiplicity states is "SV-cycling-dependent", the site is termed "SV-cycling-dependent" as well.

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
  dupl <- df[df$Significance & (df$Site_id4 %in% df$Site_id4[duplicated(df$Site_id4)]), c("Site_id4", "Site_id3", "Multiplicity", "Regulation_group", "log2FC.CaEGTA", "log2FC.BoNT")]
  dim(dupl)
  
  # Only "primary Ca-dependent" or "SV-cycling-dependent" entries are possible
  # if events have same regulation group, keep it
  # if events have different regulation groups, consider them as SV-cycling-dependent
  dupl <- dupl[, list(Regulation_group = unlist(lapply(.SD, function(x){
    
    if(length(unique(x)) == 1) return(unique(x))
    else return("SV-cycling-dependent")
    
  }))), by = "Site_id4", .SDcols = "Regulation_group"]
  
  names(dupl)[names(dupl) == "Regulation_group"] <- "Regulation_group_resolved"
  
  # combine with initial data table
  df <- merge(ph, dupl, by = "Site_id4", all.x = TRUE)
  df[is.na(Regulation_group_resolved), Regulation_group_resolved := Regulation_group]
  df[!df$Significance, Regulation_group_resolved := "not-affected"]
  
  fwrite(df, "temp\\PhPeptIntensities4.tsv", sep = "\t")
  
})
