# Do:
# prepare Uniprot domain and regions data

# INPUT:
# "external\\data_Uniprot\\Annotation_Domain_Uniprot.tsv"

# OUTPUT:
# temp\\Domains.tsv

local({
  
  library(data.table)
  library(stringr)
  
  df <- fread("external\\data_Uniprot\\Annotation_Domain_Uniprot.tsv", check.names = TRUE)
  dom <- df[, c("Entry", "Region", "Domain..FT.")]
  names(dom)[names(dom) == "Entry"] <- "Accession"
  dom
  
##### regions #####
  coord <- str_match_all(dom$Region, "REGION <?(\\d+)(\\.\\.(\\d+))?")
  coord <- lapply(coord, function(x) x[, -3, drop = FALSE])
  
  notes <- str_match_all(dom$Region, "/note=\"(.+?)\"")
  
  # combine matrices, convert to data.tables
  regions <- Map(function(x, y){
    
  colnames(x) <- c("Domain_Type", "Start", "End")
  colnames(y) <- c("Note", "Description")
  
  as.data.table(cbind(x, y))
    
  }, x = coord, y = notes)
  
  # add protein accession
  regions <- Map(function(x, y){
    
    x[, Accession := y]
    
  }, x = regions, y = as.list(dom$Accession))   
  
  regions <- regions[unlist(lapply(regions, nrow)) != 0]
  regions <- rbindlist(regions)
  
##### domain ft #####
  coord <- str_match_all(dom$Domain..FT., "DOMAIN <?(\\d+)(\\.\\.(\\d+))?")
  coord <- lapply(coord, function(x) x[, -3, drop = FALSE])
  
  notes <- str_match_all(dom$Domain..FT., "/note=\"(.+?)\"")
  
  # combine matrices, convert to data.tables
  domains <- Map(function(x, y){
    
    colnames(x) <- c("Domain_Type", "Start", "End")
    colnames(y) <- c("Note", "Description")
    
    as.data.table(cbind(x, y))
    
  }, x = coord, y = notes)
  
  # add protein accession
  domains <- Map(function(x, y){
    
    x[, Accession := y]
    
  }, x = domains, y = as.list(dom$Accession))   
  
  domains <- domains[unlist(lapply(domains, nrow)) != 0]
  domains <- rbindlist(domains)
  
  #which(unlist(lapply(coord, nrow)) != unlist(lapply(notes, nrow)))
  
##### finalize #####
  
  # combine
  domains <- rbind(domains, regions)
  domains[, Domain_Type := str_match(Domain_Type, "^\\w+")]
  dim(domains)
  
  # Break long descriptions
  domains[, Description_wrapped := paste0(Domain_Type, ": ", Description)]
  domains[, Description_wrapped := str_wrap(domains$Description_wrapped, 40)]
  
  fwrite(domains, "temp\\Domains.tsv", sep = "\t")

    
})

