# Do: Site occupancy calculation following Hoegrebe et al

# INPUT:
# "temp\\PhPeptIntensities4.tsv"
# "search_results\\Unbound_MockBoNT_AD\\combined\\txt\\peptides.txt.gz"
# "search_results\\Unbound_MockBoNT_AD\\combined\\txt\\proteinGroups.txt.gz"
# "search_results\\Unbound_MockBoNT_AD\\combined\\txt\\proteinGroups.txt.gz"
# "search_results\\Unbound_MockBoNT_CB\\combined\\txt\\proteinGroups.txt.gz"
# "search_results\\Synapt_Ph\\combined\\txt\\modificationSpecificPeptides.txt.gz"

# OUTPUT:
# "plots\\Occupancy_phpept.pdf"
# "Figures\\SupplFig_4\\site_occupancies.txt"
# "temp\\PhPeptIntensities5.tsv"

local({
  
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  
  if(!dir.exists("temp")) dir.create("temp")
  if(!dir.exists("Figures\\SupplFig_4")) dir.create("Figures\\SupplFig_4", recursive = TRUE)
  
  samples  <- c("Mock.01", "Mock.02", "Mock.03", "BoNT.01", "BoNT.02", "BoNT.03")
  
  sites <- fread("temp\\PhPeptIntensities4.tsv")
  pept1 <- fread("search_results\\Unbound_MockBoNT_AD\\combined\\txt\\peptides.txt.gz", check.names = TRUE)
  pept2 <- fread("search_results\\Unbound_MockBoNT_CB\\combined\\txt\\peptides.txt.gz", check.names = TRUE)
  pg1 <- fread("search_results\\Unbound_MockBoNT_AD\\combined\\txt\\proteinGroups.txt.gz", check.names = TRUE)
  pg1 <- fread("search_results\\Unbound_MockBoNT_CB\\combined\\txt\\proteinGroups.txt.gz", check.names = TRUE)
  ph  <- fread("search_results\\Synapt_Ph\\combined\\txt\\modificationSpecificPeptides.txt.gz", check.names = TRUE)
  
  # subset phosphorylated peptides
  ph <- ph[grepl("Phospho", Modifications)]
  
  # remove those that are NA in BTX_AD_01 or BTX_CB_02
  ph <- ph[!is.na(Experiment.BTX_AD_01) | !is.na(Experiment.BTX_CB_02)]
  
  # subset usefull columns
  take_cols <- c("id", "Sequence", "Modifications", "Protein.group.IDs", "Phospho..STY..site.IDs",
                 "Acetyl..Protein.N.term.", "Oxidation..M.", "Phospho..STY.")
  int_cols  <- grep("Reporter\\.intensity\\.corrected\\.\\d\\.BTX", names(ph), value = TRUE)
  ph <- ph[, .SD, .SDcols = c(take_cols, int_cols)]
  names(ph) <- gsub("Reporter\\.intensity\\.corrected\\.", "X", names(ph))
  
  # split into two data tables based on experiment
  int_cols1 <- names(ph)[grepl("^X", names(ph)) & grepl("AD_01", names(ph))]
  ph1 <- ph[, .SD, .SDcols = c(take_cols, int_cols1)]
  
  int_cols2 <- names(ph)[grepl("^X", names(ph)) & grepl("CB_02", names(ph))]
  ph2 <- ph[, .SD, .SDcols = c(take_cols, int_cols2)]
  
  # remove experiment names from column names
  names(ph1) <- gsub("\\.BTX_.+", "", names(ph1))
  names(ph2) <- gsub("\\.BTX_.+", "", names(ph2))
  
  # use condition in the column names
  names(ph1)[names(ph1) == "X1"] <- "Mock.01"
  names(ph1)[names(ph1) == "X2"] <- "BoNT.01"
  names(ph1)[names(ph1) == "X3"] <- "Mock.02"
  names(ph1)[names(ph1) == "X4"] <- "BoNT.02"
  names(ph1)[names(ph1) == "X5"] <- "BoNT.03"
  names(ph1)[names(ph1) == "X6"] <- "Mock.03"
  
  # use condition in the column names
  names(ph2)[names(ph2) == "X1"] <- "Mock.01"
  names(ph2)[names(ph2) == "X2"] <- "BoNT.01"
  names(ph2)[names(ph2) == "X3"] <- "Mock.02"
  names(ph2)[names(ph2) == "X4"] <- "BoNT.02"
  names(ph2)[names(ph2) == "X5"] <- "Mock.03"
  names(ph2)[names(ph2) == "X6"] <- "BoNT.03"
  
  # calculate number of 0 values in total and per condition
  int_mock <- grep("^Mock", names(ph1), value = TRUE)
  int_bont <- grep("^BoNT", names(ph1), value = TRUE)
  ph1[, n.zero := sum(unlist(.SD) == 0), by = "id", .SDcols = c(int_mock, int_bont)]
  ph1[, n.zero.mock := sum(unlist(.SD) == 0), by = "id", .SDcols = int_mock]
  ph1[, n.zero.bont := sum(unlist(.SD) == 0), by = "id", .SDcols = int_bont]
  ph2[, n.zero := sum(unlist(.SD) == 0), by = "id", .SDcols = c(int_mock, int_bont)]
  ph2[, n.zero.mock := sum(unlist(.SD) == 0), by = "id", .SDcols = int_mock]
  ph2[, n.zero.bont := sum(unlist(.SD) == 0), by = "id", .SDcols = int_bont]
  
  # delete if
  # n.zero.mock > 1 & n.zero.bont > 1  
  # n.zero.mock == 3 | n.zero.bont == 3
  
  ph1[, to.delete := FALSE]
  ph1[n.zero.mock > 1 & n.zero.bont > 1, to.delete := TRUE]
  ph1[n.zero.mock == 3 | n.zero.bont == 3, to.delete := TRUE]
  
  ph2[, to.delete := FALSE]
  ph2[n.zero.mock > 1 & n.zero.bont > 1, to.delete := TRUE]
  ph2[n.zero.mock == 3 | n.zero.bont == 3, to.delete := TRUE]
  
  ph1 <- ph1[!to.delete == TRUE]
  ph2 <- ph2[!to.delete == TRUE]
  
  #check how many 0 left
  sum(ph1$n.zero > 0)
  sum(ph2$n.zero > 0)
  
  # make 0 to NA
  ph1[n.zero > 0, c(int_mock, int_bont) := lapply(.SD, function(x) { 
    
    x[x == 0] <- NA_real_
    x
    
  }
  ), .SDcols = c(int_mock, int_bont)]
  ph2[n.zero > 0, c(int_mock, int_bont) := lapply(.SD, function(x) { 
    
    x[x == 0] <- NA_real_
    x
    
  }
  ), .SDcols = c(int_mock, int_bont)]
  
  ph1[n.zero > 0]
  ph2[n.zero > 0]
  
  # use median polishing for normalization
  dat1 <- ph1[, lapply(.SD, log2), .SDcols = samples]
  dat2 <- ph2[, lapply(.SD, log2), .SDcols = samples]
  
  dat1 <- medpolish(dat1, na.rm = TRUE)
  dat2 <- medpolish(dat2, na.rm = TRUE)
  
  dat1 <- as.data.table(dat1$residuals)
  dat2 <- as.data.table(dat2$residuals)
  
  dat1 <- dat1[, lapply(.SD, function(x) 2^x), .SDcols = samples]
  dat2 <- dat2[, lapply(.SD, function(x) 2^x), .SDcols = samples]
  
  ph1 <- cbind(ph1[, -c(..samples)], dat1)
  ph2 <- cbind(ph2[, -c(..samples)], dat2)
  
  # calculate mean and sd
  ph1[, mock_mean := mean(as.numeric(.SD), na.rm = TRUE), by = "id", .SDcols = int_mock]
  ph1[, bont_mean := mean(as.numeric(.SD), na.rm = TRUE), by = "id", .SDcols = int_bont]
  ph1[, mock_sd := sd(as.numeric(.SD), na.rm = TRUE), by = "id", .SDcols = int_mock]
  ph1[, bont_sd := sd(as.numeric(.SD), na.rm = TRUE), by = "id", .SDcols = int_bont]
  
  ph2[, mock_mean := mean(as.numeric(.SD), na.rm = TRUE), by = "id", .SDcols = int_mock]
  ph2[, bont_mean := mean(as.numeric(.SD), na.rm = TRUE), by = "id", .SDcols = int_bont]
  ph2[, mock_sd := sd(as.numeric(.SD), na.rm = TRUE), by = "id", .SDcols = int_mock]
  ph2[, bont_sd := sd(as.numeric(.SD), na.rm = TRUE), by = "id", .SDcols = int_bont]
  
  # compute X which is Phosphopeptide intensity ratio Mock vs BoNT
  ph1[, X := mock_mean / bont_mean]
  ph1[, X_sd := sqrt((mock_sd/bont_mean)^2 + (mock_mean*bont_sd/bont_mean^2)^2)]
  
  ph2[, X := mock_mean / bont_mean]
  ph2[, X_sd := sqrt((mock_sd/bont_mean)^2 + (mock_mean*bont_sd/bont_mean^2)^2)]
  
  # add peptide id indicating peptide sequence and number of modification
  ph1[, pept_id := paste0(Sequence, "_A", Acetyl..Protein.N.term., "_O", Oxidation..M.)]
  ph2[, pept_id := paste0(Sequence, "_A", Acetyl..Protein.N.term., "_O", Oxidation..M.)]
  
  ph1[, phpept_id := paste0(Sequence, "_A", Acetyl..Protein.N.term., "_O", Oxidation..M., "_P", Phospho..STY.)]
  ph2[, phpept_id := paste0(Sequence, "_A", Acetyl..Protein.N.term., "_O", Oxidation..M., "_P", Phospho..STY.)]
  
  # prepare not-phophosphorylated counterparts
  nph1 <- pept1[Sequence %in% ph1$Sequence]
  nph2 <- pept2[Sequence %in% ph2$Sequence]
  
  take_cols <- c("id", "Sequence", "Protein.group.IDs", "Missed.cleavages")
  int_cols  <- grep("Reporter\\.intensity\\.corrected\\.\\d$", names(nph1), value = TRUE)
  
  # subset useful columns
  nph1 <- nph1[, .SD, .SDcols = c(take_cols, int_cols)]
  nph2 <- nph2[, .SD, .SDcols = c(take_cols, int_cols)]
  
  # rename intensity columns
  names(nph1) <- gsub("Reporter\\.intensity\\.corrected\\.", "X", names(nph1))
  names(nph2) <- gsub("Reporter\\.intensity\\.corrected\\.", "X", names(nph2))
  
  # use condition in the column names
  names(nph1)[names(nph1) == "X1"] <- "Mock.01"
  names(nph1)[names(nph1) == "X2"] <- "BoNT.01"
  names(nph1)[names(nph1) == "X3"] <- "Mock.02"
  names(nph1)[names(nph1) == "X4"] <- "BoNT.02"
  names(nph1)[names(nph1) == "X5"] <- "Mock.03"
  names(nph1)[names(nph1) == "X6"] <- "BoNT.03"
  
  # use condition in the column names
  names(nph2)[names(nph2) == "X1"] <- "Mock.01"
  names(nph2)[names(nph2) == "X2"] <- "BoNT.01"
  names(nph2)[names(nph2) == "X3"] <- "Mock.02"
  names(nph2)[names(nph2) == "X4"] <- "BoNT.02"
  names(nph2)[names(nph2) == "X5"] <- "Mock.03"
  names(nph2)[names(nph2) == "X6"] <- "BoNT.03"
  
  # check number of 0
  nph1[, n.zero := sum(as.numeric(.SD) == 0), by = "id", .SDcols = c(int_mock, int_bont)]
  nph1[, n.zero.mock := sum(as.numeric(.SD) == 0), by = "id", .SDcols = int_mock]
  nph1[, n.zero.bont := sum(as.numeric(.SD) == 0), by = "id", .SDcols = int_bont]
  nph2[, n.zero := sum(as.numeric(.SD) == 0), by = "id", .SDcols = c(int_mock, int_bont)]
  nph2[, n.zero.mock := sum(as.numeric(.SD) == 0), by = "id", .SDcols = int_mock]
  nph2[, n.zero.bont := sum(as.numeric(.SD) == 0), by = "id", .SDcols = int_bont]
  
  # delete if
  # n.zero.mock > 1 & n.zero.bont > 1  
  # n.zero.mock == 3 | n.zero.bont == 3
  
  nph1[, to.delete := FALSE]
  nph1[n.zero.mock > 1 & n.zero.bont > 1, to.delete := TRUE]
  nph1[n.zero.mock == 3 | n.zero.bont == 3, to.delete := TRUE]
  
  nph2[, to.delete := FALSE]
  nph2[n.zero.mock > 1 & n.zero.bont > 1, to.delete := TRUE]
  nph2[n.zero.mock == 3 | n.zero.bont == 3, to.delete := TRUE]
  
  nph1 <- nph1[!to.delete == TRUE]
  nph2 <- nph2[!to.delete == TRUE]
  
  #check how many 0 left
  sum(nph1$n.zero > 0)
  sum(nph2$n.zero > 0)
  
  # make 0 to NA
  nph1[n.zero > 0, c(int_mock, int_bont) := lapply(.SD, function(x) { 
    
    x[x == 0] <- NA_real_
    x
    
  }
  ), .SDcols = c(int_mock, int_bont)]
  nph2[n.zero > 0, c(int_mock, int_bont) := lapply(.SD, function(x) { 
    
    x[x == 0] <- NA_real_
    x
    
  }
  ), .SDcols = c(int_mock, int_bont)]
  
  nph1[n.zero > 0]
  nph2[n.zero > 0]
  
  # use median polishing for normalization
  dat1 <- nph1[, lapply(.SD, log2), .SDcols = samples]
  dat2 <- nph2[, lapply(.SD, log2), .SDcols = samples]
  
  dat1 <- medpolish(dat1, na.rm = TRUE)
  dat2 <- medpolish(dat2, na.rm = TRUE)
  
  dat1 <- as.data.table(dat1$residuals)
  dat2 <- as.data.table(dat2$residuals)
  
  dat1 <- dat1[, lapply(.SD, function(x) 2^x), .SDcols = samples]
  dat2 <- dat2[, lapply(.SD, function(x) 2^x), .SDcols = samples]
  
  nph1 <- cbind(nph1[, -c(..samples)], dat1)
  nph2 <- cbind(nph2[, -c(..samples)], dat2)
  
  # calculate mean and sd
  nph1[, mock_mean := mean(unlist(.SD), na.rm = TRUE), by = "id", .SDcols = int_mock]
  nph1[, bont_mean := mean(unlist(.SD), na.rm = TRUE), by = "id", .SDcols = int_bont]
  nph1[, mock_sd := sd(unlist(.SD), na.rm = TRUE), by = "id", .SDcols = int_mock]
  nph1[, bont_sd := sd(unlist(.SD), na.rm = TRUE), by = "id", .SDcols = int_bont]
  
  nph2[, mock_mean := mean(unlist(.SD), na.rm = TRUE), by = "id", .SDcols = int_mock]
  nph2[, bont_mean := mean(unlist(.SD), na.rm = TRUE), by = "id", .SDcols = int_bont]
  nph2[, mock_sd := sd(unlist(.SD), na.rm = TRUE), by = "id", .SDcols = int_mock]
  nph2[, bont_sd := sd(unlist(.SD), na.rm = TRUE), by = "id", .SDcols = int_bont]
  
  # compute Y which is non-phosnphopeptide intensity ratio Mock vs BoNT
  nph1[, Y := mock_mean / bont_mean]
  nph1[, Y_sd := sqrt((mock_sd/bont_mean)^2 + (mock_mean*bont_sd/bont_mean^2)^2)]
  
  nph2[, Y := mock_mean / bont_mean]
  nph2[, Y_sd := sqrt((mock_sd/bont_mean)^2 + (mock_mean*bont_sd/bont_mean^2)^2)]
  
  # use only the first protein group id
  nph1[, Protein.group.ID := unlist(lapply(str_split(Protein.group.IDs, ";"), "[", 1))]
  nph2[, Protein.group.ID := unlist(lapply(str_split(Protein.group.IDs, ";"), "[", 1))]
  
  # prepare protein intensities based on peptides that are not involved into phosphorylation
  prot1 <- pept1[!Sequence %in% nph1$Sequence]
  prot2 <- pept2[!Sequence %in% nph2$Sequence]
  
  # subset useful columns
  prot1 <- prot1[, .SD, .SDcols = c(take_cols, int_cols)]
  prot2 <- prot2[, .SD, .SDcols = c(take_cols, int_cols)]
  
  # rename intensity columns
  names(prot1) <- gsub("Reporter\\.intensity\\.corrected\\.", "X", names(prot1))
  names(prot2) <- gsub("Reporter\\.intensity\\.corrected\\.", "X", names(prot2))
  
  # use condition in the column names
  names(prot1)[names(prot1) == "X1"] <- "Mock.01"
  names(prot1)[names(prot1) == "X2"] <- "BoNT.01"
  names(prot1)[names(prot1) == "X3"] <- "Mock.02"
  names(prot1)[names(prot1) == "X4"] <- "BoNT.02"
  names(prot1)[names(prot1) == "X5"] <- "Mock.03"
  names(prot1)[names(prot1) == "X6"] <- "BoNT.03"
  
  # use condition in the column names
  names(prot2)[names(prot2) == "X1"] <- "Mock.01"
  names(prot2)[names(prot2) == "X2"] <- "BoNT.01"
  names(prot2)[names(prot2) == "X3"] <- "Mock.02"
  names(prot2)[names(prot2) == "X4"] <- "BoNT.02"
  names(prot2)[names(prot2) == "X5"] <- "Mock.03"
  names(prot2)[names(prot2) == "X6"] <- "BoNT.03"
  
  # check number of 0
  prot1[, n.zero := sum(unlist(.SD) == 0), by = "id", .SDcols = c(int_mock, int_bont)]
  prot1[, n.zero.mock := sum(unlist(.SD) == 0), by = "id", .SDcols = int_mock]
  prot1[, n.zero.bont := sum(unlist(.SD) == 0), by = "id", .SDcols = int_bont]
  prot2[, n.zero := sum(unlist(.SD) == 0), by = "id", .SDcols = c(int_mock, int_bont)]
  prot2[, n.zero.mock := sum(unlist(.SD) == 0), by = "id", .SDcols = int_mock]
  prot2[, n.zero.bont := sum(unlist(.SD) == 0), by = "id", .SDcols = int_bont]
  
  # delete if
  # n.zero.mock > 1 & n.zero.bont > 1  
  # n.zero.mock == 3 | n.zero.bont == 3
  
  prot1[, to.delete := FALSE]
  prot1[n.zero.mock > 1 & n.zero.bont > 1, to.delete := TRUE]
  prot1[n.zero.mock == 3 | n.zero.bont == 3, to.delete := TRUE]
  
  prot2[, to.delete := FALSE]
  prot2[n.zero.mock > 1 & n.zero.bont > 1, to.delete := TRUE]
  prot2[n.zero.mock == 3 | n.zero.bont == 3, to.delete := TRUE]
  
  prot1 <- prot1[!to.delete == TRUE]
  prot2 <- prot2[!to.delete == TRUE]
  
  # expand protein group id column
  prot1 <- merge(prot1[, list(Protein.group.ID = unlist(str_split(Protein.group.IDs, ";"))), by = "id"],
                 prot1, by = "id", all.x = TRUE)
  
  prot2 <- merge(prot2[, list(Protein.group.ID = unlist(str_split(Protein.group.IDs, ";"))), by = "id"],
                 prot2, by = "id", all.x = TRUE)
  
  # subset protein groups for which there are ph-nph pairs
  prot1 <- prot1[Protein.group.ID %in% nph1$Protein.group.ID]
  prot2 <- prot2[Protein.group.ID %in% nph2$Protein.group.ID]
  
  # check number of peptides per protein group
  prot1[, N.pept := .N, by = "Protein.group.ID"]
  prot2[, N.pept := .N, by = "Protein.group.ID"]
  
  # check how many protein groups would have only 1 peptide for quantification
  sum(prot1[, N.pept[1], by = "Protein.group.ID"]$V1 == 1)
  sum(prot2[, N.pept[1], by = "Protein.group.ID"]$V1 == 1)
  
  # combine peptide intensities
  prot1 <- prot1[, lapply(.SD, sum), by = Protein.group.ID, .SDcols = c(int_mock, int_bont)]
  prot2 <- prot2[, lapply(.SD, sum), by = Protein.group.ID, .SDcols = c(int_mock, int_bont)]
  
  # use median polishing for normalization
  dat1 <- prot1[, lapply(.SD, log2), .SDcols = samples]
  dat2 <- prot2[, lapply(.SD, log2), .SDcols = samples]
  
  dat1 <- dat1[, lapply(.SD, function(x) {
    
    x[abs(x) == Inf] <- NA
    return(x)
    
  }), .SDcols = samples]
  
  dat2 <- dat2[, lapply(.SD, function(x){
    
    x[abs(x) == Inf] <- NA
    return(x)
    
  }), .SDcols = samples]
  
  dat1 <- medpolish(dat1, na.rm = TRUE)
  dat2 <- medpolish(dat2, na.rm = TRUE)
  
  dat1 <- as.data.table(dat1$residuals)
  dat2 <- as.data.table(dat2$residuals)
  
  dat1 <- dat1[, lapply(.SD, function(x) 2^x), .SDcols = samples]
  dat2 <- dat2[, lapply(.SD, function(x) 2^x), .SDcols = samples]
  
  prot1 <- cbind(prot1[, -c(..samples)], dat1)
  prot2 <- cbind(prot2[, -c(..samples)], dat2)
  
  # calculate mean and sd
  prot1[, mock_mean := mean(unlist(.SD), na.rm = TRUE), by = "Protein.group.ID", .SDcols = int_mock]
  prot1[, bont_mean := mean(unlist(.SD), na.rm = TRUE), by = "Protein.group.ID", .SDcols = int_bont]
  prot1[, mock_sd := sd(unlist(.SD), na.rm = TRUE), by = "Protein.group.ID", .SDcols = int_mock]
  prot1[, bont_sd := sd(unlist(.SD), na.rm = TRUE), by = "Protein.group.ID", .SDcols = int_bont]
  
  prot2[, mock_mean := mean(unlist(.SD), na.rm = TRUE), by = "Protein.group.ID", .SDcols = int_mock]
  prot2[, bont_mean := mean(unlist(.SD), na.rm = TRUE), by = "Protein.group.ID", .SDcols = int_bont]
  prot2[, mock_sd := sd(unlist(.SD), na.rm = TRUE), by = "Protein.group.ID", .SDcols = int_mock]
  prot2[, bont_sd := sd(unlist(.SD), na.rm = TRUE), by = "Protein.group.ID", .SDcols = int_bont]
  
  # compute Z which is protein intensity ratio Mock vs BoNT
  prot1[, Z := mock_mean / bont_mean]
  prot1[, Z_sd := sqrt((mock_sd/bont_mean)^2 + (mock_mean*bont_sd/bont_mean^2)^2)]
  
  prot2[, Z := mock_mean / bont_mean]
  prot2[, Z_sd := sqrt((mock_sd/bont_mean)^2 + (mock_mean*bont_sd/bont_mean^2)^2)]
  
  # combine Z values with Y values
  
  nph1 <- merge(nph1, prot1[, c("Protein.group.ID", "Mock.01", "Mock.02", "Mock.03", "BoNT.01", "BoNT.02", "BoNT.03", "Z", "Z_sd")], by = "Protein.group.ID", all.x = TRUE, suffixes = c("", ".prot"))
  nph2 <- merge(nph2, prot2[, c("Protein.group.ID", "Mock.01", "Mock.02", "Mock.03", "BoNT.01", "BoNT.02", "BoNT.03", "Z", "Z_sd")], by = "Protein.group.ID", all.x = TRUE, suffixes = c("", ".prot"))
  
  sum(is.na(nph1$Z))
  sum(is.na(nph2$Z))
  
  # assume protein intensity of 1 for cases where Z could not be calculated
  nph1[is.na(Z), paste0("Mock.0", 1:3, ".prot") := 1]
  nph1[is.na(Z), paste0("BoNT.0", 1:3, ".prot") := 1]
  nph2[is.na(Z), paste0("Mock.0", 1:3, ".prot") := 1]
  nph2[is.na(Z), paste0("BoNT.0", 1:3, ".prot") := 1]
  
  # combine X values with Y and Z
  int_ph   <- paste0(samples, ".ph")
  int_nph  <- paste0(samples, ".pept")
  int_prot <- paste0(samples, ".prot")
  
  ph1 <- merge(ph1, nph1[, c("Sequence", ..samples, "Y", "Y_sd", "Z", "Z_sd", ..int_prot)], by = "Sequence", suffixes = c(".ph", ".pept"))
  ph2 <- merge(ph2, nph2[, c("Sequence", ..samples, "Y", "Y_sd", "Z", "Z_sd", ..int_prot)], by = "Sequence", suffixes = c(".ph", ".pept"))
  
  # merge ph1 and ph2
  
  ph <- merge(ph1[, .SD, .SDcols = c("id", "Phospho..STY..site.IDs", int_ph, int_nph, int_prot)],
              ph2[, .SD, .SDcols = c("id", int_ph, int_nph, int_prot)], by = "id", suffixes = c(".1", ".2"), all = TRUE
  )
  
  # 3DMM
  
  int_ph <- grep("\\.ph", names(ph), value = TRUE)
  int_nph <- grep("\\.pept", names(ph), value = TRUE)
  int_prot <- grep("\\.prot", names(ph), value = TRUE)
  
  models    <- apply(ph[, .SD, .SDcols = c(int_ph, int_nph, int_prot)], 1, function(x) tryCatch(lm(x[1:12] ~ x[25:36] + x[13:24] + c(rep(0, 6), rep(1, 6))), error = function(e) NA))
  coef      <- unlist(lapply(models, function(x) tryCatch(coefficients(summary(x))[3, 1], error = function(e) NA)))
  p.val     <- unlist(lapply(models, function(x) tryCatch(coefficients(summary(x))[3, 4], error = function(e) NA)))
  
  hist(p.val, breaks = 100)
  
  ph[, m := ..coef]
  ph[, occ.p.val := ..p.val]
  
  occ <-  mapply(function(ph, nph, m) -ph/(nph*m)/(1-ph/(nph*m)), ph = ph[, ..int_ph], nph =ph[, ..int_nph], m = ph[, "m"])
  occ <- as.data.table(occ)
  names(occ) <- c(paste0("Occupancy.", samples, ".1"), paste0("Occupancy.", samples, ".2"))
  
  ph <- cbind(ph, occ)
  
  # mean occupancy per condition
  ph[, occupancy_mock := mean(unlist(.SD), na.rm = TRUE), by = "id", .SDcols = c(paste0("Occupancy.Mock.0", 1:3, ".1"), paste0("Occupancy.Mock.0", 1:3, ".2"))]
  ph[, occupancy_bont := mean(unlist(.SD), na.rm = TRUE), by = "id", .SDcols = c(paste0("Occupancy.BoNT.0", 1:3, ".1"), paste0("Occupancy.Mock.0", 1:3, ".2"))]
  
  # keep only plausible values
  ph <- ph[(occupancy_mock > 0 & occupancy_mock < 1) & (occupancy_bont > 0 & occupancy_bont < 1)]
  
  ph[, occupancy_mock_bont := occupancy_mock - occupancy_bont]
  
  # combine results from 2 experiments
  take_cols <- c("id", "Phospho..STY..site.IDs", "occupancy_mock", "occupancy_bont", "occupancy_mock_bont", "occ.p.val")
  res <- ph[, .SD, .SDcols = take_cols]
  
  res <- res[!is.na(occupancy_bont) & !is.na(occupancy_mock)]
  
  # expand res table base on ph sites id
  res <- merge(res[, list(Phospho..STY..site.IDs = unlist(str_split(Phospho..STY..site.IDs, ";"))), by = "id"],
               res[, -c("Phospho..STY..site.IDs")], by = "id", all.x = TRUE)
  
  res[, Phospho..STY..site.IDs := as.integer(Phospho..STY..site.IDs)]
  
  temp <- merge(sites[, c("Site_id3", "id", "log2FC.BoNT")], res[, c("Phospho..STY..site.IDs", "occupancy_bont", "occupancy_mock", "occ.p.val")],
                by.x = "id", by.y = "Phospho..STY..site.IDs", all.x = TRUE)
  temp <- temp[!is.na(occupancy_bont)]
  temp <- temp[!is.na(log2FC.BoNT)]
  temp <- temp[, list(occupancy_bont = mean(occupancy_bont),
                      occupancy_mock = mean(occupancy_mock),
                      occ.p.val      = min(occ.p.val)),
               by = "Site_id3"]
  
  sites <- merge(sites, temp, by = "Site_id3", all.x = TRUE)
  
  sites[, occupancy_mock_bont := occupancy_mock - occupancy_bont]
  sites[, label := ""]
  sites[, Site_id5 := paste0(Gene.name, "-", Position, "(", Multiplicity, ")")]
  sites[which(abs(log2FC.BoNT) > 0.5 | abs(log2FC.CaEGTA) > 0.5), label := Site_id5]
  sites[occupancy_mock < 0.4, label := ""]
  sites[Gene.name %in% c("Vamp2", "Cnr1", "Stx1a", "Ctnnd2","Rock2", "Syn1", "Dnm1",
                         "Ywhaz", "Atp1a2", "Atp1a3"), label := Site_id2]
  sites[, label := Site_id5]
  sites[Gene.name %in% c("LOC100362814"), label := ""]

  g <- ggplot(sites[!is.na(occupancy_mock) & Regulation_group_resolved != "not-affected" & occ.p.val < 0.1],
              aes(x = log2FC.BoNT, y = log2FC.CaEGTA, label = label,
                  size = occupancy_mock, color = Regulation_group_resolved))
  g <- g + coord_equal(x = c(-2.2, 3.2), y = c(-2.2, 3.2), expand = FALSE)
  g <- g + scale_color_manual(values = c("#fcb533", "#51baf4"),
                              labels = c(expression(paste("primary Ca"^"2+", "-dependent")),
                                         "SV-cycling-dependent"))
  g <- g + geom_hline(yintercept = log2(c(1.2, 1/1.2)), color = "darkred", linetype = "dashed")
  g <- g + geom_vline(xintercept = log2(c(1.2, 1/1.2)), color = "darkred", linetype = "dashed")
  g <- g + geom_point(alpha = 0.6)
  g <- g + geom_text_repel(size = 3.5,
                           color = scales::alpha("black", 0.6),
                           segment.color = scales::alpha("darkgrey", 0.8))
  g <- g + guides(color = guide_legend(title = "Regulation Group", override.aes = list(alpha = 0.6, size = 5), label.hjust = 0),
                  size = guide_legend(title = "Occupancy in Mock"))
  g <- g + xlab(expression(paste("log"[2], " (Mock / BoNT)")))
  g <- g + ylab(expression(paste("log"[2], " (Ca / EGTA)")))
  g <- g + theme_bw()

  pdf("plots\\Occupancy_phpept.pdf", width = 11)
  print(g)
  dev.off()
  
  # export figure data
  fwrite(sites[!is.na(occupancy_mock) & Regulation_group_resolved != "not-affected" & occ.p.val < 0.1,
               c("Gene.name",
                 "Position",
                 "Multiplicity",
                 "Site_id5",
                 "label",
                 "log2FC.CaEGTA",
                 "log2FC.BoNT",
                 "occupancy_mock",
                 "occupancy_bont",
                 "occ.p.val",
                 "Regulation_group_resolved")], 
         "Figures\\SupplFig_4\\site_occupancies.txt", sep = "\t")
  
  # median occupancy
  median(res[occ.p.val < 0.1, occupancy_mock[1], by = id]$V1)
  
  # save occupancy data
  fwrite(sites, "temp\\PhPeptIntensities5.tsv", sep = "\t")
  
})
