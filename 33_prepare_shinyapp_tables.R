
# DO:
# prepare tables for Shiny App

# INPUT:
# "temp\\PhPeptIntensities4.tsv"

# OUTPUT:
# "temp\\PhPeptInt_long.tsv"
# "temp\\PhPeptInt_CaEGTA.tsv"
# "temp\\PhPeptInt_BoNT.tsv"
# "temp\\PhPeptInt_CaEGTA_BoNT.tsv"

# prepare data
local({
  
  library(data.table)
  
  dat <- fread("temp\\PhPeptIntensities4.tsv")
  
  # mark not significant sites as not regulated
  dat[!dat$Significance, Regulation_group := "not-affected"]
  
  dat <- dat[order(Gene.name, Position)]
  dat <- dat[, Site_id := paste0(Accession, "-", Amino.acid, Position)]
  
  # prepare data frame with norm intensities in long format for CaEGTA experiment
  df <- dat[, c("Gene.name", "Accession", "Position", "Amino.acid", "Site_id", "Multiplicity",
                grep("^CaEGTA.+\\.norm$", names(..dat), value = TRUE),
                "log2FC.CaEGTA", "q.val.CaEGTA", "Regulation_group_resolved")]
  df <- melt(df, measure.vars = grep("\\.norm", names(df), value = TRUE), variable.name = "Experiment", value.name = "Norm.intensity", variable.factor = FALSE)
  names(df)[names(df) == "log2FC.CaEGTA"] <- "log2FC"
  names(df)[names(df) == "q.val.CaEGTA"]  <- "q.val"
  df <- df[!is.na(df$log2FC)]
  
  # select one multiplicity state based on the magnitude of log2FC
  # consider significantly regulated sites first
  temp <- df[df$Regulation_group_resolved != "not-affected",
             list(Amino.acid       = Amino.acid[which(abs(log2FC) == max(abs(log2FC)))],
                  Multiplicity     = Multiplicity[which(abs(log2FC) == max(abs(log2FC)))],
                  log2FC           = log2FC[which(abs(log2FC) == max(abs(log2FC)))],
                  q.val            = q.val[which(abs(log2FC) == max(abs(log2FC)))],
                  Experiment       = Experiment[which(abs(log2FC) == max(abs(log2FC)))],
                  Norm.intensity   = Norm.intensity[which(abs(log2FC) == max(abs(log2FC)))]
                 ), by = c("Gene.name", "Accession", "Position", "Site_id")]
  
  # same for not-affected sites
  df <- rbind(temp, df[df$Regulation_group_resolved == "not-affected" & !Site_id %in% temp$Site_id,
                       list(Amino.acid       = Amino.acid[which(abs(log2FC) == max(abs(log2FC)))],
                            Multiplicity     = Multiplicity[which(abs(log2FC) == max(abs(log2FC)))],
                            log2FC           = log2FC[which(abs(log2FC) == max(abs(log2FC)))],
                            q.val            = q.val[which(abs(log2FC) == max(abs(log2FC)))],
                            Experiment       = Experiment[which(abs(log2FC) == max(abs(log2FC)))],
                            Norm.intensity   = Norm.intensity[which(abs(log2FC) == max(abs(log2FC)))]
                       ), by = c("Gene.name", "Accession", "Position", "Site_id")])
  
  # long format data
  df_long <- df
  
  # prepare data frame with norm intensities in long format for BoNT experiment
  df <- dat[, c("Gene.name", "Accession", "Position", "Amino.acid", "Site_id", "Multiplicity",
                grep("^MockBoNT.+\\.norm$", names(..dat), value = TRUE),
                "log2FC.BoNT", "q.val.BoNT", "Regulation_group_resolved")]
  df <- melt(df, measure.vars = grep("\\.norm", names(df), value = TRUE), variable.name = "Experiment", value.name = "Norm.intensity", variable.factor = FALSE)
  names(df)[names(df) == "log2FC.BoNT"] <- "log2FC"
  names(df)[names(df) == "q.val.BoNT"]  <- "q.val"

  # select one multiplicity state based on the magnitude of log2FC
  # consider significantly regulated sites first
  temp <- df[df$Regulation_group_resolved != "not-affected",
             list(Amino.acid       = Amino.acid[which(abs(log2FC) == max(abs(log2FC)))],
                  Multiplicity     = Multiplicity[which(abs(log2FC) == max(abs(log2FC)))],
                  log2FC           = log2FC[which(abs(log2FC) == max(abs(log2FC)))],
                  q.val            = q.val[which(abs(log2FC) == max(abs(log2FC)))],
                  Experiment       = Experiment[which(abs(log2FC) == max(abs(log2FC)))],
                  Norm.intensity   = Norm.intensity[which(abs(log2FC) == max(abs(log2FC)))]
                  ), by = c("Gene.name", "Accession", "Position", "Site_id")]
  
  # same for not-affected sites
  df <- rbind(temp, df[df$Regulation_group_resolved == "not-affected" & !Site_id %in% temp$Site_id,
                       list(Amino.acid       = Amino.acid[which(abs(log2FC) == max(abs(log2FC)))],
                            Multiplicity     = Multiplicity[which(abs(log2FC) == max(abs(log2FC)))],
                            log2FC           = log2FC[which(abs(log2FC) == max(abs(log2FC)))],
                            q.val            = q.val[which(abs(log2FC) == max(abs(log2FC)))],
                            Experiment       = Experiment[which(abs(log2FC) == max(abs(log2FC)))],
                            Norm.intensity   = Norm.intensity[which(abs(log2FC) == max(abs(log2FC)))]
                       ), by = c("Gene.name", "Accession", "Position", "Site_id")])
  
  
  df_long<- rbind(df_long, df)
  
  # ExperimentID
  df_long[grepl("^CaEGTA_01", df_long$Experiment), ExperimentID := "CaEGTA_01"]
  df_long[grepl("^CaEGTA_02", df_long$Experiment), ExperimentID := "CaEGTA_02"]
  df_long[grepl("^MockBoNT_03", df_long$Experiment), ExperimentID := "BoNT_01"]
  df_long[grepl("^MockBoNT_04", df_long$Experiment), ExperimentID := "BoNT_02"]
  
  # Condition
  df_long[grepl("_Ca_", df_long$Experiment), Condition := "Ca"]
  df_long[grepl("_EGTA_", df_long$Experiment), Condition := "EGTA"]
  df_long[grepl("_Mock_", df_long$Experiment), Condition := "Mock"]
  df_long[grepl("_BoNT_", df_long$Experiment), Condition := "BoNT"]
  
  fwrite(df_long, "temp\\PhPeptInt_long.tsv")
  
  annot <- dat[, lapply(.SD, "[", 1), by = "Accession",
               .SDcols = c("Protein", "Protein.description", "Gene.name", "Length", "Function")]
  annot[, Function := gsub("\\s\\{.+?\\}", "", Function)]
  annot[, Function := gsub("\\s\\(PubMed.+?\\)", "", Function)]
  annot[, Function := gsub("\\.\\.", ".", Function)]

##### log2FC caegta #####

  df <- dat[!is.na(dat$log2FC.CaEGTA)]
  names(df)[names(df) == "log2FC.CaEGTA"] <- "log2FC"
  
  # select multiplicity with the highest magnitude
  temp <- df[df$Regulation_group_resolved != "not-affected",
             list(Amino.acid       = Amino.acid[which.max(abs(log2FC))],
                  Multiplicity     = Multiplicity[which.max(abs(log2FC))],
                  log2FC           = log2FC[which.max(abs(log2FC))],
                  Regulation_group = Regulation_group_resolved[which.max(abs(log2FC))]
                  ), by = c("Accession", "Position", "Site_id")]
  
  df <- rbind(temp, df[df$Regulation_group_resolved == "not-affected" & !Site_id %in% temp$Site_id,
                       list(Amino.acid       = Amino.acid[which.max(abs(log2FC))],
                            Multiplicity     = Multiplicity[which.max(abs(log2FC))],
                            log2FC           = log2FC[which.max(abs(log2FC))],
                            Regulation_group = Regulation_group_resolved[which.max(abs(log2FC))]
                       ),
                       by = c("Accession", "Position", "Site_id")])
  
  # create Siteid
  df[, Siteid := paste0(Amino.acid, Position)]
  
  # mark enriched sites
  df[, Enriched := FALSE]
  df[df$log2FC > 0, Enriched := TRUE]
  
  # merge with protein annotation
  df <- merge(df, annot, by = "Accession", all.x = TRUE)
  fwrite(df, "temp\\PhPeptInt_CaEGTA.tsv")
  
  

###### log2FC BoNT #####
  
  df <- dat[!is.na(dat$log2FC.BoNT)]
  names(df)[names(df) == "log2FC.BoNT"] <- "log2FC"
  
  # select multiplicity with the highest magnitude
  temp <- df[df$Regulation_group_resolved != "not-affected",
             list(Amino.acid       = Amino.acid[which.max(abs(log2FC))],
                  Multiplicity     = Multiplicity[which.max(abs(log2FC))],
                  log2FC           = log2FC[which.max(abs(log2FC))],
                  Regulation_group = Regulation_group_resolved[which.max(abs(log2FC))]
                  ), by = c("Accession", "Position", "Site_id")]
  
  df <- rbind(temp, df[df$Regulation_group_resolved == "not-affected" & !Site_id %in% temp$Site_id,
                       list(Amino.acid       = Amino.acid[which.max(abs(log2FC))],
                            Multiplicity     = Multiplicity[which.max(abs(log2FC))],
                            log2FC           = log2FC[which.max(abs(log2FC))],
                            Regulation_group = Regulation_group_resolved[which.max(abs(log2FC))]
                       ),
                       by = c("Accession", "Position", "Site_id")])
  
  # create Siteid
  df[, Siteid := paste0(Amino.acid, Position)]
  
  # mark enriched sites
  df[, Enriched := FALSE]
  df[df$log2FC > 0, Enriched := TRUE]
  
  # merge with protein annotation
  df <- merge(df, annot, by = "Accession", all.x = TRUE)
  fwrite(df, "temp\\PhPeptInt_BoNT.tsv")
  
  

##### log2FC caegta_BoNT #####
  
  df <- dat[!is.na(dat$log2FC.CaEGTA_BoNT)]
  names(df)[names(df) == "log2FC.CaEGTA_BoNT"] <- "log2FC"
  
  # select multiplicity with the highest magnitude
  temp <- df[df$Regulation_group_resolved != "not-affected",
             list(Amino.acid       = Amino.acid[which.max(abs(log2FC))],
                  Multiplicity     = Multiplicity[which.max(abs(log2FC))],
                  log2FC           = log2FC[which.max(abs(log2FC))],
                  Regulation_group = Regulation_group_resolved[which.max(abs(log2FC))]
            ), by = c("Accession", "Position", "Site_id")]
  
  df <- rbind(temp, df[df$Regulation_group == "not-affected" & !Site_id %in% temp$Site_id,
                       list(Amino.acid       = Amino.acid[which.max(abs(log2FC))],
                            Multiplicity     = Multiplicity[which.max(abs(log2FC))],
                            log2FC           = log2FC[which.max(abs(log2FC))],
                            Regulation_group = Regulation_group_resolved[which.max(abs(log2FC))]
                       ), by = c("Accession", "Position", "Site_id")])
  
  # create Siteid
  df[, Siteid := paste0(Amino.acid, Position)]
  
  # mark enriched sites
  df[, Enriched := FALSE]
  df[df$log2FC > 0, Enriched := TRUE]
  
  # merge with protein annotation
  df <- merge(df, annot, by = "Accession", all.x = TRUE)
  fwrite(df, "temp\\PhPeptInt_CaEGTA_BoNT.tsv")
  
})