
# Do:
#
# Use prepared phosphorylation site table
# from the study by Engholm-Keller et al
# to compute log2FC relatively to the 0 time point
# Intensities are normalized/imputed in the same way
# as done for our data set
# at least 2 out of 3 time points should be quantified
# for the site to be quantified in the replicate

# INPUT:
# "comparison_EK\\temp\\Phosphosites_prepared_EK.tsv"

# OUTPUT:
# comparison_EK\\plots\\Reporter_Int_distributions.pdf
# comparison_EK\\temp\\PhSiteContrasts_EK.tsv

local({

  if(!dir.exists("comparison_EK\\temp")) dir.create("comparison_EK\\temp")
  if(!dir.exists("comparison_EK\\plots")) dir.create("comparison_EK\\plots")
  
  library(data.table)
  
  # minimum localization probability
  prob.cutoff <- 0.75
    
  ph <- fread("comparison_EK\\temp\\Phosphosites_prepared_EK.tsv", check.names = TRUE, quote = "")
  
  # subset phosphosites with a given minimal localizatin probability
  message("Sites meeting the localization probability cutoff: ", sum(ph$Localization.prob >= prob.cutoff))
  
  ph  <- ph[Localization.prob >= prob.cutoff]
  dim(ph)

  # subset reporter intensity data and convert to a long format #####
  
  # column numbers corresponding to Reporter intensities
  take  <- grep("Intensity\\.[LMH]\\.Rep[0-9_]+", names(ph), value = TRUE)
  
  # additional columns to be taken
  take <- c("id", take)

  # subset corresponding columns from the data frame
  temp <- ph[, ..take]
  temp <- cbind(temp[, 1], temp[, lapply(.SD, as.numeric), .SDcols = take[-1]])

  # convert df to a long format
  temp <- melt(temp, measure.vars = grep("Intensity\\.[LMH]\\.Rep[0-9_]+", names(temp), value = TRUE), variable.factor = FALSE)

  # extract experiment name from 'variable'
  temp[, Experiment := str_match(variable, "Intensity\\.[LMH]\\.(Rep[0-9_]+)")[, 2]]

  # extract Channel from 'variable'
  temp[, Channel := str_match(variable, "Intensity\\.([LHM])\\.")[, 2]]

  # Extract Set number from 'variable'
  temp[, Set := str_match(Experiment, "[28]0_([12])")[, 2]]

  # Extract Stimulation parameter
  temp[, Stimulation := str_match(Experiment, "_([28]0)_")[, 2]]

  # Extract Replicate number
  temp[, Replicate := str_match(Experiment, "Rep([0-9])_")[, 2]]

  # Timepoint
  
  # Rep 1, 2
  temp[Replicate %in% c("1", "2") &
       Channel == "L", Timepoint := 0]
  temp[Replicate %in% c("1", "2") &
         Set == "1" &
         Channel == "M", Timepoint := 10]
  temp[Replicate %in% c("1", "2") &
         Set == "1" &
         Channel == "H", Timepoint := 90]
  temp[Replicate %in% c("1", "2") &
         Set == "2" &
         Channel == "M", Timepoint := 300]
  temp[Replicate %in% c("1", "2") &
         Set == "2" &
         Channel == "H", Timepoint := 900]
  
  # Rep 3, 4
  temp[Replicate %in% c("3", "4") &
         Channel == "M", Timepoint := 0]
  temp[Replicate %in% c("3", "4") &
         Set == "1" &
         Channel == "H", Timepoint := 10]
  temp[Replicate %in% c("3", "4") &
         Set == "1" &
         Channel == "L", Timepoint := 90]
  temp[Replicate %in% c("3", "4") &
         Set == "2" &
         Channel == "H", Timepoint := 300]
  temp[Replicate %in% c("3", "4") &
         Set == "2" &
         Channel == "L", Timepoint := 900]
  
  # Rep 5, 6
  temp[Replicate %in% c("5", "6") &
         Channel == "H", Timepoint := 0]
  temp[Replicate %in% c("5", "6") &
         Set == "1" &
         Channel == "L", Timepoint := 10]
  temp[Replicate %in% c("5", "6") &
         Set == "1" &
         Channel == "M", Timepoint := 90]
  temp[Replicate %in% c("5", "6") &
         Set == "2" &
         Channel == "L", Timepoint := 300]
  temp[Replicate %in% c("5", "6") &
         Set == "2" &
         Channel == "M", Timepoint := 900]
  
  ### remove entries having too many zero values
  n.zero <- temp[, list(n.zero = sum(value == 0)), by = c("id", "Experiment")]
  n.zero[, to.delete := n.zero > 1]
  min.value <- temp[value != 0, list(min.value = min(value, na.rm = TRUE)),by = c("Experiment", "Channel")]
  
  norm.min <- temp[value != 0, list(SD = sd(log2(value), na.rm = TRUE),
                                    percentile = quantile(log2(value), probs = 0.05, na.rm = T)), by = c("Experiment", "Channel")]
  
  temp <- merge(temp, n.zero,    all.x = TRUE)
  temp <- merge(temp, min.value, by = c("Experiment", "Channel"), all.x = TRUE)
  temp <- merge(temp, norm.min, by = c("Experiment", "Channel"), all.x = TRUE)
  
  temp[, value2 := value]
  temp[to.delete == TRUE, value2 := NA_real_]
  
  # % of missing (zero) values per channel
  temp[!is.na(value2), 100*sum(value2 == 0)/.N, by = c("Experiment", "Channel")] 
  
  # use probability distribution for missing value imputation
  set.seed(100)
  temp[value2 == 0, value2 := 2^abs(rnorm(n = .N, mean = percentile, sd = 2*SD))]
  
  # temp[value2 == 0, value2 := min.value]
  temp <- temp[!is.na(value2)]
  
  temp <- dcast(temp, Experiment + n.zero + id ~ Timepoint, value.var = "value2")
  dim(temp)
  
  # extract set, replicate number and stimulation
  temp[, Set := str_match(Experiment, "_([0-9])$")[, 2]]
  temp[, Stimulation := str_match(Experiment, "_([28]0)_")[, 2]]
  temp[, Replicate   := str_match(Experiment, "Rep([0-9])_")[, 2]]
  
  # switch to timepoint 0, 1, 2 instead of exact timing
  temp[Set == "2", `10` := `300`]
  temp[Set == "2", `90` := `900`]
  temp <- temp[, -c("300", "900")]
  
  # rename to timepoint T0, T1, T2
  names(temp)[names(temp) == "0"]  <- "T0"
  names(temp)[names(temp) == "10"] <- "T1"
  names(temp)[names(temp) == "90"] <- "T2"
  
  # log2
  temp[, c("T0", "T1", "T2") := lapply(.SD, log2), .SDcols = c("T0", "T1", "T2")]
  
  # split df based on Experiment
  temp <- split(temp, f = temp$Experiment)
  
  pdf("comparison_EK\\plots\\Reporter_Int_distributions.pdf")
  line_colors  <- c("black", "blue", "green", "cyan", "red")
  line_types   <- rep(1, 5)
  take.columns <- c("T0", "T1", "T2")
  
  for(i in seq_along(temp)){
    
    dat     <- temp[[i]][, ..take.columns]
    dat.med <- medpolish(dat, maxiter = 3, na.rm = TRUE)
    dat.res <- dat.med$residuals
    
    boxplot(dat, main = paste0("Box plot ", names(temp)[i], " before normalization"), xaxt = "n")
    axis(side = 1, at = 1:6, labels = paste0("X", 1:6))
    boxplot(dat.res, main = paste0("Box plot ", names(temp)[i], " after med polishing"), xaxt = "n")
    axis(side = 1, at = 1:6, labels = paste0("X", 1:6))
    
    list_densities <- lapply(dat, density, na.rm = TRUE)
    xmax <- max(unlist(lapply(list_densities, function(dens) max(dens$x))))
    xmin <- min(unlist(lapply(list_densities, function(dens) min(dens$x))))
    ymax <- max(unlist(lapply(list_densities, function(dens) max(dens$y))))
    ymin <- min(unlist(lapply(list_densities, function(dens) min(dens$y))))
    
    plot(list_densities[[1]], type = "n",
         xlim  = c(floor(xmin), ceiling(xmax)),
         ylim = c(floor(ymin), ymax),
         main = paste0("Density plot ", names(temp)[i]))
    for(j in seq_along(list_densities)){
      
      lines(list_densities[[j]], col = line_colors[j], lty = line_types[j], lwd = 2)
      
    }
    legend(x = "topright", legend = paste0("X", 1:6), fill = line_colors)
    
    temp[[i]][, paste0(take.columns, ".norm") := data.table(..dat.res)]
    
  }
  dev.off()
  
  # compute contrasts
  temp <- rbindlist(temp)
  temp[, T1 := T1 - T0]
  temp[, T2 := T2 - T0]
  temp[, T1.norm := T1.norm - T0.norm]
  temp[, T2.norm := T2.norm - T0.norm]
  
  temp <- temp[, -c("T0", "T0.norm")]
  names(temp)
  
  # wide format
  temp <- split(temp, f = temp$Set)
  temp[[1]] <- merge(temp[[1]][, -c("Set")], temp[[2]][, c("id", "Replicate", "Stimulation", "T1", "T2", "T1.norm", "T2.norm")], by = c("id", "Replicate", "Stimulation"), suffixes = c("", ".2"))
  temp <- temp[[1]]
  names(temp)[names(temp) == "T1.2"] <- "T3"
  names(temp)[names(temp) == "T2.2"] <- "T4"
  names(temp)[names(temp) == "T1.norm.2"] <- "T3.norm"
  names(temp)[names(temp) == "T2.norm.2"] <- "T4.norm"
  
  # take average of each replicate
  temp <- temp[, lapply(.SD, mean, na.rm = TRUE), by = c("id", "Stimulation"), .SDcols = c("T1", "T2", "T3", "T4",
                                                                                           "T1.norm", "T2.norm",
                                                                                           "T3.norm", "T4.norm")]
  
  temp <- merge(ph[, c("id", "Accession", "Gene.name", "Position",
                       "Sequence.window_7", "Sequence.window_15",
                       "Sequence.window_7_noSpace")], 
                temp, by = "id")
  
  fwrite(temp, "comparison_EK\\temp\\PhSiteContrasts_EK.tsv", sep = "\t")

})