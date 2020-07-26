# Do:
#
# Fisher's exact test of the netphorest group enrichment
# use Networkin predictions for human proteome as background
#
# INPUT:
# "temp\\PhPeptIntensities_slim.tsv"
# "external\\data_NetworKIN\\networkin_human_predictions_3.1.tsv.bz2"

# OUTPUT:
# "SupplTable_1_2_netphorest_group_human_bcgr_fisher.txt"


local({ 
  
  library(data.table)
  
  # phosphosites data
  ph <- fread("temp\\PhPeptIntensities_slim.tsv", check.names = TRUE)

  # split netphorest kinase groups into single elements
  ph[, netphorest_group := str_split(netphorest_group, ";")]
  
  # extract the second outmost kinase group
  ph[, netphorest_group2 := lapply(netphorest_group, function(x) x[max(c(length(x) - 1, 2))])]
  ph[, netphorest_group2 := unlist(netphorest_group2)]
  ph[, netphorest_group2 := gsub("_group", "", netphorest_group2)]
  
  # extract the highest level kinase group
  ph[, netphorest_group_high := unlist(lapply(netphorest_group, function(x) x[length(x)]))]
  
  # number of ph sites regulated by a particular kinase group
  npgroup_count_all <- ph[!is.na(ph$netphorest_group2), list(In_all = .N), by = netphorest_group2]
  
  # number of all sites - the number of sites regulated by a particular kinase group
  npgroup_count_all[, Out_all := ph[!is.na(ph$netphorest_group), .N] - In_all]

  ####### Mock / BTX ##############
  
  # subset candidates with annotated netphorest groups
  ph_sub <- ph[ph$Candidate.BTX & !is.na(ph$netphorest_group2)]
  ph_sub[, netphorest_group2 := sub("_group$", "", netphorest_group2)]
  
  unique(ph_sub$netphorest_group2)
  
  npgroup_count  <- ph_sub[, list(N_count = .N), by = netphorest_group2]
  npgroup_count <- npgroup_count[order(N_count), ]
  npgroup_count[, Out := ph_sub[, .N] - N_count]

  # compare netphorest group enrichment to human data set
  np <- fread("temp\\netphorest_groups.tsv", check.name = TRUE)
  
  # extract the second outmost kinase group
  np[, netphorest_group := str_split(Group_all, ";")]
  np[, netphorest_group2 := lapply(netphorest_group, function(x) x[max(c(length(x) - 1, 2))])]
  np[, netphorest_group2 := unlist(netphorest_group2)]
  np[, netphorest_group2 := gsub("_group", "", netphorest_group2)]
  
  # update netphorest group based on consensus kinase
  # use the most specific group
  # np is sorted for more specific groups being in the upper part of the table
  np <- np[!duplicated(np[, c("Gene_id")])]
  dim(np)
  
  # load networkin predictions for human proteins
  nk_pred <- fread("external\\data_NetworKIN\\networkin_human_predictions_3.1.tsv.bz2", check.names = TRUE)
  nk_pred <- nk_pred[nk_pred$tree == "KIN"]
  nk_pred <- merge(nk_pred, np[, c("Gene_id", "netphorest_group2")], by.x = c("id"), by.y = "Gene_id", all.x = TRUE)
  nk_pred <- nk_pred[nk_pred$netphorest_group2 %in% npgroup_count$netphorest_group2]
  nk_pred <- nk_pred[order(-networkin_score)]
  nk_pred <- nk_pred[nk_pred$networkin_score > 0.5]
  nk_pred <- nk_pred[!duplicated(nk_pred[, c("sequence", "X.substrate")])]
  table(nk_pred$netphorest_group2)
  round((100*table(nk_pred$netphorest_group2))/nrow(nk_pred), 2)
  
  npgroup_count_bgr <- nk_pred[, list(In_all = .N), by = netphorest_group2]
  npgroup_count_bgr[, Out_all := nk_pred[, .N] - In_all]
  
  npgroup_count <- merge(npgroup_count, npgroup_count_bgr, by = "netphorest_group2")
  npgroup_count <- npgroup_count[!is.na(netphorest_group2)]
  
  # prepare matrices for testing
  cont_tables <- list()
  for(i in 1:nrow(npgroup_count)){
    
    cont_tables[[i]] <- matrix(c(npgroup_count$N_count[i],
                                 npgroup_count$Out[i],
                                 npgroup_count$In_all[i],
                                 npgroup_count$Out_all[i]),
                               nrow = 2, byrow = TRUE,
                               dimnames = list(Candidate = c("In", "Out"), Background = c("In", "Out")))
    
  }
  names(cont_tables) <- npgroup_count$netphorest_group2
  
  # perform Fisher's exact test  
  p_value <- numeric()
  for(i in seq_along(cont_tables)){
    
    p_value[i] <- tryCatch({fisher.test(cont_tables[[i]])$p.value}, error = function(e) NA)
    
  }
  
  BTX <- npgroup_count
  BTX[, Percent_data := 100*N_count/Out]
  BTX[, Percent_all := 100*In_all/Out_all]
  BTX[, Percent_diff := Percent_data - Percent_all]
  BTX <- cbind(BTX, p.val = p_value)
  BTX[, p.adj := p.adjust(p.val, "BH")]
  BTX[, Experiment := "BTX"]

  ########## CaEGTA #################
  
  # subset candidates with annotated netphorest groups
  ph_sub <- ph[ph$Candidate.CaEGTA & !is.na(ph$netphorest_group2)]
  ph_sub[, netphorest_group2 := sub("_group$", "", netphorest_group2)]
  
  unique(ph_sub$netphorest_group2)
  
  npgroup_count  <- ph_sub[, list(N_count = .N), by = netphorest_group2]
  npgroup_count <- npgroup_count[order(N_count), ]
  npgroup_count[, Out := ph_sub[, .N] - N_count]
  
  # compare netphorest group enrichment to human data set
  np <- fread("temp\\netphorest_groups.tsv", check.name = TRUE)
  
  # extract the second outmost kinase group
  np[, netphorest_group := str_split(Group_all, ";")]
  np[, netphorest_group2 := lapply(netphorest_group, function(x) x[max(c(length(x) - 1, 2))])]
  np[, netphorest_group2 := unlist(netphorest_group2)]
  np[, netphorest_group2 := gsub("_group", "", netphorest_group2)]
  
  # update netphorest group based on consensus kinase
  # use the most specific group
  # np is sorted for more specific groups being in the upper part of the table
  np <- np[!duplicated(np[, c("Gene_id")])]
  dim(np)

  # load networkin predictions for human proteins
  nk_pred <- fread("external\\data_NetworKIN\\networkin_human_predictions_3.1.tsv.bz2", check.names = TRUE)
  nk_pred <- nk_pred[nk_pred$tree == "KIN"]
  nk_pred <- merge(nk_pred, np[, c("Gene_id", "netphorest_group2")], by.x = c("id"), by.y = "Gene_id", all.x = TRUE)
  nk_pred <- nk_pred[nk_pred$netphorest_group2 %in% npgroup_count$netphorest_group2]
  nk_pred <- nk_pred[order(-networkin_score)]
  nk_pred <- nk_pred[nk_pred$networkin_score > 0.5]
  nk_pred <- nk_pred[!duplicated(nk_pred[, c("sequence", "X.substrate")])]
  table(nk_pred$netphorest_group2)
  round((100*table(nk_pred$netphorest_group2))/nrow(nk_pred), 2)
  
  npgroup_count_bgr <- nk_pred[, list(In_all = .N), by = netphorest_group2]
  npgroup_count_bgr[, Out_all := nk_pred[, .N] - In_all]
  
  npgroup_count <- merge(npgroup_count, npgroup_count_bgr, by = "netphorest_group2")
  npgroup_count <- npgroup_count[!is.na(netphorest_group2)]
  
  # prepare matrices for testing
  cont_tables <- list()
  for(i in 1:nrow(npgroup_count)){
    
    cont_tables[[i]] <- matrix(c(npgroup_count$N_count[i],
                                 npgroup_count$Out[i],
                                 npgroup_count$In_all[i],
                                 npgroup_count$Out_all[i]),
                               nrow = 2, byrow = TRUE,
                               dimnames = list(Candidate = c("In", "Out"), Background = c("In", "Out")))
    
  }
  names(cont_tables) <- npgroup_count$netphorest_group2
  
  
  p_value <- numeric()
  for(i in seq_along(cont_tables)){
    
    p_value[i] <- tryCatch({fisher.test(cont_tables[[i]])$p.value}, error = function(e) NA)
    
  }
  
  CaEGTA <- npgroup_count
  CaEGTA[, Percent_data := 100*N_count/Out]
  CaEGTA[, Percent_all := 100*In_all/Out_all]
  CaEGTA[, Percent_diff := Percent_data - Percent_all]
  CaEGTA <- cbind(CaEGTA, p.val = p_value)
  CaEGTA[, p.adj := p.adjust(p.val, "BH")]
  CaEGTA[, Experiment := "CaEGTA"]
  
  # combine tables
  df <- rbind(CaEGTA, BTX)
  df[, p.adj := p.adjust(p.val, "BH")]
  df <- df[order(-Experiment, -Percent_diff)]
  df <- df[, c("netphorest_group2", "N_count", "Out", "In_all", "Out_all", "p.val", "p.adj", "Experiment")]
  names(df) <- c("Kinase group", "Sites in", "Sites out", "Sites in bcgr", "Sites out bcgr", "p.val", "p.adj", "Experiment")
  
  fwrite(df, "SupplTable_1_2_netphorest_group_human_bcgr_fisher.txt", sep = "\t")

})