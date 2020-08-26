# Do:
#
# Fisher's exact test of the netphorest group enrichment
# between CaEGTA and MockBoNT experiments
#
# INPUT:
# "temp\\PhPeptIntensities_slim.tsv"
#
# OUTPUT:
# "Tables\\SupplTable_3_netphorest_group_CaEGTA_MockBoNT_Fisher.txt"

local({
  
  if(!dir.exists("Tables")) dir.create("Tables")
  
  library(data.table)
  
  # Phoshpo data
  ph <- fread("temp\\PhPeptIntensities_slim.tsv")
  
  # split netphorest kinase groups into single elements
  ph[, netphorest_group := str_split(netphorest_group, ";")]
  
  # extract the second outmost kinase group
  ph[, netphorest_group2 := lapply(netphorest_group, function(x) x[max(c(length(x) - 1, 2))])]
  ph[, netphorest_group2 := unlist(netphorest_group2)]
  ph[, netphorest_group2 := gsub("_group", "", netphorest_group2)]
  
  # extract the highest level kinase group
  ph[, netphorest_group_high := unlist(lapply(netphorest_group, function(x) x[length(x)]))]
  
  ph_cand_Ca  <- ph[ph$Candidate.CaEGTA & !is.na(netphorest_group2)]
  ph_cand_BoNT <- ph[ph$Candidate.BoNT & !is.na(netphorest_group2)]
  
  count_Ca  <- ph_cand_Ca[, .N, by = netphorest_group2]
  count_Ca  <- count_Ca[, All := nrow(ph[ph$Candidate.CaEGTA & ! is.na(netphorest_group2)])]
  
  count_BoNT <- ph_cand_BoNT[, .N, by = netphorest_group2]
  count_BoNT <- count_BoNT[, All := nrow(ph[ph$Candidate.BoNT & ! is.na(netphorest_group2)])]
  
  count <- merge(count_Ca, count_BoNT, by = "netphorest_group2", all = TRUE, suffixes = c(".Ca", ".BoNT"))
  
  # replace NAs either by 0 or by the total count
  count[is.na(All.Ca), All.Ca := count$All.Ca[!is.na(count$All.Ca)][1]]
  count[is.na(All.BoNT), All.BoNT := count$All.BoNT[!is.na(count$All.BoNT)][1]]
  
  count[is.na(N.Ca), N.Ca := 0]
  count[is.na(N.BoNT), N.BoNT := 0]
  
  count <- count[count$N.Ca > 5 | count$N.BoNT > 5]
  
  # use Fisher exact test to compute p values
  for(i in 1:nrow(count)){
    
    count[i, p.val := NA_real_]
    mtx <- matrix(c(count$N.Ca[i], count$All.Ca[i],
                    count$N.BoNT[i], count$All.BoNT[i]),
                  byrow = TRUE, ncol = 2)
    ft <- tryCatch({fisher.test(x = mtx, alternative = "two.sided")$p.value}, error = function(e) NA_real_)
    count[i, p.val := ft]
    
  }
  
  count[, p.adj.BH := p.adjust(count$p.val, "BH")]
  count[, Percent_Ca  := 100*N.Ca/All.Ca]
  count[, Percent_BoNT := 100*N.BoNT/All.BoNT]
  count[, Percent_diff := Percent_Ca - Percent_BoNT]
  count <- count[order(-Percent_diff)]
  count[p.adj.BH < 0.01]
  
  count <- count[, c("netphorest_group2", "N.Ca", "All.Ca", "N.BoNT", "All.BoNT", "p.val", "p.adj.BH")]
  names(count) <- c("Kinase group", "CaEGTA", "Total CaEGTA", "MockBoNT", "Total MockBoNT", "p.val", "p.adj")

  fwrite(count, "Tables\\SupplTable_3_netphorest_group_CaEGTA_MockBoNT_Fisher.txt", sep = "\t")
  
})

