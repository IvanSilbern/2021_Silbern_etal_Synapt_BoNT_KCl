# Do:
#
# Test enrichment of kinasee-substrate relationships and GO biological function terms
#
# INPUT:
# "temp\\PhPeptIntensities_slim.tsv"
# "external\\data_DAVID\\DAVID_GO_BP_CaEGTA.txt"
# (GO BP annotation of proteins carrying regulated phosphorylation sites in CaEGTA)
#
# OUTPUT:
# "Figures\\Fig_3A\\NPh_group_GOBP.txt"
# "plots\\NPh_group_GOBP_pval.pdf"

local({
  
  if(!dir.exists("Figures\\Fig_3A")) dir.create("Figures\\Fig_3A", recursive = TRUE)
  if(!dir.exists("plots")) dir.create("plots")
  
  library(data.table)
  library(ggplot2)
  
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
  
  ph_cand <- ph[ph$Candidate.CaEGTA & !is.na(netphorest_group2)]
  
  count_groups <- unique(ph_cand[, c("Gene.name", "netphorest_group2", "netphorest_group_high")])
  count_groups <- count_groups[, .N, by = c("netphorest_group2", "netphorest_group_high")]
  count_groups <- count_groups[order(-N)]
  count_groups <- count_groups[order(netphorest_group_high, -N)]
  
  take_groups  <- count_groups[N >= 5, ]
  take_groups  <- take_groups[take_groups$netphorest_group2 != ""]

  # GOBP
  gobp    <- fread("external\\data_DAVID\\DAVID_GO_BP_CaEGTA.txt", check.names = TRUE)
  gobp[, Term := gsub("^.+~", "", Term)] # remove GO:...~annotation
    
  gobp <- gobp[order(Benjamini)]
  gobp <- gobp[gobp$Benjamini < 0.001, ]
  dim(gobp)

  gobp[, gene_names := strsplit(gobp$Genes, split = ", ")]
  
  unique(gobp$Term)

# Fisher exact test on GO groups
# test if sites regulated by particular np group and belonging to one of go terms
# are significantly significantly enriched as compared to the occurance of the sites regulated by np group
# among all candidate sites

# Nr. of sites regulated by NP kinase group  | Total Number of sites in GO term
# for a given GO Term
# -------------------------------------------|------------------------------------- 
# Nr. of sites regulated by NP kinase group  | Nr. of all candidate sites
# among all candidate sites                  |          

  temp <- data.table(netph_group = c(take_groups$netphorest_group2, "Other"))
  for(i in seq_along(gobp$gene_names)){
    
    count <- ph_cand[toupper(ph_cand$Gene.name) %in% gobp$gene_names[[i]], .N, by = netphorest_group2]
    other <- count[!count$netphorest_group2 %in% take_groups$netphorest_group2]
    count <- count[count$netphorest_group2 %in% take_groups$netphorest_group2]
    count <- rbind(count, data.table(netphorest_group2 = "Other", N = sum(other$N, na.rm = TRUE)))
    
    temp <- merge(temp, count[, c("netphorest_group2", "N")],
                  by.x = "netph_group", by.y = "netphorest_group2",
                  all.x = TRUE, suffixes = c("", paste0(".", i)))
    
  }
  
  # add GO Term names as column names
  names(temp)[-1] <- gobp$Term

  # total number of sites per GO Term
  count_total <- apply(temp[, -1], 2, sum, na.rm = TRUE)
  
  # number of candidate sites in total
  sites_total <- nrow(ph_cand)
  
  # number of sites per netphorest group
  sites_np <- ph_cand[, .N, by = netphorest_group2]
  other <- sites_np[!sites_np$netphorest_group2 %in% take_groups$netphorest_group2]
  sites_np <- sites_np[sites_np$netphorest_group2 %in% take_groups$netphorest_group2]
  sites_np <- rbind(sites_np, data.table(netphorest_group2 = "Other", N = sum(other$N, na.rm = TRUE)))

  # present as a matrix and apply Fisher exact testing
  
  counts <- data.table()
  for(i in names(temp[, -1])){
    
    for(j in temp$netph_group){
    
    counts <- rbind(counts, data.table(go_term = i,
                                       np_group = j,
                                       count = unlist(temp[netph_group == j, ..i]),
                                       count_total = count_total[names(count_total) == i],
                                       sites_np = unlist(sites_np[netphorest_group2 == j, "N"]),
                                       sites_total = sites_total))
    
    }
    
  }
  
  for(i in 1:nrow(counts)){
    
    counts[i, p.val := NA_real_]
    mtx <- matrix(c(counts$count[i], counts$count_total[i],
                    counts$sites_np[i], counts$sites_total[i]),
                    byrow = TRUE, ncol = 2)
    ft <- tryCatch({fisher.test(x = mtx, alternative = "greater")$p.value}, error = function(e) NA_real_)
    counts[i, p.val := ft]
      
  }

  counts[, np_group := factor(np_group, levels = c(take_groups$netphorest_group2, "Other"))]
  counts[, GO_Term := factor(go_term, levels = gobp$Term)]

  # provide source file
  
  fwrite(counts, "Figures\\Fig_3A\\NPh_group_GOBP.txt", sep = "\t")
  fwrite(take_groups, "Figures\\Fig_3A\\netphorest_groups.txt", sep = "\t")
  fwrite(gobp, "Figures\\Fig_3A\\GO_BP.txt", sep = "\t")
  
  pdf("plots\\NPh_group_GOBP_pval.pdf", width = 12 , heigh = 10)
  g <- ggplot(counts, aes(np_group, GO_Term))
  
  g <- g + geom_tile(aes(fill = p.val), colour = "grey20")
  g <- g + scale_fill_gradient(low = "orange", high = "steelblue", na.value = "white", limits = c(0, 0.20), oob = scales::squish,
                               breaks = c(0, 0.05, 0.1, 0.15, 0.2), label = c("  0.00", "  0.05", "  0.10", "  0.15", "> 0.20"))
  g <- g + scale_x_discrete(expand = c(0, 0), position = "top", limits = c(levels(counts$np_group)))
  g <- g + scale_y_discrete(expand = c(0, 0), limits = rev(levels(counts$GO_Term)))
  g <- g + theme(axis.text.x = element_text(angle = 30, hjust = 0, face = "bold", size = 16),
                 axis.text.y = element_text(size = 16),
                 legend.key.size = unit(1, "cm"),
                 legend.text = element_text(size = 14),
                 legend.title = element_text(size = 18))
  g <- g + xlab("") + ylab("")
  g <- g + guides(fill = guide_colorbar(title = "p-Value"))
  print(g)
  dev.off()

})
