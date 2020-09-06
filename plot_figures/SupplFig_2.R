
local({
  
  library(data.table)
  library(ggplot2)
  
  ph <- fread("Figures\\SupplFig_2\\Ph_events_netphorest.txt")
  
  # split netphorest kinase groups into single elements
  ph[, netphorest_group := str_split(netphorest_group, ";")]
  
  # extract the second outmost kinase group
  ph[, netphorest_group2 := lapply(netphorest_group, function(x) x[max(c(length(x) - 1, 2))])]
  ph[, netphorest_group2 := unlist(netphorest_group2)]
  ph[, netphorest_group2 := gsub("_group", "", netphorest_group2)]
  
  # extract the highest level kinase group
  ph[, netphorest_group_high := unlist(lapply(netphorest_group, function(x) x[length(x)]))]
  
  # number of ph events regulated by a particular kinase group
  npgroup_count_all <- ph[!is.na(ph$netphorest_group2), list(In_all = .N), by = netphorest_group2]
  
  # number of all sites - the number of sites regulated by a particular kinase group
  npgroup_count_all[, Out_all := ph[!is.na(ph$netphorest_group), .N] - In_all]
  
  #### CaEGTA Experiment ####
  
  # subset candidates from Ca experiment with annotated netphorest groups
  ph_Ca <- ph[ph$Candidate.CaEGTA & !is.na(ph$netphorest_group)]
  ph_Ca[, Site     :=  paste0(Gene.name, "_", Amino.acid, Position)]
  ph_Ca[, Regulation := "Down"]
  ph_Ca[log2FC.CaEGTA > 0, Regulation := "Up"]
  
  # number of candidate sites by netphorest group
  npgroup_count_Ca  <- ph_Ca[, list(N_count = .N,
                                    N_up    = sum(Regulation == "Up"),
                                    N_down  = sum(Regulation == "Down")), by = netphorest_group2]
  npgroup_count_Ca <- npgroup_count_Ca[order(N_count), ]
  npgroup_count_Ca[, Out := ph_Ca[, .N] - N_count]
  
  # prepare for plotting #
  
  # order netphorest groups
  npgroup_count_Ca <- npgroup_count_Ca[order(N_count)]
  npgroup_count_Ca <- npgroup_count_Ca[!is.na(netphorest_group2)]
  
  # order by high level netphorest groups
  npgroup_count_high_Ca  <- ph_Ca[, list(N_count = .N,
                                         N_up    = sum(Regulation == "Up"),
                                         N_down  = sum(Regulation == "Down")), by = c("netphorest_group_high", "netphorest_group2")]
  npgroup_count_high_Ca <- npgroup_count_high_Ca[order(-N_count), ]
  unique(npgroup_count_high_Ca$netphorest_group_high)
  temp <- split(npgroup_count_high_Ca, npgroup_count_high_Ca$netphorest_group_high, sorted = FALSE)[unique(npgroup_count_high_Ca$netphorest_group_high)]
  npgroup_count_high_Ca <- rbindlist(temp)
  npgroup_count_high_Ca <- npgroup_count_high_Ca[nrow(npgroup_count_high_Ca) : 1]
  
  # sort netphorest groups based on the sorting of the higher kinase groups
  ph_Ca[, netphorest_group2 := factor(netphorest_group2, levels = unique(npgroup_count_high_Ca$netphorest_group2))]
  
  #### Mock/BoNT experiment ####
  
  # subset candidates from Ca experiment with annotated netphorest groups
  ph_BoNT <- ph[ph$Candidate.BoNT & !is.na(ph$netphorest_group)]
  ph_BoNT[, Site     :=  paste0(Gene.name, "_", Amino.acid, Position)]
  ph_BoNT[, Regulation := "Down"]
  ph_BoNT[log2FC.BoNT > 0, Regulation := "Up"]
  
  # number of candidate sites by netphorest group
  npgroup_count_BoNT  <- ph_BoNT[, list(N_count = .N,
                                      N_up    = sum(Regulation == "Up"),
                                      N_down  = sum(Regulation == "Down")), by = netphorest_group2]
  npgroup_count_BoNT <- npgroup_count_BoNT[order(N_count), ]
  npgroup_count_BoNT[, Out := ph_BoNT[, .N] - N_count]
  
  # prepare for plotting #
  
  # order netphorest groups
  npgroup_count_BoNT <- npgroup_count_BoNT[order(N_count)]
  npgroup_count_BoNT <- npgroup_count_BoNT[!is.na(netphorest_group2)]
  
  # order by high level netphorest groups
  npgroup_count_high_BoNT  <- ph_BoNT[, list(N_count = .N,
                                           N_up    = sum(Regulation == "Up"),
                                           N_down  = sum(Regulation == "Down")), by = c("netphorest_group_high", "netphorest_group2")]
  npgroup_count_high_BoNT <- npgroup_count_high_BoNT[order(-N_count), ]
  unique(npgroup_count_high_BoNT$netphorest_group_high)
  temp <- split(npgroup_count_high_BoNT, npgroup_count_high_BoNT$netphorest_group_high, sorted = FALSE)[unique(npgroup_count_high_BoNT$netphorest_group_high)]
  npgroup_count_high_BoNT <- rbindlist(temp)
  npgroup_count_high_BoNT <- npgroup_count_high_BoNT[nrow(npgroup_count_high_BoNT) : 1]
  
  # sort netphorest groups based on the sorting of the higher kinase groups
  ph_BoNT[, netphorest_group2 := factor(netphorest_group2, levels = unique(npgroup_count_high_BoNT$netphorest_group2))]
  
  ##### plot CaEGTA #####
  
  max_y_limit <- 1.1*max(c(npgroup_count_Ca$N_count, npgroup_count_BoNT$N_count))
  x_side <- round(max_y_limit/nrow(npgroup_count_Ca), 0) / 2
  
  # higher group annotation
  temp_Ca <- data.table(label = c("PDHK",
                                  "STE",
                                  "CK",
                                  "AGC group",
                                  "CAMK group", 
                                  "CMGC group"),
                        xmin = c(5, 7, 10, 12, 25, 34),
                        xmax = c(6, 9, 11, 24, 33, 39),
                        ymin = rep(max(npgroup_count_Ca$N_count) -10, 6),
                        ymax = rep(max(npgroup_count_Ca$N_count) -10, 6))
  temp_Ca[, col := c("#ababab", "#71c0a8", "#ae6161", "#a7a4cc", "#ec79b4", "#9ec473")]
  
  pdf("Figures\\SupplFig_2\\SupplFig_2A.pdf", width = 13, heigh = 8)
  
  g <- ggplot(ph_Ca[!(is.na(netphorest_group) | netphorest_group == "")], aes(x = netphorest_group2))
  g <- g + geom_bar(aes(fill = Regulation))
  g <- g + scale_fill_manual(breaks = c("Up", "Down"), values = c("#51baf4", "orange"))
  g <- g + scale_y_continuous(expand = c(0, 1), limits = c(-x_side, max_y_limit))
  g <- g + scale_x_discrete(labels = npgroup_count_high_Ca$netphorest_group2)
  g <- g + coord_flip()
  g <- g + labs(title = "Netphorest Groups")
  g <- g + xlab("") + ylab("Number of regulated phosphorylation sites")
  g <- g + theme_bw()
  g <- g + theme(axis.text.y = element_text(face = "bold"),
                 plot.margin = unit(c(1, 1, 1, 3.4), "cm")
  )
  
  g <- g + annotate("rect",
                    xmin = temp_Ca$xmin-0.4,
                    xmax = temp_Ca$xmax+0.4,
                    ymin = -x_side,
                    ymax = -1,
                    fill = temp_Ca$col)
  
  g <- g + annotate("rect",
                    xmin = temp_Ca$xmin-0.4,
                    xmax = temp_Ca$xmax+0.4,
                    ymin = max_y_limit - 10,
                    ymax = max_y_limit - 5,  
                    fill = scales::alpha(temp_Ca$col, 0.6))
  
  
  g <- g + annotate("text", x = (temp_Ca$xmin + temp_Ca$xmax)/2,
                    y = max_y_limit - 10,
                    label = temp_Ca$label,
                    size = 3,
                    fontface = "bold",
                    angle = -90,
                    hjust = 0.5,
                    vjust = -0.5)
  
  print(g)
  
  dev.off()
  
  ##### plot NoTox/Tox #####
  
  temp_BoNT <- data.table(label = c("STE",
                                   "CAMK group",
                                   "CK",
                                   "AGC group", 
                                   "CMGC group"),
                         xmin = c(5, 8, 14, 16, 27),
                         xmax = c(7, 13, 15, 26, 33),
                         ymin = rep(max(npgroup_count_BoNT$N_count) -10, 5),
                         ymax = rep(max(npgroup_count_BoNT$N_count) -10, 5))
  temp_BoNT <- merge(temp_BoNT, temp_Ca[, c("label", "col")], by = "label", all.x = TRUE)
  
  
  pdf("Figures\\SupplFig_2\\SupplFig_2B.pdf", width = 13, heigh = 8)
  
  g <- ggplot(ph_BoNT[!(is.na(netphorest_group) | netphorest_group == "")], aes(x = netphorest_group2))
  g <- g + geom_bar(aes(fill = Regulation))
  g <- g + scale_fill_manual(breaks = c("Up", "Down"), values = c("#51baf4", "orange"))
  g <- g + scale_y_continuous(expand = c(0, 1), limits = c(-x_side, max_y_limit))
  g <- g + scale_x_discrete(labels = npgroup_count_high_BoNT$netphorest_group2)
  g <- g + coord_flip()
  g <- g + labs(title = "Netphorest Groups")
  g <- g + xlab("") + ylab("Number of regulated phosphorylation sites")
  g <- g + theme_bw()
  g <- g + theme(axis.text.y = element_text(face = "bold"),
                 plot.margin = unit(c(1, 1, 1, 0), "cm")
  )
  g <- g + annotate("rect",
                    xmin = temp_BoNT$xmin-0.4,
                    xmax = temp_BoNT$xmax+0.4,
                    ymin = -x_side,
                    ymax = -1,
                    fill = temp_BoNT$col)
  
  g <- g + annotate("rect",
                    xmin = temp_BoNT$xmin-0.4,
                    xmax = temp_BoNT$xmax+0.4,
                    ymin = max_y_limit - 10,
                    ymax = max_y_limit - 5,  
                    fill = scales::alpha(temp_BoNT$col, 0.6))
  
  
  g <- g + annotate("text", x = (temp_BoNT$xmin + temp_BoNT$xmax)/2,
                    y = max_y_limit - 10,
                    label = temp_BoNT$label,
                    size = 3,
                    fontface = "bold",
                    angle = -90,
                    hjust = 0.5,
                    vjust = -0.5)
  
  print(g)
  
  dev.off()

})
