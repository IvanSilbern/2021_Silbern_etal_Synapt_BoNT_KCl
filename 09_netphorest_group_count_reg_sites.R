
# DO:
# calculate number of regulated events per netphorest group
# plot
# assess singificance in the number of down/up regulated sites in each group

# INPUT:
# "temp\\PhPeptIntensities_slim.tsv"

# OUTPUT:
# "Figures\\SupplFig_6\\Ph_events_netphorest.txt"
# "Figures\\Fig_2EF\\Ph_events_netphorest.txt"
# "plots\\NetphorestGroups_Count_CaEGTA.pdf"
# "plots\\NetphorestGroups_Count_BTX.pdf"
# "plots\\NetphorestGroups_Count_CaEGTA_red.pdf"
# "plots\\NetphorestGroups_Count_BTX_red.pdf"

local({
  
  if(!dir.exists("plots")) dir.create("plots")
  if(!dir.exists("Figures\\SupplFig_6")) dir.create("Figures\\SupplFig_6", recursive = TRUE)
  if(!dir.exists("Figures\\Fig_2EF")) dir.create("Figures\\Fig_2EF", recursive = TRUE)
  
  library(data.table)
  library(ggplot2)
  
  # phosphosites data
  ph <- fread("temp\\PhPeptIntensities_slim.tsv", check.names = TRUE)
  
  # provide source data
  fwrite(ph[, c("id", "Gene.name", "Accession",
                "Amino.acid", "Position", "Multiplicity",
                "q.val.CaEGTA", "log2FC.CaEGTA", "Candidate.CaEGTA",
                "q.val.BTX", "log2FC.BTX", "Candidate.BTX", "netphorest_group")],
         "Figures\\SupplFig_6\\Ph_events_netphorest.txt",
         sep = "\t")
  
  fwrite(ph[, c("id", "Gene.name", "Accession",
                "Amino.acid", "Position", "Multiplicity",
                "q.val.CaEGTA", "log2FC.CaEGTA", "Candidate.CaEGTA",
                "q.val.BTX", "log2FC.BTX", "Candidate.BTX", "netphorest_group")],
         "Figures\\Fig_2EF\\Ph_events_netphorest.txt",
         sep = "\t")

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

  # test changes in the number of up/down regulated sites
  
  # prepare matrices (contigency tables) for testing
  cont_tables <- list()
  for(i in 1:nrow(npgroup_count_Ca)){
    
    cont_tables[[i]] <- matrix(c(npgroup_count_Ca$N_up[i],
                                 npgroup_count_Ca$N_down[i],
                                 sum(ph$log2FC.CaEGTA[!is.na(ph$log2FC.CaEGTA)] > 0),
                                 sum(ph$log2FC.CaEGTA[!is.na(ph$log2FC.CaEGTA)] < 0)),
                               nrow = 2, byrow = TRUE,
                               dimnames = list(Candidate = c("In", "Out"), Background = c("In", "Out")))
    
  }
  names(cont_tables) <- npgroup_count_Ca$netphorest_group2
  
  # Fisher's exact test; store p-values
  p_value <- numeric()
  for(i in seq_along(cont_tables)){
    
    p_value[i] <- fisher.test(cont_tables[[i]])$p.value
    
  }
  
  # Benjamini-Hochberg adjustment
  p.adj.updown <- p.adjust(p_value, "BH")
  
  # significant
  deregulated_updown <- names(cont_tables[which(p.adj.updown < 0.01)])

  # prepare for plotting #
  
  # order netphorest groups
  npgroup_count_Ca <- npgroup_count_Ca[order(N_count)]
  npgroup_count_Ca <- npgroup_count_Ca[!is.na(netphorest_group2)]
  
  # mark kinase groups that show significant deregulation of up or down-phosphorylated sites
  npgroup_count_Ca[, np_marked := netphorest_group2]
  npgroup_count_Ca[npgroup_count_Ca$netphorest_group2 %in% deregulated_updown, np_marked := paste0("* ", netphorest_group2)]
  
  # order by high level netphorest groups
  npgroup_count_high_Ca  <- ph_Ca[, list(N_count = .N,
                                       N_up    = sum(Regulation == "Up"),
                                       N_down  = sum(Regulation == "Down")), by = c("netphorest_group_high", "netphorest_group2")]
  npgroup_count_high_Ca <- npgroup_count_high_Ca[order(-N_count), ]
  unique(npgroup_count_high_Ca$netphorest_group_high)
  temp <- split(npgroup_count_high_Ca, npgroup_count_high_Ca$netphorest_group_high, sorted = FALSE)[unique(npgroup_count_high_Ca$netphorest_group_high)]
  npgroup_count_high_Ca <- rbindlist(temp)
  npgroup_count_high_Ca <- npgroup_count_high_Ca[nrow(npgroup_count_high_Ca) : 1]
  npgroup_count_high_Ca[, np_marked := netphorest_group2]
  npgroup_count_high_Ca[npgroup_count_high_Ca$netphorest_group2 %in% deregulated_updown, np_marked := paste0("* ", netphorest_group2)]
  
  # sort netphorest groups based on the sorting of the higher kinase groups
  ph_Ca[, netphorest_group2 := factor(netphorest_group2, levels = unique(npgroup_count_high_Ca$netphorest_group2))]
  

#### Mock/BTX experiment ####
  
  # subset candidates from Ca experiment with annotated netphorest groups
  ph_BTX <- ph[ph$Candidate.BTX & !is.na(ph$netphorest_group)]
  ph_BTX[, z_scores := (log2FC.BTX - mean(log2FC.BTX)) / sd(log2FC.BTX)]
  ph_BTX[, Site     :=  paste0(Gene.name, "_", Amino.acid, Position)]
  ph_BTX[, Regulation := "Down"]
  ph_BTX[log2FC.BTX > 0, Regulation := "Up"]
  
  # number of candidate sites by netphorest group
  npgroup_count_BTX  <- ph_BTX[, list(N_count = .N,
                                    N_up    = sum(Regulation == "Up"),
                                    N_down  = sum(Regulation == "Down")), by = netphorest_group2]
  npgroup_count_BTX <- npgroup_count_BTX[order(N_count), ]
  npgroup_count_BTX[, Out := ph_BTX[, .N] - N_count]


  # test changes in number of up/down regulated sites #####
  
  cont_tables <- list()
  for(i in 1:nrow(npgroup_count_BTX)){
    
    cont_tables[[i]] <- matrix(c(npgroup_count_BTX$N_up[i],
                                 npgroup_count_BTX$N_down[i],
                                 sum(ph$log2FC.BTX[!is.na(ph$log2FC.BTX)] > 0),
                                 sum(ph$log2FC.BTX[!is.na(ph$log2FC.BTX)] < 0)),
                               nrow = 2, byrow = TRUE,
                               dimnames = list(Candidate = c("In", "Out"), Background = c("In", "Out")))
    
  }
  names(cont_tables) <- npgroup_count_BTX$netphorest_group2
  cont_tables[1]
  
  p_value <- numeric()
  for(i in seq_along(cont_tables)){
    
    p_value[i] <- fisher.test(cont_tables[[i]])$p.value
    
  }
  
  cont_tables[which(p_value < 0.01)]
  p_value[which(p_value < 0.01)]
  p.adj.updown <- p.adjust(p_value, "BH")
  
  deregulated_updown <- names(cont_tables[which(p.adj.updown < 0.01)])

#### prepare for plotting ####
  
  # order netphorest groups
  npgroup_count_BTX <- npgroup_count_BTX[order(N_count)]
  npgroup_count_BTX <- npgroup_count_BTX[!is.na(netphorest_group2)]
  
  # mark kinase groups that show significant deregulation of up or down-phosphorylated sites
  npgroup_count_BTX[, np_marked := netphorest_group2]
  npgroup_count_BTX[npgroup_count_BTX$netphorest_group2 %in% deregulated_updown, np_marked := paste0("* ", netphorest_group2)]
  
  # order by high level netphorest groups
  npgroup_count_high_BTX  <- ph_BTX[, list(N_count = .N,
                                         N_up    = sum(Regulation == "Up"),
                                         N_down  = sum(Regulation == "Down")), by = c("netphorest_group_high", "netphorest_group2")]
  npgroup_count_high_BTX <- npgroup_count_high_BTX[order(-N_count), ]
  unique(npgroup_count_high_BTX$netphorest_group_high)
  temp <- split(npgroup_count_high_BTX, npgroup_count_high_BTX$netphorest_group_high, sorted = FALSE)[unique(npgroup_count_high_BTX$netphorest_group_high)]
  npgroup_count_high_BTX <- rbindlist(temp)
  npgroup_count_high_BTX <- npgroup_count_high_BTX[nrow(npgroup_count_high_BTX) : 1]
  npgroup_count_high_BTX[, np_marked := netphorest_group2]
  npgroup_count_high_BTX[npgroup_count_high_BTX$netphorest_group2 %in% deregulated_updown, np_marked := paste0("* ", netphorest_group2)]

  # sort netphorest groups based on the sorting of the higher kinase groups
  ph_BTX[, netphorest_group2 := factor(netphorest_group2, levels = unique(npgroup_count_high_BTX$netphorest_group2))]

##### plot CaEGTA #####
  
  max_y_limit <- 1.1*max(c(npgroup_count_Ca$N_count, npgroup_count_BTX$N_count))
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
  
  pdf("plots\\NetphorestGroups_Count_CaEGTA.pdf", width = 13, heigh = 8)
  
  g <- ggplot(ph_Ca[!(is.na(netphorest_group) | netphorest_group == "")], aes(x = netphorest_group2))
  g <- g + geom_bar(aes(fill = Regulation))
  g <- g + scale_fill_manual(breaks = c("Up", "Down"), values = c("lightblue", "orange"))
  g <- g + scale_y_continuous(expand = c(0, 1), limits = c(-x_side, max_y_limit))
  g <- g + scale_x_discrete(labels = npgroup_count_high_Ca$np_marked)
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

##### plot Mock/BTX #####
  
  temp_BTX <- data.table(label = c("STE",
                               "CAMK group",
                               "CK",
                               "AGC group", 
                               "CMGC group"),
                     xmin = c(5, 8, 14, 16, 27),
                     xmax = c(7, 13, 15, 26, 33),
                     ymin = rep(max(npgroup_count_BTX$N_count) -10, 5),
                     ymax = rep(max(npgroup_count_BTX$N_count) -10, 5))
  temp_BTX <- merge(temp_BTX, temp_Ca[, c("label", "col")], by = "label", all.x = TRUE)
  
  
  pdf("plots\\NetphorestGroups_Count_BTX.pdf", width = 13, heigh = 8)
  
  g <- ggplot(ph_BTX[!(is.na(netphorest_group) | netphorest_group == "")], aes(x = netphorest_group2))
  g <- g + geom_bar(aes(fill = Regulation))
  g <- g + scale_fill_manual(breaks = c("Up", "Down"), values = c("lightblue", "orange"))
  g <- g + scale_y_continuous(expand = c(0, 1), limits = c(-x_side, max_y_limit))
  g <- g + scale_x_discrete(labels = npgroup_count_high_BTX$np_marked)
  g <- g + coord_flip()
  g <- g + labs(title = "Netphorest Groups")
  g <- g + xlab("") + ylab("Number of regulated phosphorylation sites")
  g <- g + theme_bw()
  g <- g + theme(axis.text.y = element_text(face = "bold"),
                 plot.margin = unit(c(1, 1, 1, 0), "cm")
                  )
  g <- g + annotate("rect",
                    xmin = temp_BTX$xmin-0.4,
                    xmax = temp_BTX$xmax+0.4,
                    ymin = -x_side,
                    ymax = -1,
                    fill = temp_BTX$col)
  
  g <- g + annotate("rect",
                    xmin = temp_BTX$xmin-0.4,
                    xmax = temp_BTX$xmax+0.4,
                    ymin = max_y_limit - 10,
                    ymax = max_y_limit - 5,  
                    fill = scales::alpha(temp_BTX$col, 0.6))
  
  
  g <- g + annotate("text", x = (temp_BTX$xmin + temp_BTX$xmax)/2,
                    y = max_y_limit - 10,
                    label = temp_BTX$label,
                    size = 3,
                    fontface = "bold",
                    angle = -90,
                    hjust = 0.5,
                    vjust = -0.5)
  
  print(g)
  
  dev.off()
  
  })

### reduced graphs

# Ca

ph_Ca <- ph_Ca[ph_Ca$netphorest_group2 %in% npgroup_count_Ca$netphorest_group2[npgroup_count_Ca$N_count > 10]]
npgroup_count_high_Ca <- npgroup_count_high_Ca[npgroup_count_high_Ca$netphorest_group2 %in% npgroup_count_Ca$netphorest_group2[npgroup_count_Ca$N_count > 10]]

max_y_limit <- 1.1*max(c(npgroup_count_Ca$N_count, npgroup_count_BTX$N_count))
x_side <- round(max_y_limit/nrow(npgroup_count_Ca), 0) / 2


temp_Ca <- data.table(label = c("STE",
                                "CK",
                                "AGC",
                                "CAMK", 
                                "CMGC"),
                      xmin = c(1, 2, 3, 4, 6),
                      xmax = c(1, 2, 3, 5, 9),
                      ymin = rep(max(npgroup_count_Ca$N_count) -10, 5),
                      ymax = rep(max(npgroup_count_Ca$N_count) -10, 5))
temp_Ca[, col := c("#71c0a8", "#ae6161", "#a7a4cc", "#ec79b4", "#9ec473")]

pdf("plots\\NetphorestGroups_Count_CaEGTA_red.pdf", width = 10, heigh = 7)

g <- ggplot(ph_Ca, aes(x = netphorest_group2))
g <- g + geom_bar(aes(fill = Regulation))
g <- g + scale_fill_manual(breaks = c("Up", "Down"), values = c("lightblue", "orange"))
g <- g + scale_y_continuous(expand = c(0, 1), limits = c(-x_side, max_y_limit))
g <- g + scale_x_discrete(labels = npgroup_count_high_Ca$netphorest_group2)
g <- g + coord_flip()
g <- g + ggtitle("Netphorest Groups")
g <- g + xlab("") + ylab("# regulated phosphorylation sites")
g <- g + theme_light()
g <- g + theme(axis.text.y = element_text(face = "bold",
                                          size = 18),
               axis.text.x = element_text(size = 18),
               axis.title.x = element_text(size = 20),
               plot.title = element_text(size = 20),
               legend.key.size = unit(1.5, "cm"),
               legend.text = element_text(size = 16),
               legend.title = element_text(size = 18),
               plot.margin = unit(c(1, 1, 1, 3.4), "cm")
               )

g <- g + annotate("rect",
                  xmin = temp_Ca$xmin-0.45,
                  xmax = temp_Ca$xmax+0.45,
                  ymin = -x_side,
                  ymax = -1,
                  fill = temp_Ca$col)

g <- g + annotate("rect",
                  xmin = temp_Ca$xmin-0.45,
                  xmax = temp_Ca$xmax+0.45,
                  ymin = max_y_limit - 0.03*max_y_limit,
                  ymax = max_y_limit - 0.01*max_y_limit,  
                  fill = scales::alpha(temp_Ca$col, 0.6))

g <- g + guides(color = guide_legend("Color"))

print(g)

dev.off()

# BTX

ph_BTX <- ph_BTX[ph_BTX$netphorest_group2 %in% npgroup_count_BTX$netphorest_group2[npgroup_count_BTX$N_count > 10]]
npgroup_count_high_BTX <- npgroup_count_high_BTX[npgroup_count_high_BTX$netphorest_group2 %in% npgroup_count_BTX$netphorest_group2[npgroup_count_BTX$N_count > 10]]

temp_BTX <- data.table(label = c("STE",
                                 "CAMK",
                                 "CK",
                                 "AGC", 
                                 "CMGC"),
                      xmin = c(1, 2, 4, 5, 7),
                      xmax = c(1, 3, 4, 6, 11),
                      ymin = rep(max(npgroup_count_BTX$N_count) -10, 5),
                      ymax = rep(max(npgroup_count_BTX$N_count) -10, 5))
temp_BTX <- merge(temp_BTX, temp_Ca[, c("label", "col")], by = "label", all.x = TRUE)

pdf("plots\\NetphorestGroups_Count_BTX_red.pdf", width = 10, heigh = 7)

g <- ggplot(ph_BTX[ph_BTX$netphorest_group2 != ""], aes(x = netphorest_group2))
g <- g + geom_bar(aes(fill = Regulation))
g <- g + scale_fill_manual(breaks = c("Up", "Down"), values = c("lightblue", "orange"))
g <- g + scale_y_continuous(expand = c(0, 1), limits = c(-x_side, max_y_limit))
g <- g + scale_x_discrete(labels = npgroup_count_high_BTX$netphorest_group2)
g <- g + coord_flip()
g <- g + ggtitle("Netphorest Groups")
g <- g + xlab("") + ylab("# regulated phosphorylation sites")
g <- g + theme_light()
g <- g + theme(axis.text.y = element_text(face = "bold",
                                          size = 18),
               axis.text.x = element_text(size = 18),
               axis.title.x = element_text(size = 20),
               plot.title = element_text(size = 20),
               legend.key.size = unit(1.5, "cm"),
               legend.text = element_text(size = 16),
               legend.title = element_text(size = 18),
               plot.margin = unit(c(1, 1, 1, 3.4), "cm")
               )

g <- g + annotate("rect",
                  xmin = temp_BTX$xmin-0.45,
                  xmax = temp_BTX$xmax+0.45,
                  ymin = -x_side,
                  ymax = -1,
                  fill = temp_BTX$col)

g <- g + annotate("rect",
                  xmin = temp_BTX$xmin-0.45,
                  xmax = temp_BTX$xmax+0.45,
                  ymin = max_y_limit - 0.03*max_y_limit,
                  ymax = max_y_limit - 0.01*max_y_limit,  
                  fill = scales::alpha(temp_BTX$col, 0.6))

print(g)

dev.off()