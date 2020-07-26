local({
  
  library(data.table)
  library(ggplot2)
  
  counts <- fread("Figures\\Fig_3A\\NPh_group_GOBP.txt")
  take_groups <- fread("Figures\\Fig_3A\\netphorest_groups.txt")
  gobp <- fread("Figures\\Fig_3A\\GO_BP.txt")
  
  #take_groups <- take_groups[order(-N)]
  gobp <- gobp[order(Benjamini)]
  
  counts[, np_group := factor(np_group, levels = c(take_groups$netphorest_group2, "Other"))]
  counts[, GO_Term := factor(go_term, levels = gobp$Term)]
  
  pdf("Figures\\Fig_3A\\Fig_3A.pdf", width = 12 , heigh = 10)
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