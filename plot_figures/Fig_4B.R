local({
  
  library(data.table)
  
  # minimal percentage of the sites belonging to the same regulation group
  # to consider the protein being mostly regulated the same way
  min_per <- 75
  
  prot_groups_count <- fread("Figures\\Fig_4B\\Gene_RegulationGroups.txt")
  
  n_ca_dep  <- prot_groups_count[Ca_dependent_per > min_per, .N]
  n_cycling <- prot_groups_count[Cycling_dependent_per > min_per, .N]
  other <- prot_groups_count[, .N] - n_ca_dep - n_cycling
  
  pt <- data.table(type = c("mostly primary Ca-dependent",
                            "mostly SV-cycling-dependent",
                            "Mixed"),
                   N = c(n_ca_dep, 
                         n_cycling,
                         other))
  pt[, Percent := 100*N/nrow(prot_groups_count)]
  pt[, type := factor(type, pt$type)]
  
  pt[, type := factor(type, levels = pt$type)]
  
  colors <- c("#fcb533",
              "#51baf4",
              "#c79a8d")
  
  pdf("Figures\\Fig_4B\\Fig_4B.pdf", width = 8)
  
  g <- ggplot(pt, aes(x = "", y = Percent, fill = type))
  g <- g + geom_bar(stat = "identity", width = 0.5)
  g <- g + annotate("segment", x = 0.5, xend = 1.25, y = 100-sum(pt$Percent[1:2]), yend = 100-sum(pt$Percent[1:2]), color = "white", size = 2)
  g <- g + annotate("segment", x = 0.5, xend = 1.25, y = 0, yend = 0, color = "white", size = 2)
  g <- g + annotate("segment", x = 0.5, xend = 1.25, y = sum(pt$Percent[2:3]), yend = sum(pt$Percent[2:3]), color = "white", size = 2)
  g <- g + scale_fill_manual(values = scales::alpha(colors, 0.6))
  g <- g + coord_polar("y", start=0)
  g <- g + guides(fill = guide_legend(""))
  g <- g + theme_void()
  g <- g + theme(legend.text = element_text(size = 20),
                 legend.title = element_text(size = 20),
                 legend.key.size = unit(1.5, "cm"))
  print(g)
  
  dev.off()
  
  })