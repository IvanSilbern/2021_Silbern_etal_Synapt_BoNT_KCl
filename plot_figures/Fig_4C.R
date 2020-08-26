local({
  
  
  library(data.table)
  library(ggplot2)
  
  sites_by_np_reg <- fread("Figures\\Fig_4C\\np_reg_sites_count.txt")
  df_tests        <- fread("Figures\\Fig_4C\\np_reg_sites_pval.txt") 
  
  # order Regulation groups
  sites_by_np_reg[, Regulation_group_resolved := factor(Regulation_group_resolved, levels = c("primary Ca-dependent", "SV-cycling-dependent"))]
  
  # order netphorest groups at the highest level
  np_groups_high <- sites_by_np_reg[, sum(max_N), by = netphorest_group_high]
  np_groups_high <- np_groups_high[order(-V1)]
  
  sites_by_np_reg[, netphorest_group_high := factor(netphorest_group_high, levels = np_groups_high$netphorest_group_high)]
  
  # order Netphorest groups
  np_groups_ordered <- unique(sites_by_np_reg$netphorest_group2[order(sites_by_np_reg$netphorest_group_high, -sites_by_np_reg$max_N)])
  sites_by_np_reg[, netphorest_group2 := factor(netphorest_group2, levels = np_groups_ordered)]
  
  # order np groups
  df_tests[, netphorest_group2 := factor(netphorest_group2, levels = np_groups_ordered)]

  # order "comparisons"
  df_tests[, comparison := gsub("_", " vs\n", comparison)]
  df_tests[, comparison := factor(comparison, levels = c("primary Ca-dependent vs\nSV-cycling-dependent"))]
  
    ##### plots #####  
  
  pdf("Figures\\Fig_4C\\Fig_4C.pdf", width = 8, height = 10)
  
  g <- ggplot(sites_by_np_reg, aes(Regulation_group_resolved, netphorest_group2))
  
  g <- g + geom_tile(aes(fill = N_per), color = "white", size = 1)
  g <- g + coord_equal()
  g <- g + geom_text(aes(label = paste0(round(N_per, 1), "%\n(", N, ")")), color = "grey40")
  g <- g + annotate("text", x = c(1, 2), y = c(0.5, 0.5), label = c("", ""))
  g <- g + annotate("text", x = c(1, 2), y = c(0.95, 0.95),
                    label = c(sites_by_np_reg$Total[sites_by_np_reg$Regulation_group_resolved == "primary Ca-dependent"][1],
                              sites_by_np_reg$Total[sites_by_np_reg$Regulation_group_resolved == "SV-cycling-dependent"][1]),
                    fontface = "bold")
  g <- g + annotate("segment", x = 0.5, xend = 2.5, y = 1.4, yend = 1.4, color = "black", size = 1.5)
  #g <- g + scale_fill_gradient(high = "#fcb533", low = "#ebf7e9", na.value = "white", limits = c(0, 30), oob = scales::squish)
  g <- g + scale_fill_gradient(high = "#eeaaff", low = "#ebf7e9", na.value = "white", limits = c(0, 30), oob = scales::squish)
  g <- g + scale_x_discrete(expand = c(0, 0), position = "top", breaks = levels(sites_by_np_reg$Regulation_group_resolved))
  g <- g + scale_y_discrete(expand = c(0, 0), limits = rev(c(levels(sites_by_np_reg$netphorest_group2), "Total")))
  g <- g + theme(plot.margin = unit(c(3,3,3,3), "cm"),
                 axis.text.x = element_text(angle = 30, hjust = 0, face = "bold", size = 11),
                 axis.text.y = element_text(size = 16),
                 legend.key.size = unit(1, "cm"),
                 legend.text = element_text(size = 14),
                 legend.title = element_text(size = 18),
                 panel.grid = element_blank(),
                 panel.background = element_rect(fill = "white"),
                 axis.ticks.y = element_blank()
  )
  g <- g + xlab("") + ylab("")
  g <- g + guides(fill = guide_colorbar(title = "% Sites"))
  print(g)
  
  dev.off()
  
  pdf("Figures\\Fig_4C\\Fig_4C_2.pdf", width = 8, height = 10)
  
  g <- ggplot(df_tests, aes(comparison, netphorest_group2))
  g <- g + geom_tile(aes(fill = p.adj), color = "white", size = 1)
  g <- g + coord_equal()
  g <- g + geom_text(aes(label = round(p.adj, 3)), color = "grey40")
  #g <- g + scale_fill_gradient(low = "#51baf4", high = "#ebf7e9", na.value = "white", limits = c(0, 0.2), oob = scales::squish)
  g <- g + scale_fill_gradient(low = "#aad400", high = "#ebf7e9", na.value = "white", limits = c(0, 0.2), oob = scales::squish)
  g <- g + scale_x_discrete(expand = c(0, 0), position = "top", breaks = levels(df_tests$comparison))
  g <- g + scale_y_discrete(expand = c(0, 0), limits = rev(levels(df_tests$netphorest_group2)))
  g <- g + theme(plot.margin = unit(c(3,3,3,3), "cm"),
                 axis.text.x = element_text(angle = 30, hjust = 0, face = "bold", size = 10),
                 axis.text.y = element_text(size = 16),
                 legend.key.size = unit(1, "cm"),
                 legend.text = element_text(size = 14),
                 legend.title = element_text(size = 18),
                 panel.grid = element_blank(),
                 panel.background = element_rect(fill = "white"),
                 axis.ticks.y = element_blank()
  )
  g <- g + xlab("") + ylab("")
  g <- g + guides(fill = guide_colorbar(title = "adj. p-Value"))
  print(g)
  
  dev.off()
  
  })