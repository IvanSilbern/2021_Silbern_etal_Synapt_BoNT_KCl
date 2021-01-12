local({
  
  library(data.table)
  library(ggplot2)
  
  path_long <- fread("Figures\\Fig_4D\\reactome_count.txt")
  df_tests  <- fread("Figures\\Fig_4D\\reactome_pvalue.txt")

  # order np groups
  df_tests[, Pathway.name := factor(Pathway.name, levels = unique(path_long$Pathway.name))]
  path_long[, Pathway.name := factor(Pathway.name, levels = unique(path_long$Pathway.name))]
  
  g <- ggplot(path_long, aes(y = Pathway.name, x = Regulation.group))
  g <- g + geom_tile(aes(fill = percentage), color = "white", size = 1)
  g <- g + coord_equal()
  #g <- g + scale_fill_gradient(high = "#fcb533", low = "#ebf7e9", na.value = "white", limits = c(0, 10), oob = scales::squish)
  g <- g + scale_fill_gradient(high = "#eeaaff", low = "#ebf7e9", na.value = "white", limits = c(0, 10), oob = scales::squish)
  g <- g + theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 11),
                 axis.text.y = element_text(size = 16),
  )
  g <- g + xlab("") + ylab("")
  g <- g + guides(fill = guide_colorbar(title = "% Sites"))
  
  pdf("Figures\\Fig_4D\\Fig_4D_Reactome_count.pdf", width = 12, height = 9)
  print(g)
  dev.off()
  
  g <- ggplot(df_tests, aes(comparison, Pathway.name))
  g <- g + geom_tile(aes(fill = p.adj), color = "white", size = 1)
  g <- g + coord_equal()
  g <- g + scale_fill_gradient(low = "#aad400", high = "#ebf7e9", na.value = "white", limits = c(0, 0.2), oob = scales::squish)
  #g <- g + scale_fill_gradient(low = "#51baf4", high = "#ebf7e9", na.value = "white", limits = c(0, 0.2), oob = scales::squish)
  g <- g + scale_x_discrete(expand = c(0, 0), position = "top", breaks = levels(df_tests$comparison))
  g <- g + scale_y_discrete(expand = c(0, 0), limits = rev(levels(df_tests$netphorest_group2)))
  g <- g + theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold", size = 11),
                 axis.text.y = element_text(size = 16),
  )
  g <- g + xlab("") + ylab("")
  g <- g + guides(fill = guide_colorbar(title = "adj. p-Value"))
  
  pdf("Figures\\Fig_4D\\Fig_4D_Reactome_pvalue.pdf", width = 12, height = 9)
  print(g)
  dev.off()
  })