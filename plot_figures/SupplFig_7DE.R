local({
  
  library(data.table)
  library(ggplot2)
  
  colors <- c("#fcb533",
              "#51baf4")
  
  per_pp1 <- fread("Figures\\SupplFig_7\\sites_pp1_per.txt")
  count_pp1 <- fread("Figures\\SupplFig_7\\sites_pp1_count.txt")
  
  per_cn <- fread("Figures\\SupplFig_7\\sites_cn_per.txt")
  count_cn <- fread("Figures\\SupplFig_7\\sites_cn_count.txt")
  
  g <- ggplot(per_pp1, aes(x = "", y = V1, fill = Regulation_group_resolved))
  g <- g + geom_bar(stat = "identity")
  g <- g + annotate("segment", x = 0.55, xend = 1.45,
                    y = per_pp1[Regulation_group_resolved == "SV-cycling-dependent"]$V1,
                    yend = per_pp1[Regulation_group_resolved == "SV-cycling-dependent"]$V1,
                    color = "white", size = 1.25)
  g <- g + annotate("segment", x = 0.55, xend = 1.45,
                    y = 0, yend = 0,
                    color = "white", size = 1.25)
  g <- g + annotate("text", label = paste0(count_pp1[Regulation_group_resolved == "SV-cycling-dependent"]$N,
                                           "\n(",
                                           round(100*per_pp1[Regulation_group_resolved == "SV-cycling-dependent"]$V1, 1),
                                           "%)"),
                    x = 1, y = 0.25, color = "white", size = 6)
  
  g <- g + annotate("text", label = paste0(count_pp1[Regulation_group_resolved == "primary Ca-dependent"]$N,
                                           "\n(",
                                           round(100*per_pp1[Regulation_group_resolved == "primary Ca-dependent"]$V1, 1),
                                           "%)"),
                    x = 1, y = 0.75, color = "white", size = 6)
  
  g <- g + coord_polar("y")
  g <- g + theme_void()
  g <- g + scale_fill_manual(values = scales::alpha(colors, 0.8),
                             labels = c(expression(paste("primary Ca"^"2+", "-dependent")),
                                        "SV-cycling-dependent"))
  g <- g + guides(fill = guide_legend(title = "Regulation group", label.hjust = 0))
  print(g)
  
  pdf("Figures\\SupplFig_7\\SupplFig_7D_PP1_sites_pie.pdf", width = 6, height = 4)
  print(g)
  dev.off()
  
  g <- ggplot(per_cn, aes(x = "", y = V1, fill = Regulation_group_resolved))
  g <- g + geom_bar(stat = "identity")
  g <- g + annotate("segment", x = 0.55, xend = 1.45,
                    y = per_cn[Regulation_group_resolved == "SV-cycling-dependent"]$V1,
                    yend = per_cn[Regulation_group_resolved == "SV-cycling-dependent"]$V1,
                    color = "white", size = 1.25)
  g <- g + annotate("segment", x = 0.55, xend = 1.45,
                    y = 0, yend = 0,
                    color = "white", size = 1.25)
  g <- g + annotate("text", label = paste0(count_cn[Regulation_group_resolved == "SV-cycling-dependent"]$N,
                                           "\n(",
                                           round(100*per_cn[Regulation_group_resolved == "SV-cycling-dependent"]$V1, 1),
                                           "%)"),
                    x = 1, y = 0.25, color = "white", size = 6)
  
  g <- g + annotate("text", label = paste0(count_cn[Regulation_group_resolved == "primary Ca-dependent"]$N,
                                           "\n(",
                                           round(100*per_cn[Regulation_group_resolved == "primary Ca-dependent"]$V1, 1),
                                           "%)"),
                    x = 1, y = 0.75, color = "white", size = 6)
  
  g <- g + coord_polar("y")
  g <- g + theme_void()
  g <- g + scale_fill_manual(values = scales::alpha(colors, 0.8),
                             labels = c(expression(paste("primary Ca"^"2+", "-dependent")),
                                        "SV-cycling-dependent"))
  g <- g + guides(fill = guide_legend(title = "Regulation group", label.hjust = 0))
  
  pdf("Figures\\SupplFig_7\\SupplFig_7E_CN_sites_pie.pdf", width = 6, height = 4)
  print(g)
  dev.off()
  
})