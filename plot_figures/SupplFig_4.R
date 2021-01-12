local({
  
  library(data.table)
  library(ggplot2)
  
  sites <- fread("Figures\\SupplFig_4\\site_occupancies.txt")
  
  g <- ggplot(sites,
              aes(x = log2FC.BoNT, y = log2FC.CaEGTA, label = label,
                  size = occupancy_mock, color = Regulation_group_resolved))
  g <- g + coord_equal(x = c(-2.2, 3.2), y = c(-2.2, 3.2), expand = FALSE)
  g <- g + scale_color_manual(values = c("#fcb533", "#51baf4"),
                              labels = c(expression(paste("primary Ca"^"2+", "-dependent")),
                                         "SV-cycling-dependent"))
  g <- g + geom_hline(yintercept = log2(c(1.2, 1/1.2)), color = "darkred", linetype = "dashed")
  g <- g + geom_vline(xintercept = log2(c(1.2, 1/1.2)), color = "darkred", linetype = "dashed")
  g <- g + geom_point(alpha = 0.6)
  g <- g + geom_text_repel(size = 3.5,
                           color = scales::alpha("black", 0.6),
                           segment.color = scales::alpha("darkgrey", 0.8))
  g <- g + guides(color = guide_legend(title = "Regulation Group", override.aes = list(alpha = 0.6, size = 5), label.hjust = 0),
                  size = guide_legend(title = "Occupancy in Mock"))
  g <- g + xlab(expression(paste("log"[2], " (Mock / BoNT)")))
  g <- g + ylab(expression(paste("log"[2], " (Ca / EGTA)")))
  g <- g + theme_bw()
  
  pdf("Figures\\SupplFig_4\\SupplFig_4.pdf", width = 11)
  print(g)
  dev.off()
  
  })