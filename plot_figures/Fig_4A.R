local({
  
  library(data.table)
  library(ggplot2)
  
  ph_cand <- fread("Figures\\Fig_4A\\Calcium_vs_Cycling_all.txt")
  ph_cand[, Regulation_group := factor(Regulation_group, levels = c("primary Ca-dependent", "SV-cycling-dependent", "not-affected"))]
  
  coord_min <- min(c(ph_cand$log2FC.CaEGTA, ph_cand$log2FC.BoNT), na.rm = TRUE)
  coord_max <- max(c(ph_cand$log2FC.CaEGTA, ph_cand$log2FC.BoNT), na.rm = TRUE)
  
  pdf("Figures\\Fig_4A\\Fig_4A.pdf", width = 10 , height = 7)
  g <- ggplot(ph_cand, aes(x = log2FC.BoNT, y = log2FC.CaEGTA, color = Regulation_group))
  g <- g + geom_point(alpha = 0.4)
  g <- g + coord_equal(xlim = c(coord_min, coord_max), ylim = c(coord_min, coord_max))
  g <- g + scale_color_manual(values = c("#fcb533", "#51baf4ff", "grey"))
  g <- g + xlab("log2FC (Mock / BoNT)") + ylab("log2FC (Ca / EGTA)")
  g <- g + geom_hline(yintercept = log2(c(1.2, 1/1.2)), color = "darkred", linetype = "dashed")
  g <- g + geom_vline(xintercept = log2(c(1.2, 1/1.2)), color = "darkred", linetype = "dashed")
  g <- g + guides(color = guide_legend(override.aes = list(alpha = 1, size = 4), title = "Regulation Group"))
  g <- g + theme_bw()
  g <- g + theme(axis.text = element_text(size = 18),
                 axis.title = element_text(size = 20),
                 legend.text = element_text(size = 20),
                 legend.title = element_text(size = 20),
                 legend.key.size = unit(1.5, "cm"))
  print(g)
  dev.off()
  
  })