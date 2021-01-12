local({
  
  
  library(data.table)
  library(stringr)
  library(ggplot2)
  
  fc_cut <- log2(1.2)
  
  ph <- fread("Figures\\Fig_4A\\Calcium_vs_Cycling_all.txt")
  ph <- ph[Regulation_group != "not-affected"]
  
  n_sites <- rbindlist(list(
    
    ph[log2FC.CaEGTA > log2(1.2)      & log2FC.BoNT < log2(1/1.2),    list(id = 1, N = .N)],
    ph[log2FC.CaEGTA > log2(1.2)      & abs(log2FC.BoNT) < log2(1.2), list(id = 2, N = .N)],
    ph[log2FC.CaEGTA > log2(1.2)      & log2FC.BoNT > log2(1.2),      list(id = 3, N = .N)],
    ph[abs(log2FC.CaEGTA) < log2(1.2) & log2FC.BoNT < log2(1/1.2),    list(id = 4, N = .N)],
    ph[abs(log2FC.CaEGTA) < log2(1.2) & abs(log2FC.BoNT) < log2(1.2), list(id = 5, N = .N)],
    ph[abs(log2FC.CaEGTA) < log2(1.2) & log2FC.BoNT > log2(1.2),      list(id = 6, N = .N)],
    ph[log2FC.CaEGTA < log2(1/1.2)    & log2FC.BoNT < log2(1/1.2),    list(id = 7, N = .N)],
    ph[log2FC.CaEGTA < log2(1/1.2)    & abs(log2FC.BoNT) < log2(1.2), list(id = 8, N = .N)],
    ph[log2FC.CaEGTA < log2(1/1.2)    & log2FC.BoNT > log2(1.2),      list(id = 9, N = .N)])
    
  )
  n_sites[, label := paste0("n = ", N)]
  n_sites[id %in% c(1, 4, 7), x := -2]
  n_sites[id %in% c(2, 5, 8), x := -0.2]
  n_sites[id %in% c(3, 6, 9), x := 2.5]
  
  n_sites[id %in% c(1, 2, 3), y := 3]
  n_sites[id %in% c(4, 5, 6), y := 0]
  n_sites[id %in% c(7, 8, 9), y := -2]
  

  
  g <- ggplot(ph, aes(x = log2FC.BoNT, y = log2FC.CaEGTA, color = Regulation_group))
  #g <- g + facet_wrap(~netphorest_group)
  g <- g + scale_color_manual(values = c("#fcb533", "#51baf4", "grey"),
                              labels = c(expression(paste("primary Ca"^"2+", "-dependent")),
                                         "SV-cycling-dependent",
                                         "not-affected"))
  g <- g + coord_equal(x = c(-2.2, 3.2), y = c(-2.2, 3.2), expand = FALSE)
  g <- g + geom_point(alpha = 0.4, size = 1)
  g <- g + geom_hline(yintercept = log2(c(1.2, 1/1.2)), color = "darkred", linetype = "dashed")
  g <- g + geom_vline(xintercept = log2(c(1.2, 1/1.2)), color = "darkred", linetype = "dashed")
  g <- g + theme_bw()
  g <- g + geom_text(data = n_sites,
                     aes(x = x, y = y, label = label),
                     color = scales::alpha("black", 0.8), hjust = 0)
  g <- g + xlab(expression(paste("log"[2], " (Mock / BoNT)")))
  g <- g + ylab(expression(paste("log"[2], " (Ca / EGTA)")))
  #g <- g + guides(color = guide_legend(title = "Regulation Group", override.aes = list(alpha = 1), label.hjust = 0))
  g <- g + guides(color = FALSE)
  
  pdf("Figures//Fig_4A//Fig4A_CaEGTA_vs_BoNT_sign.pdf", width = 5, height = 5)
  print(g)
  dev.off()
  
})