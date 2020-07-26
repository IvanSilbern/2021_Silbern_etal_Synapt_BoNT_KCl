local({
  
  library(ggrepel)
  library(data.table)
  library(ggplot2)
  library(stringr)
  
  lm_eqn <- function(m){
    eq <- substitute(italic(y) == b %.% italic(x)* + italic(a)*","~~italic(r)^2~"="~r2, 
                     list(b  = format(coef(m)[[1]], digits = 3, nsmall = 3),
                          a  = format(coef(m)[[2]], digits = 3, nsmall = 3),
                          r2 = format(summary(m)$r.squared,   digits = 3, nsmall = 3)))
    as.character(as.expression(eq))                 
  }
  
  ph_sub <- fread("Figures\\Fig_2D\\Log2FC_comparison_EK_IS.txt")
  ph_sub2 <- fread("Figures\\Fig_2D\\Log2FC_comparison_EK_IS_selected.txt")
  
  lm_mod <- lm(data = ph_sub, T1.norm ~ log2FC.CaEGTA)
  
  pdf("Figures\\Fig_2D\\Fig_2D.pdf", width = 7, height = 7)
  g <- ggplot(ph_sub, aes(x = log2FC.CaEGTA, y = T1.norm, label = Site_id))
  g <- g + geom_point(alpha = 0.6, color = "steelblue", size = 3)
  g <- g + geom_point(data = ph_sub2, color = "red", size = 3)
  g <- g + stat_smooth(method = "lm", se = FALSE, color = scales::alpha("orange", 0.8))
  g <- g + xlab("log2 Ca / EGTA, 2 min stimulation")
  g <- g + ylab("log2 KCl / Mock-Stimulated,\n10 s stimulation (Engholm-Keller et al)")
  #g <- g + ggtitle("Comparision log2FC to Engholm-Keller et al")
  g <- g + annotate("text", y = 2.7, x = -2.5, color = "steelblue", label = lm_eqn(lm_mod), parse = TRUE, size = 5)
  g <- g + geom_text_repel(data = ph_sub2[which(ph_sub2$log2FC.CaEGTA < -3)],
                           segment.color = "darkgrey",
                           segment.alpha = 0.5,
                           xlim = c(-2.5, -5),
                           ylim = c(-1, -5),
                           size = 4,
                           nudge_x = -0.4,
                           nudge_y = -0.8)
  g <- g + geom_text_repel(data = ph_sub2[which(ph_sub2$log2FC.CaEGTA < -2.2 & ph_sub2$log2FC.CaEGTA > -3)],
                           segment.color = "darkgrey",
                           segment.alpha = 0.5,
                           xlim = c(-2.5, -5),
                           ylim = c(-1, 2.5),
                           size = 4,
                           nudge_x = -0.2,
                           nudge_y = 0.8)
  g <- g + geom_text_repel(data = ph_sub2[which(ph_sub2$log2FC.CaEGTA > 2.2)],
                           segment.color = "darkgrey",
                           segment.alpha = 0.5,
                           xlim = c(2, 5),
                           #ylim = c(-2, -3),
                           size = 4,
                           nudge_x = 0.2)
  g <- g + geom_text_repel(data = ph_sub2[which(ph_sub2$Site_id == "Syn1-S603")],
                           segment.color = "darkgrey",
                           segment.alpha = 0.5,
                           #xlim = c(2, 5),
                           ylim = c(2.2, 2.5),
                           size = 4,
                           nudge_x = 0.8)
  g <- g + geom_text_repel(data = ph_sub2[which(ph_sub2$T1.norm < -1.5)],
                           segment.color = "darkgrey",
                           segment.alpha = 0.5,
                           #xlim = c(2, 5),
                           ylim = c(-1.8, -4.5),
                           size = 4,
                           nudge_x = 0.2)
  g <- g + theme_bw()
  g <- g + theme(axis.text  = element_text(size = 14),
                 axis.title = element_text(size = 18))
  g <- g + scale_x_continuous(limits = c(-5, 5))
  g <- g + scale_y_continuous(limits = c(-5, 5))
  g <- g + coord_equal()
  print(g)
  dev.off()
  
  })