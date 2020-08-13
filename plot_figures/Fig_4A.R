local({
  
  library(data.table)
  library(ggplot2)
  
  # custom function
  lm_eqn <- function(m){
    
    if(coef(m)[[2]] >= 0) a <- paste0(" + ", format(abs(coef(m)[[2]]), digits = 3, nsmall = 3))
    if(coef(m)[[2]] < 0)  a <- paste0(" - ", format(abs(coef(m)[[2]]), digits = 3, nsmall = 3)) 
    
    eq <- substitute(italic(y) == b %.% italic(x)* italic(a)*","~~italic(r)^2~"="~r2, 
                     list(b  = format(coef(m)[[1]], digits = 3, nsmall = 3),
                          a  = a,
                          r2 = format(summary(m)$r.squared,   digits = 3, nsmall = 3)))
    as.character(as.expression(eq))                 
  }  
  
  ph_cand <- fread("Figures\\Fig_4A\\Calcium_vs_Cycling_all.txt")
  
  ph_cand[, Regulation_group := factor(Regulation_group, levels = c("Ca-dependent", "Ca-compensating", "Ca-enhancing", "Cycling-dependent", "not-regulated"))]
  
  lm_mod <- lm(data= ph_cand, log2FC.CaEGTA_BoNT ~ log2FC.BoNT)
  pdf("Figures\\Fig_4A\\Fig_4A.pdf", width = 10 , height = 7)
  
  g <- ggplot(ph_cand, aes(x = log2FC.BoNT, y = log2FC.CaEGTA_BoNT, color = Regulation_group))
  g <- g + geom_point(alpha = 0.6)
  g <- g + scale_color_manual(values = c("#ffa555", "#66c837", "#ff7369", "#30a8ff", "grey"))
  g <- g + geom_abline(slope = lm_mod$coefficients[2], intercept = lm_mod$coefficients[1], color = "steelblue", linetype = "dashed")
  g <- g + xlab("log2 (Mock / BoNT)") + ylab("log2 (Ca / EGTA) - log2 (Mock / BoNT)")
  g <- g + annotate("text", y = 2.25, x = 1.2, color = "steelblue", label = lm_eqn(lm_mod), parse = TRUE, size = 7)
  g <- g + guides(color = guide_legend(override.aes = list(alpha = 1, size = 4), title = "Regulation Group"))
  g <- g + stat_smooth(color = "darkgrey", method = "lm", linetype = "dotted", se = FALSE)
  g <- g + theme_bw()
  g <- g + theme(axis.text = element_text(size = 18),
                 axis.title = element_text(size = 20),
                 legend.text = element_text(size = 20),
                 legend.title = element_text(size = 20),
                 legend.key.size = unit(1.5, "cm"))
  
  print(g)
  
  dev.off()
  
  })