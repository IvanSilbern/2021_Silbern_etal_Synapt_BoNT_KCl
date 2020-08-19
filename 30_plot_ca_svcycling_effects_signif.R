
local({

  library(data.table)
  library(stringr)
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
  
  df <- fread("temp\\PhPeptIntensities4.tsv")
  names(df)
  
  unique(df$Regulation_group_resolved)
  df[, Regulation_group := factor(Regulation_group, levels = c("Ca-dependent", "Ca-compensating", "Ca-enhancing", "Cycling-dependent", "not-regulated"))]
  reg_group <- df$Regulation_group
  levels(reg_group) <-c("primary Ca-dependent", "Ca-opposing", "Ca-synergistic", "Ca-undefined", "not-regulated")
  df[, Regulation_group := reg_group]
  df$Regulation_group_resolved
  
  lm_mod <- lm(data= df, log2FC.CaEGTA_BoNT ~ log2FC.BoNT)
  
  g <- ggplot(df[Regulation_group_resolved != "not-regulated"], aes(x = log2FC.BoNT, y = log2FC.CaEGTA_BoNT, color = Regulation_group))
  g <- g + geom_point(alpha = 0.6)
  g <- g + scale_color_manual(values = c("#ffa555", "#66c837", "#ff7369", "#30a8ff", "grey"))
  g <- g + geom_abline(slope = lm_mod$coefficients[2], intercept = lm_mod$coefficients[1], color = "steelblue", linetype = "dashed")
  g <- g + xlab("log2 (Mock / BoNT)") + ylab("log2 (Ca / EGTA) - log2 (Mock / BoNT)")
  g <- g + annotate("text", y = 2.25, x = 1.2, color = "steelblue", label = lm_eqn(lm_mod), parse = TRUE, size = 7)
  g <- g + guides(color = guide_legend(override.aes = list(alpha = 1, size = 4), title = "Regulation Group"))
  g <- g + theme_bw()
  g <- g + theme(axis.text = element_text(size = 18),
                 axis.title = element_text(size = 20),
                 legend.text = element_text(size = 20),
                 legend.title = element_text(size = 20),
                 legend.key.size = unit(1.5, "cm"))
  print(g)
  
  pdf("plots\\Ca-effect_SV-effect_sign.pdf", width = 10 , height = 7)
  print(g)
  dev.off()
  
  g <- ggplot(df[Regulation_group_resolved != "not-regulated"], aes(x = log2FC.BoNT, y = log2FC.CaEGTA, color = Regulation_group))
  g <- g + geom_point(alpha = 0.6)
  g <- g + scale_color_manual(values = c("#ffa555", "#66c837", "#ff7369", "#30a8ff", "grey"))
  g <- g + xlab("log2 (Mock / BoNT)") + ylab("log2 (Ca / EGTA) - log2 (Mock / BoNT)")
  g <- g + guides(color = guide_legend(override.aes = list(alpha = 1, size = 4), title = "Regulation Group"))
  g <- g + geom_hline(yintercept = log2(c(1.2, 1/1.2)), color = "grey", linetype = "dashed")
  g <- g + geom_vline(xintercept = log2(c(1.2, 1/1.2)), color = "grey", linetype = "dashed")
  g <- g + theme_bw()
  g <- g + theme(axis.text = element_text(size = 18),
                 axis.title = element_text(size = 20),
                 legend.text = element_text(size = 20),
                 legend.title = element_text(size = 20),
                 legend.key.size = unit(1.5, "cm"))
  
  pdf("plots\\CaEGTA_MockBoNT_sign.pdf", width = 10 , height = 7)
  print(g)
  dev.off()
  
})