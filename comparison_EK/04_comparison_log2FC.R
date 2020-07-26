
# Do:
# compare log2FC of phosphorylation sites
# match based on +/- 7 aa sequence windows
# log2FC of different multiplicity-states
# for the same site are summed.

# INPUT:
# "comparison_EK\\PhSiteContrasts_EK.tsv"
# "temp\\PhPeptIntensities_slim.tsv"

# OUTPUT:
# "comparison_EK//Scatterplot_comparison_sign_76mM_10s.pdf"
# "Figures\\Fig_2D\\Log2FC_comparison_EK_IS.txt"
# "Figures\\Fig_2D\\Log2FC_comparison_EK_IS_selected".txt


local({
  
  if(!dir.exists("Figures\\Fig_2D")) dir.create("Figures\\Fig_2D", recursive = TRUE)
    
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
  
  # Load data
  ph_r <- fread("comparison_EK\\PhSiteContrasts_EK.tsv")
  ph_s <- fread("temp\\PhPeptIntensities_slim.tsv")
  
  # Use sum to combine different multiplicities; keep the minimum q value in each case
  temp <- ph_s[, list(log2FC.CaEGTA = sum(log2FC.CaEGTA),
                      q.val.CaEGTA = min(q.val.CaEGTA),
                      log2FC.BTX = sum(log2FC.BTX),
                      q.val.BTX = min(q.val.BTX)),
               by = "id"]

  # remove duplicated values due to different multiplicities
  ph_s <- ph_s[!duplicated(ph_s$id)]
  ph_s <- merge(ph_s[, c("id", "Accession", "Gene.name", "Amino.acid", "Position",
                         "Protein.description", "Accession.noIso",
                         "Sequence.window_7", "Sequence.window_15",
                         "Sequence.window_7_noSpace")],
                 temp, by = "id")
  
  ph <- merge(ph_r[, -c("Sequence.window_15", "Sequence.window_7_noSpace")],
              ph_s[, -c("Sequence.window_15", "Sequence.window_7_noSpace")],
              by = "Sequence.window_7", all = TRUE, suffixes = c(".EK", ".IS"))
  
  ph[, Site_id := paste0(Gene.name.IS, "-", Amino.acid, Position.IS)]


  lm_mod <- lm(data = ph[which(Stimulation == "80" & q.val.CaEGTA < 0.01)], T1.norm ~ log2FC.CaEGTA)
  summary(lm_mod)

  ph_sub <- ph[which(Stimulation == "80" & q.val.CaEGTA < 0.01)]

  selected_sites <- c("Bsn-S235",
                      "Bsn-S241",
                      "Bsn-S245",
                      "Syn1-S62",
                      "Dnm1-S774",
                      "Cadps-S89",
                      "Syn1-S67",
                      "Pclo-S591",
                      "Bsn-S821",
                      "Syn1-S603",
                      "Dpsl2-S536",
                      "Mapt-S661",
                      "Dpysl5-S532")
  ph_sub2 <- ph_sub[Site_id %in% selected_sites]

  pdf("comparison_EK\\Scatterplot_comparison_sign_76mM_10s.pdf", width = 7, height = 7)
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
  
  # provide figure source
  fwrite(ph_sub, "Figures\\Fig_2D\\Log2FC_comparison_EK_IS.txt", sep = "\t")
  fwrite(ph_sub2, "Figures\\Fig_2D\\Log2FC_comparison_EK_IS_selected.txt", sep = "\t")


})
