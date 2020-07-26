local({
  
  library(data.table)
  
  # minimal percentage of the sites belonging to the same regulation group
  # to consider the protein being mostly regulated the same way
  min_per <- 60
  
  prot_groups_count <- fread("Figures\\Fig_4B\\Gene_RegulationGroups.txt")
  
  n_ca_dep  <- prot_groups_count[Ca_dependent_per > min_per, .N]
  n_ca_enh  <- prot_groups_count[Ca_enhancing_per > min_per, .N]
  n_ca_comp <- prot_groups_count[Ca_compensating_per > min_per, .N]
  n_cycling <- prot_groups_count[Cycling_dependent_per > min_per, .N]
  n_ca_comp_cycling <- prot_groups_count[Ca_compensating_per + Cycling_dependent_per > min_per, .N]
  n_ca_dep_enh      <- prot_groups_count[Ca_dependent_per + Ca_enhancing_per > min_per, .N]
  other <- prot_groups_count[, .N] - n_ca_comp_cycling - n_ca_dep_enh
  
  pt <- data.table(type = c("Ca_dependent > 60%",
                            "Ca_dep + Ca_enh > 60%",
                            "Ca_compensating > 60%",
                            "Cycling_dependent > 60%",
                            "Ca_comp + Cycling_dep > 60%",
                            "Other"),
                   N = c(n_ca_dep, 
                         n_ca_dep_enh - n_ca_dep,
                         n_ca_comp, n_cycling,
                         n_ca_comp_cycling - n_ca_comp - n_cycling,
                         other))
  pt[, Percent := 100*N/nrow(prot_groups_count)]
  pt[, type := factor(type, pt$type)]
  
  pt[, type2 := c("mostly Ca-dep.",
                  "Ca-enh. or Ca-dep.",
                  "mostly Ca-comp.",
                  "mostly SV-cycling-dep.",
                  "Ca-comp. or SV-cycling-dep.",
                  "Mixed")]
  pt[, type2 := factor(type2, levels = pt$type2)]
  
  colors <- c(# "ff7469",
    "#ffa655",
    "#ff7469",
    "#66c837",
    "#31a9ff",
    "cyan",
    "darkgrey")
  
  pdf("Figures\\Fig_4B\\Fig_4B.pdf", width = 8)
  g <- ggplot(pt, aes(x = "", y = Percent, fill = type2))
  g <- g + geom_bar(stat = "identity", width = 0.5)
  g <- g + annotate("segment", x = 0.5, xend = 1.45, y = 100-sum(pt$Percent[1:2]), yend = 100-sum(pt$Percent[1:2]), color = "white", size = 2)
  g <- g + annotate("segment", x = 0.5, xend = 1.45, y = 0, yend = 0, color = "white", size = 2)
  g <- g + annotate("segment", x = 0.5, xend = 1.45, y = sum(pt$Percent[6]), yend = sum(pt$Percent[6]), color = "white", size = 2)
  g <- g + scale_fill_manual(values = scales::alpha(colors, 0.6))
  g <- g + coord_polar("y", start=0)
  g <- g + guides(fill = guide_legend(""))
  g <- g + theme_void()
  print(g)
  dev.off()
  
  })