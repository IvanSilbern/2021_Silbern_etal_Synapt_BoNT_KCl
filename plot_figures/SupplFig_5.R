local({
  
  count_cols <- c("primary Ca-dependent", "SV-cycling-dependent")
  reg_colors <- c("#fcb533",
                  "#51baf4")
  
  groups_long <- fread("Figures\\SupplFig_5\\sites_keyword_count.txt")
  groups_long[, Group_single := factor(Group_single, levels = c(groups$Group_single[groups$Group_single != "Other"], "Other"))]
  groups_long[, Regulation := factor(Regulation, levels = count_cols)]
  
  
  g <- ggplot(groups_long, aes(x = Group_single, y = Count, fill = Regulation))
  g <- g + geom_col(position = "stack", alpha = 0.75)
  g <- g + theme_classic()
  g <- g + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  g <- g + scale_fill_manual(values = reg_colors)
  g <- g + scale_y_continuous(expand = expand_scale(mult = c(0, .1)))
  g <- g + ylab("Number of regulated sites") + xlab("Keyword")
  g <- g + guides(fill = guide_legend(title = "Regulation group"))
  print(g)
  
  pdf("Figures\\SupplFig_6\\SupplFig_5_ManualAnnot.pdf", width = 12, height = 7)
  print(g)
  dev.off()
  
  })