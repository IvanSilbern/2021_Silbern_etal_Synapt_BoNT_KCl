
# DO:
# plot number of regulated sites per manually annotated key-words

# INPUT:
# "temp\\Protein_classification.txt"

# OUTPUT:
# "plots\\SupplFig_ManualAnnot.pdf"
# "Figures\\SupplFig_5\\sites_keyword_count.txt"

local({
  
  if(!dir.exists("Figures\\SupplFig_5")) dir.create("Figures\\Fig_6", recursive = TRUE)
  
  library(data.table)
  library(stringr)
  library(ggplot2)
  
  df <- fread("temp\\Protein_classification.txt")
  names(df)
  
  count_cols <- c("Ca_dependent", "Cycling_dependent")
  
  reg_colors <- c("#fcb533",
                  "#51baf4"
                  )
  
  groups <- merge(df[, list(Group_single = unlist(str_split(Keyword_manual, ";"))), by = "Gene.name"],
                  df[, c("Gene.name", ..count_cols)],
                  by = "Gene.name")
  #groups[, Group_single := gsub(";", "", Group_single)]
  groups
  
  groups <- groups[, lapply(.SD, sum), by = "Group_single", .SDcols = count_cols]
  groups[, Total := Ca_dependent + Cycling_dependent]
  groups <- groups[order(-Total)]
  groups[Group_single == "", Group_single := "Other"]
  
  
  groups_long <- melt(groups[Total > 20, -c("Total")], measure.vars = count_cols, variable.name = "Regulation", value.name = "Count")
  groups_long[, Group_single := factor(Group_single, levels = c(groups$Group_single[groups$Group_single != "Other"], "Other"))]
  groups_long[, Regulation := factor(Regulation, levels = count_cols)]
  regulation <- groups_long$Regulation
  levels(regulation) <- c("primary Ca-dependent", "SV-cycling-dependent")
  groups_long[, Regulation := regulation]
  
  g <- ggplot(groups_long, aes(x = Group_single, y = Count, fill = Regulation))
  g <- g + geom_col(position = "stack", alpha = 0.75)
  g <- g + theme_classic()
  g <- g + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  g <- g + scale_fill_manual(values = reg_colors)
  g <- g + scale_y_continuous(expand = expand_scale(mult = c(0, .1)))
  g <- g + ylab("Number of regulated sites") + xlab("Keyword")
  g <- g + guides(fill = guide_legend(title = "Regulation group"))
  print(g)
  
  pdf("plots\\SupplFig_ManualAnnot.pdf", width = 12, height = 7)
  print(g)
  dev.off()
  
  # provide source data
  
  fwrite(groups_long, "Figures\\SupplFig_5\\sites_keyword_count.txt", sep = "\t")

})