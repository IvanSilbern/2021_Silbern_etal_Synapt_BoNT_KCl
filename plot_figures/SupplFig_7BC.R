local({
  
  library(igraph)
  
  load("Figures\\SupplFig_7\\PP1_graph.Rdata")
  load("Figures\\SupplFig_7\\PP1_slices.Rdata")
  
  pdf("Figures\\SupplFig_7\\SupplFig_7B.pdf", width = 5, height = 5)
  plot(sg,
       vertex.label = V(sg)$Label,
       vertex.label.cex = 0.6,
       vertex.color = V(sg)$col,
       vertex.frame.color = V(sg)$frame_color,
       vertex.size = V(sg)$size,
       vertex.label = "",
       vertex.shape = V(sg)$shape,
       vertex.pie = slices,
       vertex.pie.color = list(c(scales::alpha("#fcb533", 0.8), scales::alpha("#51baf4", 0.8))))
  dev.off()
  
  load("Figures\\SupplFig_7\\CN_graph.Rdata")
  load("Figures\\SupplFig_7\\CN_slices.Rdata")
  
  pdf("Figures\\SupplFig_7\\SupplFig_7C.pdf", width = 5, height = 5)
  plot(sg,
       vertex.label = V(sg)$Label,
       vertex.label.cex = 0.6,
       vertex.color = V(sg)$col,
       vertex.frame.color = V(sg)$frame_color,
       vertex.size = V(sg)$size,
       vertex.label = "",
       vertex.shape = V(sg)$shape,
       vertex.pie = slices,
       vertex.pie.color = list(c(scales::alpha("#fcb533", 0.8), scales::alpha("#51baf4", 0.8))))
  dev.off()
  
})