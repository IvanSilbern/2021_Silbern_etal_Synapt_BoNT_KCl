# DO:
# Read String db for rat proteins
# Compute connection score that is based only experimental evidence
# and evidence derived from data bases
# define a minimal score cut-off
# build a graph
# save under temp/stringDB_graph.RData
#
# INPUT:
# "external\\data_Stringdb\\10116.protein.links.full.v10.5.txt.gz" (String data base, downloaded from web page)
# 
# OUTPUT:
# "temp\\stringDB_graph.Rdata"

local({
  
  if(!dir.exists("temp")) dir.create("temp")
  
  library(igraph)
  library(data.table)
  
  links_full <- fread("external\\data_Stringdb\\10116.protein.links.full.v10.5.txt.gz", check.names = TRUE)
  head(links_full)
  
  #Clean Data
  # keep only entries with confidence more than the threshold
  
  # direct and transferred experimental evidence, without prior
  links_full$exp_direct_nopr <- (links_full$experiments/1000 - 0.041) / (1 - 0.041) 
  links_full$exp_transf_nopr <- (links_full$experiments_transferred/1000 - 0.041) / (1 - 0.041) 
  
  # direct and transferred evidence from data bases, without prior
  links_full$db_direct_nopr <- (links_full$database/1000 - 0.041) / (1 - 0.041) 
  links_full$db_transf_nopr <- (links_full$database_transferred/1000 - 0.041) / (1 - 0.041) 
  
  # combined direct and transferred evidence for experiments and data bases
  links_full$exp_comb <- 1 - (1 - links_full$exp_direct_nopr) * (1 - links_full$exp_transf_nopr)
  links_full$db_comb  <- 1 - (1 - links_full$db_direct_nopr) * (1 - links_full$db_transf_nopr)
  
  # combined evidence for experiments and data bases
  links_full$exp_db_comb <- 1 - (1 - links_full$exp_comb) * (1 - links_full$db_comb) 
  
  # add prior
  links_full$exp_db_comb <- links_full$exp_db_comb + 0.041 * (1 - links_full$exp_db_comb)
  links_full$exp_db_comb[links_full$exp_db_comb < 0] <- 0
  links_full$exp_db_comb <- links_full$exp_db_comb * 1000
  
  # head(links_full[, c("combined_score", "exp_db_comb")])
  # sum(links_full$exp_db_comb > 400)
  
  links_full <- links_full[links_full$exp_db_comb > 400, ]
  links_full <- links_full[order(links_full$exp_db_comb, decreasing = TRUE), ]
  # sum(duplicated(links_full[, c("protein1", "protein2")]))
  
  graph <- graph_from_data_frame(links_full, directed = F)
  head(graph)
  
  graph <- simplify(graph, edge.attr.comb = "max")
  E(graph)$weight <- E(graph)$exp_db_comb / 1000
  
  save(graph, file = "temp\\stringDB_graph.Rdata")
  
})