# DO: comparison on MockBoNT and CaEGTA experiments

# INPUT: 
# "temp\\PhPeptIntensities_slim.tsv"

#OUTPUT:
# "temp\\Candidate_GeneNames_both.tsv"
# "temp\\Candidate_Accession_CaEGTA.tsv"
# "temp\\Candidate_Accession_BoNT.tsv"
# "temp\\Candidate_Accession_Both.tsv"
# "plots\\CaEGTA_vs_MockBoNT_all.pdf"
# "temp\\PhIntensities_cand.tsv"
# "temp\\Gene_RegulationGroups_per.tsv"



local({

  if(!dir.exists("Figures\\Fig_4A")) dir.create("Figures\\Fig_4A", recursive = TRUE)
  if(!dir.exists("plots")) dir.create("plots")
  
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  set.seed(500)
  
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

  # Compare BoNT and CaEGTA candidates
  ph <- fread("temp\\PhPeptIntensities_slim.tsv", check.names = TRUE)

  # calculate log2FC difference CaEGTA and BoNT
  ph[, log2FC.CaEGTA_BoNT := log2FC.CaEGTA - log2FC.BoNT]
  
  # add entry id
  names(ph)[names(ph) == "id"] <- "phosphosite_table_id"
  ph[, id := 1:.N]

  # add site annotation
  ph[, Ph_site := paste0(Position, Amino.acid)]
  ph[, Site_id := paste0(Gene.name, "-", Ph_site)]
  ph[, Site_id2 := paste0(Gene.name, "-", Ph_site, "_", Multiplicity)]

  # first compare all sites found in both experiments
  ph_cand <- ph[!is.na(ph$q.val.CaEGTA) & !is.na(ph$q.val.BoNT)]

  # candidate genes in CaEGTA and BoNT experiment
  gn_caegta <- unique(ph$Gene.name[ph$Candidate.CaEGTA])
  gn_bont    <- unique(ph$Gene.name[ph$Candidate.BoNT])

  # intersecting genes
  gn <- intersect(gn_caegta, gn_bont)
  gn <- gn[!is.na(gn)]
  gn <- gn[!grepl('"', gn)]

  writeLines(gn, "temp\\Candidate_GeneNames_both.tsv")

  # candidate uniprot accessions (without isoform number)
  acc_caegta <- unique(ph$Accession.noIso[ph$Candidate.CaEGTA])
  acc_bont    <- unique(ph$Accession.noIso[ph$Candidate.BoNT])
  
  # intersecting accessions
  acc <- intersect(acc_caegta, acc_bont)
  acc <- acc[!is.na(acc)]
  
  writeLines(acc_caegta, "temp\\Candidate_Accession_CaEGTA.tsv")
  writeLines(acc_bont,   "temp\\Candidate_Accession_BoNT.tsv")
  writeLines(acc,        "temp\\Candidate_Accession_Both.tsv")

  # phosphoevent is significant in:
  ph_cand[, Significant_in := "None"]
  ph_cand[ph_cand$Candidate.CaEGTA & !ph_cand$Candidate.BoNT, Significant_in := "Ca vs EGTA"]
  ph_cand[!ph_cand$Candidate.CaEGTA & ph_cand$Candidate.BoNT, Significant_in := "No Toxin vs Toxin"]
  ph_cand[ph_cand$Candidate.CaEGTA & ph_cand$Candidate.BoNT,  Significant_in := "Both"]
  
  ph_cand[, Significant_in := factor(Significant_in, levels = c("Ca vs EGTA", "No Toxin vs Toxin", "Both", "None"))]

  ##### Group based on log2FC #####
  
  # first compare all sites found in both experiments
  ph_cand <- ph[!is.na(ph$q.val.CaEGTA) & !is.na(ph$q.val.BoNT)]
  
  # distribute sites in different regulation groups
  fc_cut <- log2(1.2)
  
  # not regulated
  ph_cand[, Regulation_group := "not-affected"]
  
  # Ca-dependent
  ph_cand[abs(log2FC.CaEGTA) >  fc_cut & abs(log2FC.BoNT) <  fc_cut, Regulation_group := "primary Ca-dependent"]
  
  # SV-cycling-dependent
  ph_cand[abs(log2FC.BoNT) > fc_cut, Regulation_group := "SV-cycling-dependent"]
  
  # plot
  ph_cand[, Regulation_group := factor(Regulation_group, levels = c("primary Ca-dependent", "SV-cycling-dependent", "not-affected"))]
  
  coord_min <- min(c(ph_cand$log2FC.CaEGTA, ph_cand$log2FC.BoNT), na.rm = TRUE)
  coord_max <- max(c(ph_cand$log2FC.CaEGTA, ph_cand$log2FC.BoNT), na.rm = TRUE)
  
  pdf("plots\\CaEGTA_vs_MockBoNT_all.pdf", width = 10 , height = 7)
  g <- ggplot(ph_cand, aes(x = log2FC.BoNT, y = log2FC.CaEGTA, color = Regulation_group))
  g <- g + geom_point(alpha = 0.4)
  g <- g + coord_equal(xlim = c(coord_min, coord_max), ylim = c(coord_min, coord_max))
  g <- g + scale_color_manual(values = c("#fcb533", "#51baf4", "grey"))
  g <- g + xlab("log2 (Mock / BoNT)") + ylab("log2 (Ca / EGTA)")
  g <- g + guides(color = guide_legend(override.aes = list(alpha = 1, size = 4), title = "Regulation Group"))
  g <- g + geom_hline(yintercept = log2(c(1.2, 1/1.2)), color = "darkred", linetype = "dashed")
  g <- g + geom_vline(xintercept = log2(c(1.2, 1/1.2)), color = "darkred", linetype = "dashed")
  g <- g + theme_bw()
  g <- g + theme(axis.text = element_text(size = 18),
                 axis.title = element_text(size = 20),
                 legend.text = element_text(size = 20),
                 legend.title = element_text(size = 20),
                 legend.key.size = unit(1.5, "cm"))
  
  print(g)
  dev.off()

  # provide figure source
  fwrite(ph_cand[, c("phosphosite_table_id", "Gene.name", "Accession", "Amino.acid",
                     "Position", "Multiplicity", "q.val.CaEGTA", "log2FC.CaEGTA", "Candidate.CaEGTA",
                     "q.val.BoNT", "log2FC.BoNT", "Candidate.BoNT", "log2FC.CaEGTA_BoNT",
                     "Regulation_group")],
         "Figures\\Fig_4A\\Calcium_vs_Cycling_all.txt", sep = "\t")
  

  # SV-Cycling-dependent group should show q.val.BoNT < 0.01
  # primary Ca-dependent group should show q.val.CaEGTA < 0.01
  
  ph_cand <- ph_cand[ph_cand$Regulation_group != "not-regulated"]
  
  list_cand <- list(Cycling_dependent = ph_cand[ph_cand$Regulation_group == "SV-cycling-dependent" & (q.val.BoNT < 0.01 | Candidate.CaEGTA)],
                    Ca_dependent = ph_cand[ph_cand$Regulation_group == "primary Ca-dependent" & q.val.CaEGTA < 0.01])
  lapply(list_cand, nrow)
  
  for(i in seq_along(list_cand)){
    
    fwrite(list_cand[[i]], paste0("temp\\PhPept_Intensities_", names(list_cand)[i], ".tsv"), sep = "\t")
    
  }

  ph_cand <- rbindlist(list_cand)
  dim(ph_cand)

  
  # Ca-dependent if signifcantly regulated in CaEGTA but not found in BoNT experiment
  ph[which(ph$Candidate.CaEGTA & is.na(ph$Candidate.BoNT)), Regulation_group := "primary Ca-dependent"]
  
  # Cycling-dependent if signifcantly regulated in BoNT experiment but not found in CaEGTA experiment
  ph[which(ph$Candidate.BoNT & is.na(ph$Candidate.CaEGTA)), Regulation_group := "SV-cycling-dependent"]
  
  ph_cand <- rbind(ph_cand, ph[which(ph$Candidate.CaEGTA & is.na(ph$Candidate.BoNT))])
  ph_cand <- rbind(ph_cand, ph[which(ph$Candidate.BoNT & is.na(ph$Candidate.CaEGTA))])
  dim(ph_cand)
  table(ph_cand$Regulation_group)
  
  pdf("plots\\Ca_vs_BoNT_sign.pdf", width = 10 , height = 8)
  g <- ggplot(ph_cand, aes(x = log2FC.BoNT, y = log2FC.CaEGTA, color = Regulation_group))
  g <- g + geom_point(alpha = 0.4)
  g <- g + coord_equal(xlim = c(coord_min, coord_max), ylim = c(coord_min, coord_max))
  g <- g + scale_color_manual(values = c("#fcb533", "#51baf4ff", "grey"))
  g <- g + xlab("log2FC (Mock / BoNT)") + ylab("log2FC (Ca / EGTA)")
  g <- g + geom_hline(yintercept = log2(c(1.2, 1/1.2)), color = "darkred", linetype = "dashed")
  g <- g + geom_vline(xintercept = log2(c(1.2, 1/1.2)), color = "darkred", linetype = "dashed")
  g <- g + guides(color = guide_legend(override.aes = list(alpha = 1, size = 4), title = "Regulation Group"))
  g <- g + theme_bw()
  g <- g + theme(axis.text = element_text(size = 18),
                 axis.title = element_text(size = 20),
                 legend.text = element_text(size = 20),
                 legend.title = element_text(size = 20),
                 legend.key.size = unit(1.5, "cm"))
  print(g)
  dev.off()

  table(ph_cand$Regulation_group)

  fwrite(ph_cand, "temp\\PhIntensities_cand.tsv", sep = "\t")

  dim(ph_cand)
  length(unique(ph_cand$Site_id))
  
  unique(ph_cand$Regulation_group)

  prot_groups_count <- ph_cand[, list(Ca_dependent = sum(Regulation_group == "primary Ca-dependent"),
                                      Cycling_dependent = sum(Regulation_group == "SV-cycling-dependent")), by = Gene.name]
  prot_groups_per <- prot_groups_count[, lapply(.SD, function(x) 100*x/sum(apply(.SD, 1, sum))),
                                       by = Gene.name, .SDcols = c("Ca_dependent",
                                                                   "Cycling_dependent")]
  fwrite(prot_groups_count, "temp\\Gene_RegulationGroups_count.tsv", sep = "\t")
  fwrite(prot_groups_per, "temp\\Gene_RegulationGroups_per.tsv", sep = "\t")

})
