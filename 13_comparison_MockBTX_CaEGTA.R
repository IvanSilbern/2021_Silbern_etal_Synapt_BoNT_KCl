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

  # Compare BTX and CaEGTA candidates
  ph <- fread("temp\\PhPeptIntensities_slim.tsv", check.names = TRUE)

  # calculate log2FC difference CaEGTA and BTX
  ph[, log2FC.CaEGTA_BTX := log2FC.CaEGTA - log2FC.BTX]
  
  # add entry id
  names(ph)[names(ph) == "id"] <- "phosphosite_table_id"
  ph[, id := 1:.N]

  # add site annotation
  ph[, Ph_site := paste0(Position, Amino.acid)]
  ph[, Site_id := paste0(Gene.name, "-", Ph_site)]
  ph[, Site_id2 := paste0(Gene.name, "-", Ph_site, "_", Multiplicity)]

  # first compare all sites found in both experiments
  ph_cand <- ph[!is.na(ph$q.val.CaEGTA) & !is.na(ph$q.val.BTX)]

  # candidate genes in CaEGTA and BTX experiment
  gn_caegta <- unique(ph$Gene.name[ph$Candidate.CaEGTA])
  gn_btx    <- unique(ph$Gene.name[ph$Candidate.BTX])

  # intersecting genes
  gn <- intersect(gn_caegta, gn_btx)
  gn <- gn[!is.na(gn)]
  gn <- gn[!grepl('"', gn)]

  writeLines(gn, "temp\\Candidate_GeneNames_both.tsv")

  # candidate uniprot accessions (without isoform number)
  acc_caegta <- unique(ph$Accession.noIso[ph$Candidate.CaEGTA])
  acc_btx    <- unique(ph$Accession.noIso[ph$Candidate.BTX])
  
  # intersecting accessions
  acc <- intersect(acc_caegta, acc_btx)
  acc <- acc[!is.na(acc)]
  
  writeLines(acc_caegta, "temp\\Candidate_Accession_CaEGTA.tsv")
  writeLines(acc_btx,    "temp\\Candidate_Accession_BTX.tsv")
  writeLines(acc,        "temp\\Candidate_Accession_Both.tsv")

  # phosphoevent is significant in:
  ph_cand[, Significant_in := "None"]
  ph_cand[ph_cand$Candidate.CaEGTA & !ph_cand$Candidate.BTX, Significant_in := "Ca vs EGTA"]
  ph_cand[!ph_cand$Candidate.CaEGTA & ph_cand$Candidate.BTX, Significant_in := "No Toxin vs Toxin"]
  ph_cand[ph_cand$Candidate.CaEGTA & ph_cand$Candidate.BTX,  Significant_in := "Both"]
  
  ph_cand[, Significant_in := factor(Significant_in, levels = c("Ca vs EGTA", "No Toxin vs Toxin", "Both", "None"))]

  ##### Group based on log2FC #####
  
  # first compare all sites found in both experiments
  ph_cand <- ph[!is.na(ph$q.val.CaEGTA) & !is.na(ph$q.val.BTX)]
  
  # distribute sites in different regulation groups
  fc_cut <- log2(1.2)
  # not regulated
  ph_cand[(log2FC.CaEGTA_BTX <  fc_cut & log2FC.CaEGTA_BTX > -fc_cut) & (log2FC.BTX <  fc_cut & log2FC.BTX > -fc_cut), Regulation_group := "not-regulated"]
  # cycling-dependent
  ph_cand[(log2FC.CaEGTA_BTX <  fc_cut & log2FC.CaEGTA_BTX > -fc_cut) & (log2FC.BTX >  fc_cut), Regulation_group := "Cycling-dependent"]
  ph_cand[(log2FC.CaEGTA_BTX <  fc_cut & log2FC.CaEGTA_BTX > -fc_cut) & (log2FC.BTX < -fc_cut), Regulation_group := "Cycling-dependent"]
  # Ca-compensating
  ph_cand[(log2FC.CaEGTA_BTX < -fc_cut) & (log2FC.BTX >  fc_cut), Regulation_group := "Ca-compensating"]
  ph_cand[(log2FC.CaEGTA_BTX >  fc_cut) & (log2FC.BTX < -fc_cut), Regulation_group := "Ca-compensating"]
  # Ca-enhancing
  ph_cand[(log2FC.CaEGTA_BTX >  fc_cut) & (log2FC.BTX >  fc_cut), Regulation_group := "Ca-enhancing"]
  ph_cand[(log2FC.CaEGTA_BTX < -fc_cut) & (log2FC.BTX < -fc_cut), Regulation_group := "Ca-enhancing"]
  # Ca-dependent
  ph_cand[(log2FC.CaEGTA_BTX >  fc_cut) & (log2FC.BTX <  fc_cut & log2FC.BTX > -fc_cut), Regulation_group := "Ca-dependent"]
  ph_cand[(log2FC.CaEGTA_BTX < -fc_cut) & (log2FC.BTX <  fc_cut & log2FC.BTX > -fc_cut), Regulation_group := "Ca-dependent"]
  
  # direction
  ph_cand[log2FC.CaEGTA_BTX >  fc_cut, Regulation_direction := "pos"]
  ph_cand[log2FC.CaEGTA_BTX < -fc_cut, Regulation_direction := "neg"]
  ph_cand[Regulation_group == "Cycling-dependent" & log2FC.BTX >  fc_cut, Regulation_direction := "pos"]
  ph_cand[Regulation_group == "Cycling-dependent" & log2FC.BTX < -fc_cut, Regulation_direction := "neg"]

  # plot
  
  
  ph_cand[, Regulation_group := factor(Regulation_group, levels = c("Ca-dependent", "Ca-compensating", "Ca-enhancing", "Cycling-dependent", "not-regulated"))]
  
  lm_mod <- lm(data= ph_cand, log2FC.CaEGTA_BTX ~ log2FC.BTX)
  pdf("plots\\Ca_vs_BTX_all.pdf", width = 10 , height = 7)
  g <- ggplot(ph_cand, aes(x = log2FC.BTX, y = log2FC.CaEGTA_BTX, color = Regulation_group))
  g <- g + geom_point(alpha = 0.4)
  g <- g + scale_color_manual(values = c("#ffa555", "#66c837", "#ff7369", "#30a8ff", "grey"))
  g <- g + geom_abline(slope = lm_mod$coefficients[2], intercept = lm_mod$coefficients[1], color = "steelblue", linetype = "dashed")
  g <- g + xlab("log2 (Mock / BTX)") + ylab("log2 (Ca / EGTA) - log2 (Mock / BTX)")
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

  # provide figure source
  fwrite(ph_cand[, c("phosphosite_table_id", "Gene.name", "Accession", "Amino.acid",
                     "Position", "Multiplicity", "q.val.CaEGTA", "log2FC.CaEGTA", "Candidate.CaEGTA",
                     "q.val.BTX", "log2FC.BTX", "Candidate.BTX", "log2FC.CaEGTA_BTX",
                     "Regulation_group")],
         "Figures\\Fig_4A\\Calcium_vs_Cycling_all.txt", sep = "\t")
  

  # Cycling-dependent group should show q.val.BTX < 0.01
  # Ca-dependent group should show q.val.CaEGTA < 0.01
  # Ca-compensating or Ca-enhancing group should show q.val < 0.01 in either CaEGTA or BTX
  
  dim(ph_cand)
  ph_cand <- ph_cand[ph_cand$Regulation_group != "not-regulated"]
  
  list_cand <- list(Cycling_dependent = ph_cand[ph_cand$Regulation_group == "Cycling-dependent" & q.val.BTX < 0.01],
                    Ca_dependent = ph_cand[ph_cand$Regulation_group == "Ca-dependent" & q.val.CaEGTA < 0.01],
                    Ca_compensating = ph_cand[ph_cand$Regulation_group == "Ca-compensating" & (q.val.CaEGTA < 0.01 | q.val.BTX < 0.01)],
                    Ca_enhancing = ph_cand[ph_cand$Regulation_group == "Ca-enhancing" & (q.val.CaEGTA < 0.01 | q.val.BTX < 0.01)])
  lapply(list_cand, nrow)
  
  for(i in seq_along(list_cand)){
    
    fwrite(list_cand[[i]], paste0("temp\\PhPept_Intensities_", names(list_cand)[i], ".tsv"), sep = "\t")
    
  }

  ph_cand <- rbindlist(list_cand)
  dim(ph_cand)

  
  # Ca-dependent if signifcantly regulated in CaEGTA but not found in BTX experiment
  ph[which(ph$Candidate.CaEGTA & is.na(ph$Candidate.BTX)), Regulation_group := "Ca-dependent"]
  ph[ph$Regulation_group == "Ca-dependent" & ph$log2FC.CaEGTA > 0, Regulation_direction := "pos"]
  ph[ph$Regulation_group == "Ca-dependent" & ph$log2FC.CaEGTA < 0, Regulation_direction := "neg"]

  # Cycling-dependent if signifcantly regulated in BTX experiment but not found in CaEGTA experiment
  ph[which(ph$Candidate.BTX & is.na(ph$Candidate.CaEGTA)), Regulation_group := "Cycling-dependent"]
  ph[ph$Regulation_group == "Cycling-dependent" & ph$log2FC.BTX > 0, Regulation_direction := "pos"]
  ph[ph$Regulation_group == "Cycling-dependent" & ph$log2FC.BTX < 0, Regulation_direction := "neg"]

  
  ph_cand <- rbind(ph_cand, ph[which(ph$Candidate.CaEGTA & is.na(ph$Candidate.BTX))])
  ph_cand <- rbind(ph_cand, ph[which(ph$Candidate.BTX & is.na(ph$Candidate.CaEGTA))])

  pdf("plots\\Ca_vs_BTX_sign.pdf", width = 10 , height = 8)
  g <- ggplot(ph_cand, aes(x = log2FC.BTX, y = log2FC.CaEGTA_BTX, color = Regulation_group))
  g <- g + geom_point(alpha = 0.4)
  g <- g + scale_color_manual(values = c("gold", "green3", "salmon", "steelblue", "grey"))
  g <- g + geom_abline(slope = lm_mod$coefficients[2], intercept = lm_mod$coefficients[1], color = "steelblue", linetype = "dashed")
  g <- g + xlab("log2FC.BTX") + ylab("log2FC.CaEGTA - log2FC.BTX")
  g <- g + annotate("text", y = 2.25, x = 2, color = "steelblue", label = lm_eqn(lm_mod), parse = TRUE)
  print(g)
  dev.off()

  table(ph_cand$Regulation_group)

  fwrite(ph_cand, "temp\\PhIntensities_cand.tsv", sep = "\t")

  dim(ph_cand)
  length(unique(ph_cand$Site_id))
  
  unique(ph_cand$Regulation_group)

  prot_groups_count <- ph_cand[, list(Ca_dependent = sum(Regulation_group == "Ca-dependent"),
                                      Ca_compensating = sum(Regulation_group == "Ca-compensating"),
                                      Ca_enhancing = sum(Regulation_group == "Ca-enhancing"),
                                      Cycling_dependent = sum(Regulation_group == "Cycling-dependent")), by = Gene.name]
  prot_groups_per <- prot_groups_count[, lapply(.SD, function(x) 100*x/sum(apply(.SD, 1, sum))),
                                       by = Gene.name, .SDcols = c("Ca_dependent",
                                                                   "Ca_compensating",
                                                                   "Ca_enhancing",
                                                                   "Cycling_dependent")]
  fwrite(prot_groups_count, "temp\\Gene_RegulationGroups_count.tsv", sep = "\t")
  fwrite(prot_groups_per, "temp\\Gene_RegulationGroups_per.tsv", sep = "\t")

})
