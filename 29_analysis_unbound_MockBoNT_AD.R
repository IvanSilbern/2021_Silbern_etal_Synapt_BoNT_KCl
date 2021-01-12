
# Do:
# Analyse protein group intensities
# remove reverse sequences and potential contaminants
# keep protein groups with 2 or more unique/razor peptides

# INPUT:
# "search_results\\Unbound_MockBoNT_AD\\combined\\txt\\proteinGroups.txt"

# OUTPUT:
# "plots\\MockBoNT_AD_log2_raw_intensities.pdf"
# "plots\\ProteinGroups_MockBoNT_AD_Volcano_pmod_log2FC.pdf"
# "plots\\MockBoNT_AD_selected_proteins.pdf"
# "SupplData\\SupplData02_Proteingroups_Intensities_MockBoNT_AD.txt"
# "Figures\\SupplFig_21ABC\\ProteinGroups_MockBoNT_AD.txt"

local({

  if(!dir.exists("Figures\\SupplFig_21ABC")) dir.create("Figures\\SupplFig_21ABC", recursive = TRUE)
  if(!dir.exists("SupplData")) dir.create("SupplData")
        
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(limma)
  library(qvalue)

  pg <- fread("search_results\\Unbound_MockBoNT_AD\\combined\\txt\\proteinGroups.txt.gz", check.names = TRUE)

  # remove potential contaminants and reversed sequences
  # exclude proteins with < 2 uniqute/razor peptides from quantification
  # 
  pg <- pg[Reverse != "+"]
  pg <- pg[Potential.contaminant != "+"]
  pg <- pg[Peptide.counts..razor.unique. > 1]

  int_cols <- grep("^Reporter\\.intensity\\.corrected\\.[1-6]\\.", names(pg), value = TRUE)
  temp <- pg[, c("id", ..int_cols)]

  # convert to numeric format        
  temp[, c(int_cols) := lapply(.SD, as.numeric), .SDcols = int_cols]

  # convert df to a long format
  temp <- melt(temp, measure.vars = int_cols, variable.factor = FALSE)

  # extract experiment name from 'variable'
  temp[, Experiment := "Prep31"]
  temp[, Channel := str_match(variable, "\\.corrected\\.([1-6])")[, 2]]

  n.zero <- temp[, list(n.zero = sum(value == 0)), by = c("id")]
  n.zero[, to.delete := n.zero > 3]
  sum(n.zero$to.delete)

  min.value <- temp[value != 0, list(min.value = min(value, na.rm = TRUE)), 
                    by = c("Channel")]

  norm.min <- temp[value != 0, list(SD = sd(log2(value), na.rm = TRUE),
                                    percentile = quantile(log2(value), probs = 0.05, na.rm = T)),
                   by = c("Channel")]

  temp <- merge(temp, n.zero, by = "id", all.x = TRUE)
  temp <- merge(temp, min.value, by = c("Channel"), all.x = TRUE)
  temp <- merge(temp, norm.min, by = c("Channel"), all.x = TRUE)

  temp[, value2 := value]
  temp[to.delete == TRUE, value2 := NA_real_]

  # % of missing (zero) values per replicate
  temp[!is.na(value2), 100*sum(value2 == 0)/.N, by = c("Channel")] 

  # use probability distribution for missing value imputation
  set.seed(100)
  temp[value2 == 0, value2 := 2^abs(rnorm(n = .N, mean = percentile, sd = 2*SD))]

  temp <- temp[!is.na(value2)]
  temp <- dcast(temp, id + n.zero ~ variable, value.var = "value2")
  dim(temp)

  int_cols <- c("Mock_01", "BoNT_01", "Mock_02", "BoNT_02", "Mock_03", "BoNT_03")
  names(temp) <- c("id", "n.zero", int_cols)

  pdf("plots\\MockBoNT_AD_log2_raw_intensities.pdf")
  boxplot(log2(temp[, c("Mock_01", "Mock_02", "Mock_03", "BoNT_01", "BoNT_02", "BoNT_03")]))
  dev.off()

  # log transform
  temp[, c(int_cols):= lapply(.SD, log2), .SDcols = int_cols]
  
  # Tukey median polishing
  dat     <- temp[, ..int_cols]
  dat.med <- medpolish(dat, maxiter = 3, na.rm = TRUE)
  dat.res <- dat.med$residuals
 
  colnames(dat.res) <- paste0(colnames(dat.res), ".norm")
  temp <- cbind(temp, dat.res)

  int_cols_norm <- paste0(int_cols, ".norm")


  tr <- c("Mock", "BoNT", "Mock", "BoNT", "Mock", "BoNT")
  tr <- factor(tr, levels = c("BoNT", "Mock"))
  design <- model.matrix(~tr)

  fit <- lmFit(temp[, ..int_cols_norm], design = design)
  fit.eb <- eBayes(fit)

  n      <- dim(temp[, ..int_cols_norm])[1]
  log2FC <- fit.eb$coef[, "trMock"]
  df.0   <- rep(fit.eb$df.prior, n)
  df.r   <- fit.eb$df.residual
  s2.0   <- rep(fit.eb$s2.prior, n)
  s2     <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coef[, "trMock"]/fit.eb$sigma/fit.eb$stdev.unscaled[, "trMock"]
  t.mod <- fit.eb$t[, "trMock"]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, "trMock"]
  p.adj.BH <- p.adjust(p.mod, "BH")
  sum(p.adj.BH < 0.01)

  res.eb.mult <- data.table(log2FC, t.ord, t.mod, p.ord, p.mod, df.r, df.0, s2.0, s2, s2.post, p.adj.BH)
  dim(res.eb.mult)

  temp <- cbind(temp, res.eb.mult)

  take.columns <- c("id",
                    "Protein.IDs",
                    "Majority.protein.IDs",
                    "Fasta.headers",
                    "Number.of.proteins",
                    "Peptides",
                    "Razor...unique.peptides",
                    "Unique.peptides",
                    "Sequence.coverage....",
                    "Mol..weight..kDa.",
                    "Sequence.length",
                    "Q.value", "Score", "Intensity", "MS.MS.count")

  temp <- merge(pg[, ..take.columns], temp, by = "id", all.y = TRUE)
  dim(temp)

  candidates <- temp$p.adj.BH < 0.01 &
                      abs(temp$log2FC) > log2(1.5)
  
  temp[, Candidate := FALSE]
  temp[candidates, Candidate := TRUE]

  temp[, Enriched := FALSE]
  temp[which(temp$log2FC > 0), Enriched := TRUE]

  temp[, Gene.name := str_match(Fasta.headers, "GN=([^\\s]+)")[, 2]]
  temp_sub <- temp[Gene.name %in% c("Vamp2", "Stx1a", "Snap25", "Actb")]
  temp_sub2 <- temp[abs(log2FC) > 0.7]

  pdf("plots\\ProteinGroups_MockBoNT_AD_Volcano_pmod_log2FC.pdf")
  plot(y = -log10(temp$p.mod),
       x = temp$log2FC,
       type = "n",
       xlab = "log2-Ratio",
       ylab = "-log10 modarated p-value",
       xlim = c(-3, 2),
       ylim = c(0, 7))
  points(y = -log10(temp$p.mod[!temp$Candidate]),
         x = temp$log2FC[!temp$Candidate],
         pch = 21,
         bg = "lightgrey")
  points(y = -log10(temp$p.mod[temp$Candidate]),
         x = temp$log2FC[temp$Candidate],
         bg = "orange",
         pch = 21)

  abline(v = log2(1/1.5), col = "blue", lty = 3, lwd = 2)
  abline(v = log2(1.5), col = "blue", lty = 3, lwd = 2)

  dev.off()

  temp_sub <- melt(temp_sub[, c("Gene.name", ..int_cols_norm)], measure.vars = int_cols_norm)
  temp_sub[, variable := gsub("\\.norm$", "", variable)]
  temp_sub[, variable := factor(variable, levels = c("Mock_01", "Mock_02", "Mock_03", "BoNT_01", "BoNT_02", "BoNT_03"))]

  pdf("plots\\MockBoNT_AD_selected_proteins.pdf", width = 6, height = 4)
  g <- ggplot(temp_sub, aes(x = variable, y = value, group = Gene.name, color = Gene.name))
  g <- g + geom_line()
  g <- g + geom_point()
  g <- g + scale_y_continuous(limits = c(-0.6, 0.6))
  g <- g + geom_hline(yintercept = c(-0.3, 0.3), color = "blue", linetype = 2)
  g <- g + theme_bw()
  g <- g + xlab("Replicate") + ylab("log2 normalized intensity") 
  g <- g + guides(color = guide_legend(title = "Protein"))
  print(g)
  dev.off()

  temp[, Accession := unlist(lapply(str_split(Majority.protein.IDs, ";"), "[", 1))]
  temp[, Accession := str_match(Accession, "\\|([^\\s]+)\\|")[, 2]]
  temp[, Accession.noIso := gsub("-[0-9]+$", "", Accession)]

  temp[, Gene.name := str_match(Fasta.headers, "\\sGN=([^\\s]+)")[, 2]]
  temp[, Protein.name := str_match(Fasta.headers, "\\s(.*?)\\sOS=")[,2]]

  take_cols <- c("id",
                 "Accession",
                 "Gene.name",
                 "Protein.name",
                 "Protein.IDs",
                 "Majority.protein.IDs",
                 "Fasta.headers",
                 "Number.of.proteins",
                 "Peptides",
                 "Razor...unique.peptides",
                 "Unique.peptides",
                 "Sequence.coverage....",
                 "Mol..weight..kDa.",
                 "Sequence.length",
                 "Q.value",
                 "Score",
                 "Intensity",
                 "MS.MS.count",
                 "Mock_01",
                 "Mock_02",
                 "Mock_03",
                 "BoNT_01",
                 "BoNT_02",
                 "BoNT_03",
                 "Mock_01.norm",
                 "Mock_02.norm",
                 "Mock_03.norm",
                 "BoNT_01.norm",
                 "BoNT_02.norm",
                 "BoNT_03.norm",
                 "log2FC",
                 "p.mod",
                 "p.adj.BH",
                 "Candidate")
  
  df_sub <- temp[, ..take_cols]
  
  names(df_sub)[names(df_sub) == "Protein.name"] <- "Protein.description"
  names(df_sub)[names(df_sub) == "Razor...unique.peptides"] <- "Razor.and.unique.peptides"
  names(df_sub)[names(df_sub) == "Sequence.coverage...."] <- "Sequence.coverage.percent"
  names(df_sub)[names(df_sub) == "Mol..weight..kDa."] <- "Mol.weight.kDa"
  
  names(df_sub)
  
  fwrite(df_sub, "SupplData\\SupplData02_Proteingroups_Intensities_MockBoNT_AD.txt", sep = "\t")
  
  # provide figure source data
  fwrite(df_sub, "Figures\\SupplFig_21ABC\\ProteinGroups_MockBoNT_AD.txt", sep = "\t")

})