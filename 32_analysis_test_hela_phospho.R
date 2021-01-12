# Do:
# analyse phosphorylation sites intensities
# test HeLa nuclear extract sample
# treated with C.botulinum cell culture supernatants

# INPUT:
# "search_results\\Test_Hela_Ph\\combined\\txt\\Phospho (STY)Sites.txt"

# OUTPUT:
# "plots\\Test_HeLA_PhSites_volcano_qVal_log2FC_BoNT_UT.pdf"
# "plots\\Test_HeLA_PhSites_volcano_qVal_log2FC_BoNT_Mock.pdf"
# "plots\\Test_HeLA_PhSites_volcano_qVal_log2FC_Mock_UT.pdf"
# "Figures\\SupplFig_20\\ProteinGroups_HeLa_PhSites_Contrasts.tsv"


local({
        
        if(!dir.exists("Figures\\SupplFig_20")) dir.create("Figures\\SupplFig_20", recursive = TRUE)
                
        library(data.table)
        library(stringr)
        library(ggplot2)
        library(limma)
        library(qvalue)
        
        ph <- fread("search_results\\Test_Hela_Ph\\combined\\txt\\Phospho (STY)Sites.txt", check.names = TRUE)
        dim(ph)
        
        ph <- ph[Reverse != "+"]
        ph <- ph[Potential.contaminant != "+"]
        dim(ph)
        
        names(ph)
        
        int_cols <- grep("^Intensity\\..+___[1-3]$", names(ph), value = TRUE)
        temp <- ph[, c("id", ..int_cols)]
        
        # convert to numeric format
        temp[, c(int_cols) := lapply(.SD, as.numeric), .SDcols = int_cols]
        
        # convert df to a long format
        temp <- melt(temp, measure.vars = int_cols, variable.factor = FALSE)
        
        # extract experiment name from 'variable'
        temp[, Experiment := str_match(variable, "^Intensity\\.([^_]+)")[, 2]]
        temp[, Treatment  := str_match(Experiment, "^([^\\d]+)")[, 2]]
        temp[, Replicate  := str_match(Experiment, "(\\d+)$")[, 2]]
        temp[, Multiplicity  := str_match(variable, "___(\\d)$")[, 2]]
        
        # valid values; imputation
        n.zero <- temp[, list(n.zero = sum(value == 0)), by = c("id", "Multiplicity")]
        n.zero[, to.delete := n.zero > 3]
        
        min.value <- temp[value != 0, list(min.value = min(value, na.rm = TRUE)), 
                          by = c("Experiment")]
        
        norm.min <- temp[value != 0, list(SD = sd(log2(value), na.rm = TRUE),
                                          percentile = quantile(log2(value), probs = 0.05, na.rm = T)),
                         by = c("Experiment")]
        
        temp <- merge(temp, n.zero,    all.x = TRUE)
        temp <- merge(temp, min.value, by = c("Experiment"), all.x = TRUE)
        temp <- merge(temp, norm.min, by = c("Experiment"), all.x = TRUE)
        
        temp[, value2 := value]
        temp[to.delete == TRUE, value2 := NA_real_]
        
        # use probability distribution for missing value imputation
        set.seed(100)
        temp[value2 == 0, value2 := 2^abs(rnorm(n = .N, mean = percentile, sd = 2*SD))]
        
        temp <- temp[!is.na(value2)]
        
        int_cols <- unique(temp$Experiment)
        temp <- dcast(temp, id + n.zero + Multiplicity ~ Experiment, value.var = "value2")
        dim(temp)
        
        pdf("plots\\Test_Hela_PhSites_Log2_raw_intensities.pdf", width = 9)
        boxplot(log2(temp[, ..int_cols]))
        dev.off()
        
        # log transform
        temp[, c(int_cols):= lapply(.SD, log2), .SDcols = int_cols]
        
        # Tukey median polish
        dat     <- temp[, ..int_cols]
        dat.med <- medpolish(dat, maxiter = 3, na.rm = TRUE)
        dat.res <- dat.med$residuals
        
        colnames(dat.res) <- paste0(colnames(dat.res), ".norm")
        temp <- cbind(temp, dat.res)
        
        # columns to use for testing
        int_cols_norm <- paste0(int_cols, ".norm")
        
        # define samples
        samples <- factor(c("BoNT", "BoNT", "BoNT", "Mock", "Mock", "Mock", "UT", "UT", "UT"), levels = c("BoNT", "Mock", "UT"))
        design <- model.matrix(~0 + samples)
        colnames(design) <- levels(samples)
        
        # testing
        fit <- lmFit(temp[, ..int_cols_norm], design = design)
        
        # specify contrasts
        cont.matrix <- makeContrasts(BoNT_UT = BoNT - UT,
                                     BoNT_Mock = BoNT - Mock,
                                     Mock_UT = Mock - UT,
                                     levels = design)
        
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit.eb <- eBayes(fit2)
        
        log2FC <- as.data.table(fit.eb$coef, stringsAsFactors = FALSE)
        log2FC <- melt(log2FC, measure.vars = names(log2FC), value.name = "log2FC", variable.name = "Contrast")
        
        
        t.mod <- as.data.table(fit.eb$t, stringsAsFactors = FALSE)
        t.mod <- melt(t.mod, measure.vars = names(t.mod), value.name = "t.mod", variable.name = "Contrast")
        
        
        p.mod <- as.data.table(fit.eb$p.value, stringsAsFactors = FALSE)
        p.mod <- melt(p.mod, measure.vars = names(p.mod), value.name = "p.mod", variable.name = "Contrast")
        
        q.obj <- qvalue(p.mod$p.mod)
        q.val <- q.obj$qvalues
        
        
        p.adj.BH <- p.adjust(p.mod$p.mod, "BH")
        
        res.eb.mult <- data.table(log2FC, t.mod = t.mod$t.mod, p.mod = p.mod$p.mod, q.val, p.adj.BH, stringsAsFactors = FALSE)
        dim(res.eb.mult)
        
        res.eb.mult[, id := rep(temp$id, ncol(cont.matrix))]
        res.eb.mult[, Multiplicity := rep(temp$Multiplicity, ncol(cont.matrix))]
        
        take.columns <- c("id",
                          "Proteins",
                          "Leading.proteins",
                          "Fasta.headers",
                          "Localization.prob",
                          "Amino.acid",
                          "Position",
                          "Number.of.Phospho..STY.",
                          "Sequence.window",
                          "Modification.window"
                          )
        
        temp <- merge(ph[, ..take.columns], temp, by = "id", all.y = TRUE)
        temp <- merge(temp, res.eb.mult, by = c("id", "Multiplicity"), all = TRUE)
        dim(temp)
        
        candidates <- which(temp$q.val < 0.01 &
                            abs(temp$log2FC) > log2(1.5))
        
        temp[, Candidate := FALSE]
        temp[candidates, Candidate := TRUE]
        
        temp[, Enriched := FALSE]
        temp[which(temp$log2FC > 0), Enriched := TRUE]
        
        pdf("plots\\Test_HeLA_PhSites_volcano_qVal_log2FC_BoNT_UT.pdf")
        plot(y = -log10(temp$q.val),
             x = temp$log2FC,
             type = "n",
             xlab = "log2-Ratio",
             ylab = "-log10 qValue")
        points(y = -log10(temp$q.val[!temp$Candidate & temp$Contrast == "BoNT_UT"]),
               x = temp$log2FC[!temp$Candidate & temp$Contrast == "BoNT_UT"],
               pch = 21,
               bg = "lightgrey")
        points(y = -log10(temp$q.val[temp$Candidate & temp$Contrast == "BoNT_UT"]),
               x = temp$log2FC[temp$Candidate & temp$Contrast == "BoNT_UT"],
               bg = "orange",
               pch = 21)
        abline(v = log2(1/1.5), col = "blue", lty = 3, lwd = 2)
        abline(v = log2(1.5), col = "blue", lty = 3, lwd = 2)
        abline(h = -log10(0.01), col = "blue", lty = 3, lwd = 2)
        dev.off()
        
        pdf("plots\\Test_HeLA_PhSites_volcano_qVal_log2FC_BoNT_Mock.pdf")
        plot(y = -log10(temp$q.val),
             x = temp$log2FC,
             type = "n",
             xlab = "log2-Ratio",
             ylab = "-log10 qValue")
        points(y = -log10(temp$q.val[!temp$Candidate & temp$Contrast == "BoNT_Mock"]),
               x = temp$log2FC[!temp$Candidate & temp$Contrast == "BoNT_Mock"],
               pch = 21,
               bg = "lightgrey")
        points(y = -log10(temp$q.val[temp$Candidate & temp$Contrast == "BoNT_Mock"]),
               x = temp$log2FC[temp$Candidate & temp$Contrast == "BoNT_Mock"],
               bg = "orange",
               pch = 21)
        abline(v = log2(1/1.5), col = "blue", lty = 3, lwd = 2)
        abline(v = log2(1.5), col = "blue", lty = 3, lwd = 2)
        abline(h = -log10(0.01), col = "blue", lty = 3, lwd = 2)
        
        dev.off()
        
        pdf("plots\\Test_HeLA_PhSites_volcano_qVal_log2FC_Mock_UT.pdf")
        plot(y = -log10(temp$q.val),
             x = temp$log2FC,
             type = "n",
             xlab = "log2-Ratio",
             ylab = "-log10 qValue")
        points(y = -log10(temp$q.val[!temp$Candidate & temp$Contrast == "Mock_UT"]),
               x = temp$log2FC[!temp$Candidate & temp$Contrast == "Mock_UT"],
               pch = 21,
               bg = "lightgrey")
        points(y = -log10(temp$q.val[temp$Candidate & temp$Contrast == "Mock_UT"]),
               x = temp$log2FC[temp$Candidate & temp$Contrast == "Mock_UT"],
               bg = "orange",
               pch = 21)
        abline(v = log2(1/1.5), col = "blue", lty = 3, lwd = 2)
        abline(v = log2(1.5), col = "blue", lty = 3, lwd = 2)
        abline(h = -log10(0.01), col = "blue", lty = 3, lwd = 2)
        
        dev.off()
        
        temp[, Accession := unlist(lapply(str_split(Leading.proteins, " "), "[", 1))]
        temp[, Accession := str_match(Accession, "\\|([^\\s])+\\|")[, 2]]
        temp[, Accession.noIso := gsub("-[0-9]+$", "", Accession)]
        temp[, Gene.name := str_match(Fasta.headers, "\\sGN=([^\\s]+)")[, 2]]
        temp[, Protein.name := str_match(Fasta.headers, "\\s(.*?)\\sOS=")[,2]]
        
        
        fwrite(temp, "Figures\\SupplFig_20\\ProteinGroups_HeLa_PhSites_Contrasts.txt", sep = "\t")

})
