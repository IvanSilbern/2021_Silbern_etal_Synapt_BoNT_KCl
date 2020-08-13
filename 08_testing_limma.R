# DO:
# Imputation, normalization and statistical testing

# INPUT:
# "temp\\Phosphosites_prepared2.tsv"

# OUTPUT:
# Figure source files
# "Figures\\Fig_3BC\\Fig_2B_source.txt" and 
# "Figures\\Fig_3BC\\Fig_2C_source.txt"
# 
# "plots\\Reporter_Int_distributions.pdf"
# "plots\\Pvalue_histogram_CaEGTA.pdf"
# "plots\\Pvalue_histogram_BoNT.pdf"
#
# "temp\\PhPeptIntensities.tsv"
# "temp\\PhPeptIntensities_slim.tsv"
# "temp\\PhPeptIntensities_log2FC_long.tsv"
# "temp\\PhPeptIntensities_int_long.tsv"
# "temp\\Candidate_GeneNames_BoNT.tsv"
# "temp\\Candidate_GeneNames_CaEGTA.tsv"

local({
  
  if(!dir.exists("Figures\\Fig_2BC")) dir.create("Figures\\Fig_2BC", recursive = TRUE)
  if(!dir.exists("temp")) dir.create("temp")
  if(!dir.exists("plots")) dir.create("plots")
  
  library(data.table)
  library(stringr)
  library(limma)
  library(qvalue)
  
  # define variables
  prob.cutoff      <- 0.75
  pvalue.threshold <- 0.01
  fc.threshold     <- 1.2
  
  ph <- fread("temp\\Phosphosites_prepared2.tsv", check.names = TRUE, quote = "")
  
  # subset phosphosites with a given minimal localizatin probability
  message("Sites meeting the localization probability cutoff: ", sum(ph$Localization.prob >= prob.cutoff))
  
  ph  <- ph[Localization.prob >= prob.cutoff]
  dim(ph)

  # subset reporter intensity data and convert to a long format #
  
  # column numbers corresponding to Reporter intensities
  int_cols  <- grep("Reporter\\.intensity\\.corrected\\.[1-6]\\.", names(ph), value = TRUE)
  
  # subset columns from the data frame
  temp <- ph[, c("id", ..int_cols)]

  # convert df to a long format
  temp <- melt(temp, measure.vars = int_cols, variable.factor = FALSE)

  # extract experiment name from 'variable'
  temp[, Experiment := str_match(variable, "[1-6]\\.(.+)___[1-3]$")[, 2]]

  # extract Multiplicity from 'variable'
  temp[, Multiplicity := str_match(variable, "___([1-3])$")[, 2]]

  # extract Channel from 'variable'
  temp[, Channel := str_match(variable, "corrected\\.([0-6])\\.")[, 2]]

  ### remove entries having too many zero values
  n.zero <- temp[, list(n.zero = sum(value == 0)), by = c("id", "Experiment", "Multiplicity")]
  n.zero[, to.delete := n.zero > 3]
  
  # minimal intensiy per experiment/channel
  min.value <- temp[value != 0, list(min.value = min(value, na.rm = TRUE)),by = c("Experiment", "Channel")]
  
  # standard deviation of intensity and intensity at 5% quantile per experiment/channel
  norm.min <- temp[value != 0, list(SD = sd(log2(value), na.rm = TRUE),
                                    percentile = quantile(log2(value), probs = 0.05, na.rm = T)),
                   by = c("Experiment", "Channel")]
  
  temp <- merge(temp, n.zero,    all.x = TRUE)
  temp <- merge(temp, min.value, by = c("Experiment", "Channel"), all.x = TRUE)
  temp <- merge(temp, norm.min, by = c("Experiment", "Channel"), all.x = TRUE)

  # duplicate intensity column
  temp[, value2 := value]
  
  # replace entries excluded from quantification by NAs 
  temp[to.delete == TRUE, value2 := NA_real_]

  # % of missing (zero) values per channel
  temp[!is.na(value2), 100*sum(value2 == 0)/.N, by = c("Experiment", "Channel")] 

  # use probability distribution for missing value imputation
  set.seed(100)
  temp[value2 == 0, value2 := 2^abs(rnorm(n = .N, mean = percentile, sd = 2*SD))]
  
  # remove NAs
  temp <- temp[!is.na(value2)]
  
  # convert channal to a character that can be used as column names
  temp[, variable := paste0("X", Channel)]

  # return to a wide format
  temp <- dcast(temp, Experiment + Multiplicity + n.zero + id ~ variable, value.var = "value2")
  
  # log2-transform
  temp[, paste0("X", 1:6):= lapply(.SD, log2), .SDcols = paste0("X", 1:6)]

  # split df based on Experiment
  temp <- split(temp, f = temp$Experiment)

  # apply Tukey median polishing for normalization (per Experiment)
  # plot intensity distributions before and after normalization
  
  pdf("plots\\Reporter_Int_distributions.pdf")
  line_colors  <- c("black", "blue", "green", "cyan", "red", "orange")
  line_types   <- rep(1, 6)
  take.columns <- grep("^X", names(temp[[1]]), value = TRUE)
  
  for(i in seq_along(temp)){
    
    dat     <- temp[[i]][, ..take.columns]
    dat.med <- medpolish(dat, maxiter = 3, na.rm = TRUE)
    dat.res <- dat.med$residuals
    
    boxplot(dat, main = paste0("Box plot ", names(temp)[i], " before normalization"), xaxt = "n")
    axis(side = 1, at = 1:6, labels = paste0("X", 1:6))
    boxplot(dat.res, main = paste0("Box plot ", names(temp)[i], " after med polishing"), xaxt = "n")
    axis(side = 1, at = 1:6, labels = paste0("X", 1:6))
    
    list_densities <- lapply(dat, density, na.rm = TRUE)
    xmax <- max(unlist(lapply(list_densities, function(dens) max(dens$x))))
    xmin <- min(unlist(lapply(list_densities, function(dens) min(dens$x))))
    ymax <- max(unlist(lapply(list_densities, function(dens) max(dens$y))))
    ymin <- min(unlist(lapply(list_densities, function(dens) min(dens$y))))
    
    plot(list_densities[[1]], type = "n",
         xlim  = c(floor(xmin), ceiling(xmax)),
         ylim = c(floor(ymin), ymax),
         main = paste0("Density plot ", names(temp)[i]))
    for(j in seq_along(list_densities)){
      
      lines(list_densities[[j]], col = line_colors[j], lty = line_types[j], lwd = 2)
      
    }
    legend(x = "topright", legend = paste0("X", 1:6), fill = line_colors)
    
    temp[[i]][, paste0(take.columns, ".norm") := data.table(..dat.res)]
    
  }
  dev.off()
  
  names(temp[["BTX_AD_01"]]) <- gsub("X1", "MockBoNT_03_Mock_01", names(temp[["BTX_AD_01"]]))
  names(temp[["BTX_AD_01"]]) <- gsub("X2", "MockBoNT_03_BoNT_01",  names(temp[["BTX_AD_01"]]))
  names(temp[["BTX_AD_01"]]) <- gsub("X3", "MockBoNT_03_Mock_02", names(temp[["BTX_AD_01"]]))
  names(temp[["BTX_AD_01"]]) <- gsub("X4", "MockBoNT_03_BoNT_02",  names(temp[["BTX_AD_01"]]))
  names(temp[["BTX_AD_01"]]) <- gsub("X5", "MockBoNT_03_BoNT_03",  names(temp[["BTX_AD_01"]]))
  names(temp[["BTX_AD_01"]]) <- gsub("X6", "MockBoNT_03_Mock_03", names(temp[["BTX_AD_01"]]))
  
  names(temp[["BTX_CB_02"]]) <- gsub("X1", "MockBoNT_04_Mock_01", names(temp[["BTX_CB_02"]]))
  names(temp[["BTX_CB_02"]]) <- gsub("X2", "MockBoNT_04_BoNT_01",  names(temp[["BTX_CB_02"]]))
  names(temp[["BTX_CB_02"]]) <- gsub("X3", "MockBoNT_04_Mock_02", names(temp[["BTX_CB_02"]]))
  names(temp[["BTX_CB_02"]]) <- gsub("X4", "MockBoNT_04_BoNT_02",  names(temp[["BTX_CB_02"]]))
  names(temp[["BTX_CB_02"]]) <- gsub("X5", "MockBoNT_04_Mock_03", names(temp[["BTX_CB_02"]]))
  names(temp[["BTX_CB_02"]]) <- gsub("X6", "MockBoNT_04_BoNT_03",  names(temp[["BTX_CB_02"]]))
  
  names(temp[["Ca_EGTA_01"]]) <- gsub("X1", "CaEGTA_01_Ca_01",   names(temp[["Ca_EGTA_01"]]))
  names(temp[["Ca_EGTA_01"]]) <- gsub("X2", "CaEGTA_01_Ca_02",   names(temp[["Ca_EGTA_01"]]))
  names(temp[["Ca_EGTA_01"]]) <- gsub("X3", "CaEGTA_01_Ca_03",   names(temp[["Ca_EGTA_01"]]))
  names(temp[["Ca_EGTA_01"]]) <- gsub("X4", "CaEGTA_01_EGTA_01", names(temp[["Ca_EGTA_01"]]))
  names(temp[["Ca_EGTA_01"]]) <- gsub("X5", "CaEGTA_01_EGTA_02", names(temp[["Ca_EGTA_01"]]))
  names(temp[["Ca_EGTA_01"]]) <- gsub("X6", "CaEGTA_01_EGTA_03", names(temp[["Ca_EGTA_01"]]))
  
  names(temp[["Ca_EGTA_02"]]) <- gsub("X1", "CaEGTA_02_Ca_01",   names(temp[["Ca_EGTA_02"]]))
  names(temp[["Ca_EGTA_02"]]) <- gsub("X2", "CaEGTA_02_EGTA_01", names(temp[["Ca_EGTA_02"]]))
  names(temp[["Ca_EGTA_02"]]) <- gsub("X3", "CaEGTA_02_Ca_02",   names(temp[["Ca_EGTA_02"]]))
  names(temp[["Ca_EGTA_02"]]) <- gsub("X4", "CaEGTA_02_EGTA_02", names(temp[["Ca_EGTA_02"]]))
  names(temp[["Ca_EGTA_02"]]) <- gsub("X5", "CaEGTA_02_Ca_03",   names(temp[["Ca_EGTA_02"]]))
  names(temp[["Ca_EGTA_02"]]) <- gsub("X6", "CaEGTA_02_EGTA_03", names(temp[["Ca_EGTA_02"]]))
  
  exp   <- names(temp)

  if(length(temp) > 1){
    
    for(i in 2:length(temp)){
      
      temp[[1]] <- merge(temp[[1]], temp[[i]][,-c("Experiment", "n.zero")], by = c("id", "Multiplicity"), all = TRUE)
      
    }
    
  }
  dim(temp[[1]])

temp[[1]] <- temp[[1]][, -c("Experiment", "n.zero")]
temp <- temp[[1]]

int_cols <- names(temp)
int_cols <- int_cols[!int_cols %in% c("id", "Multiplicity")]  

int_cols_norm <- sort(int_cols[grepl("\\.norm", int_cols)])
int_cols      <- sort(int_cols[!grepl("\\.norm", int_cols)])

### analysis
# separate Ca_EGTA Experiment and BoNT experiment

# Ca_EGTA
CaEGTA <- temp[, c("id", "Multiplicity", grep("CaEGTA", names(..temp), value = TRUE))]
to.delete <- CaEGTA[, apply(.SD, 1, function(x) any(is.na(x))), .SDcols = c(grep("CaEGTA", names(temp), value = TRUE))]

CaEGTA <- CaEGTA[!to.delete]
dim(CaEGTA)

intensity.names <- grep("CaEGTA", int_cols_norm, value = TRUE)

tr <- c("Ca", "Ca", "Ca", "EGTA", "EGTA", "EGTA", "Ca", "Ca", "Ca", "EGTA", "EGTA", "EGTA")
tr <- factor(tr, levels = c("EGTA", "Ca"))
ex <- as.factor(c(rep(1, 6), rep(2, 6)))

design <- model.matrix(~ ex + tr)

# fit models
fit <- lmFit(CaEGTA[, ..intensity.names], design = design)

n      <- dim(CaEGTA[, ..intensity.names])[1]
fit.eb    <- eBayes(fit)
log2FC <- fit.eb$coef[, "trCa"]
df.0   <- rep(fit.eb$df.prior, n)
df.r   <- fit.eb$df.residual
s2.0   <- rep(fit.eb$s2.prior, n)
s2     <- (fit.eb$sigma)^2
s2.post <- fit.eb$s2.post
t.ord <- fit.eb$coef[, "trCa"]/fit.eb$sigma/fit.eb$stdev.unscaled[, "trCa"]
t.mod <- fit.eb$t[, "trCa"]
p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
p.mod <- fit.eb$p.value[, "trCa"]
q.obj <- qvalue(p.mod)
plot(q.obj)
q.val <- qvalue(p.mod)$qvalues  

sum(q.val < 0.01, na.rm = T)

res.eb.mult <- data.table(log2FC, t.ord, t.mod, p.ord, p.mod, df.r, df.0, s2.0, s2, s2.post, q.val)

# plot p value distribution
pdf("plots\\Pvalue_histogram_CaEGTA.pdf")
hist(res.eb.mult$p.mod, breaks = 100)
dev.off()

# compute adjusted p-values using Benjamini-Hochberg procedure
p.adj <- p.adjust(p.mod, "BH")
hist(p.adj, breaks = 100)

sum(p.adj < 0.01)

CaEGTA <- cbind(CaEGTA, res.eb.mult)
CaEGTA[, p.adj.BH := ..p.adj]

# determine candidates
candidates <- CaEGTA$q.val   < pvalue.threshold &
  (CaEGTA$log2FC  < log2(1/fc.threshold) |
     CaEGTA$log2FC > log2(fc.threshold))
sum(candidates)

CaEGTA[, Candidate := ..candidates]
CaEGTA[, Enriched  := FALSE]
CaEGTA[log2FC > 0, Enriched := TRUE]

# do volcano plot
pdf("plots\\Volcano_qVal_log2FC_CaEGTA.pdf")

plot(y = -log10(CaEGTA$q.val),
     x = CaEGTA$log2FC,
     # pch = 21,
     # bg = "lightgrey",
     #xlim = c(-2, 2),
     type = "n",
     #main = "-log10 Q-value vs log2 Ratio",
     xlab = "log2-Ratio (Ca / EGTA)",
     ylab = "-log10 q-value", cex = 2)
points(y = -log10(CaEGTA$q.val[!candidates]),
       x =  CaEGTA$log2FC[!candidates],
       pch = 21,
       bg  = "lightgrey")
points(y = -log10(CaEGTA$q.val[candidates]),
       x =  CaEGTA$log2FC[candidates],
       bg = "orange",
       pch = 21)

abline(v = log2(1/fc.threshold), col = "blue", lty = 3, lwd = 3)
abline(v = log2(fc.threshold),   col = "blue", lty = 3, lwd = 3)
abline(h = -log10(pvalue.threshold), col = "blue", lty = 3, lwd = 3)

dev.off()

# noTox_BoNT
BoNT <- temp[, c("id", "Multiplicity", grep("MockBoNT", names(..temp), value = TRUE))]
to.delete <- BoNT[, apply(.SD, 1, function(x) any(is.na(x))), .SDcols = c(grep("MockBoNT", names(temp), value = TRUE))]
BoNT <- BoNT[!to.delete]
dim(BoNT)

intensity.names <- grep("MockBoNT", int_cols_norm, value = TRUE)

tr  <- c("BoNT", "BoNT", "BoNT", "Mock", "Mock", "Mock", "BoNT", "BoNT", "BoNT", "Mock", "Mock", "Mock")
tr  <- factor(tr, levels = c("BoNT", "Mock"))
ex  <- as.factor(c(rep(1, 6), rep(2, 6)))

design <- model.matrix(~ ex + tr)

# fit models
fit <- lmFit(BoNT[, ..intensity.names], design = design)

n      <- dim(BoNT[, ..intensity.names])[1]
fit.eb <- eBayes(fit)
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
q.obj <- qvalue(p.mod)
plot(q.obj)
q.val <- qvalue(p.mod)$qvalues  

sum(q.val < 0.01, na.rm = T)

res.eb.mult <- data.table(log2FC, t.ord, t.mod, p.ord, p.mod, df.r, df.0, s2.0, s2, s2.post, q.val)

# plot p value distribution
pdf("plots\\Pvalue_histogram_BoNT.pdf")
hist(res.eb.mult$p.mod, breaks = 100)
dev.off()

# compute adjusted p-values using Benjamini-Hochberg procedure
p.adj <- p.adjust(p.mod, "BH")
hist(p.adj, breaks = 100)

sum(p.adj < 0.01)

BoNT <- cbind(BoNT, res.eb.mult)
BoNT[, p.adj.BH := ..p.adj]

# determine candidates
candidates <- BoNT$q.val   < pvalue.threshold &
  (BoNT$log2FC  < log2(1/fc.threshold) |
   BoNT$log2FC > log2(fc.threshold))
sum(candidates)

BoNT[, Candidate := ..candidates]
BoNT[, Enriched  := FALSE]
BoNT[log2FC > 0, Enriched := TRUE]

# do volcano plot
pdf("plots\\Volcano_qVal_log2FC_BoNT.pdf")

plot(y = -log10(CaEGTA$q.val),
     x = CaEGTA$log2FC,
     # pch = 21,
     # bg = "lightgrey",
     #xlim = c(-2, 2),
     type = "n",
     #main = "-log10 Q-value vs log2 Ratio",
     xlab = "log2-Ratio (Mock / BoNT)",
     ylab = "-log10 q-value")
points(y = -log10(BoNT$q.val[!BoNT$Candidate]),
       x =  BoNT$log2FC[!BoNT$Candidate],
       pch = 21,
       bg  = "lightgrey")
points(y = -log10(BoNT$q.val[BoNT$Candidate]),
       x =  BoNT$log2FC[BoNT$Candidate],
       bg = "orange",
       pch = 21)

abline(v = log2(1/fc.threshold), col = "blue", lty = 3, lwd = 3)
abline(v = log2(fc.threshold),   col = "blue", lty = 3, lwd = 3)
abline(h = -log10(pvalue.threshold), col = "blue", lty = 3, lwd = 3)

dev.off()

# provide figure source data
fwrite(merge(ph[, c("id", "Accession", "Gene.name", "Amino.acid", "Position")],
             CaEGTA[, c("id", "Multiplicity", "log2FC", "p.mod", "q.val", "Candidate")],
             by = "id"),
       "Figures\\Fig_2BC\\Fig_2B_source.txt", sep = "\t")

fwrite(merge(ph[, c("id", "Accession", "Gene.name", "Amino.acid", "Position")],
             BoNT[, c("id", "Multiplicity", "log2FC", "p.mod", "q.val", "Candidate")],
             by = "id"),
       "Figures\\Fig_2BC\\Fig_2C_source.txt", sep = "\t")

take_columns <- c("id", "Protein", "Position", "Amino.acid", "Accession", "Protein.description",
                  "Gene.name", "Fasta.header", "Accession.noIso", "Localization.prob",
                  "Human_Protein", "Human_Position",  "Proteins_sorted",              
                  "Accessions_sorted", "Positions_sorted", "Gene.names_sorted",
                  "Fasta.headers_sorted",  "Modification.window", "Unique.site",                  
                  "Sequence.window_7", "Sequence.window_15",           
                  "Sequence.window_7_noSpace", "REVIEWED", "Stringid", "Putative.kinase",              
                  "GENE", "KIN_ACC_ID", "KIN_ORGANISM", "SUB_GENE", "SUB_ACC_ID", "SUB_ORGANISM",
                  "SUB_MOD_RSD", "DOMAIN", "Human_Residue", "Rat_Stringid_kin",
                  "Rat_string_path", "networkin_score", "nlinks_missing", "netphorest_group",
                  "Rat_Gene_kin", "Kin_mapping", "Kin_gen", "Kin_stringid")

ph <- ph[, ..take_columns]

exclude_columns <- c("t.ord", "t.mod", "p.ord", "df.r", "df.0", "s2.0", "s2", "s2.post")

ph_comb <- merge(CaEGTA[, -c(..exclude_columns)],
                 BoNT[, -c(..exclude_columns)],
                 by = c("id", "Multiplicity"),
                 suffixes = c(".CaEGTA", ".BoNT"), all = TRUE)
ph_comb <- merge(ph, ph_comb, by = "id", all.y = TRUE)
dim(ph_comb)

ph_comb_slim <- ph_comb[, -c(..int_cols, ..int_cols_norm)]
dim(ph_comb_slim)

# long data table format based on log2FC
ph_comb_long <- melt(ph_comb_slim, measure.vars = c("log2FC.CaEGTA", "log2FC.BoNT"), variable.name = "Experiment", value.name = "log2FC")
ph_comb_long[, Experiment := gsub("log2FC\\.", "", Experiment)]
dim(ph_comb_long)

# long data table format based on normalized intensities
ph_comb_int_long <- ph_comb[, -c(..int_cols)]
ph_comb_int_long <- melt(ph_comb_int_long, measure.vars = int_cols_norm,
                         variable.name = "ReplicateID", value.name = "Norm.intensity")

ph_comb_int_long[, Experiment := "CaEGTA"]
ph_comb_int_long[grep("^MockBoNT", ReplicateID), Experiment := "MockBoNT"]
ph_comb_int_long[, Condition := str_match(ReplicateID, "_([^\\d]+)_")[, 2]]
ph_comb_int_long[, Replicate := str_match(ReplicateID, "_(0[123])\\.norm$")[, 2]]

ph_comb_int_long <- ph_comb_int_long[, c("id", "Protein", "Protein.description", "Position", "Multiplicity", "Amino.acid", "Gene.name",
                                         "log2FC.CaEGTA", "log2FC.BoNT", "q.val.CaEGTA", "q.val.BoNT",
                                         "Candidate.CaEGTA", "Candidate.BoNT", "ReplicateID", "Norm.intensity",
                                         "Experiment", "Replicate", "Condition")]

gen_CaEGTA <- unique(ph_comb[Candidate.CaEGTA == TRUE, Gene.name])
gen_CaEGTA <- gen_CaEGTA[!is.na(gen_CaEGTA)]
length(gen_CaEGTA)

gen_BoNT <- unique(ph_comb[Candidate.BoNT == TRUE, Gene.name])
gen_BoNT <- gen_BoNT[!is.na(gen_BoNT)]
length(gen_BoNT)

fwrite(ph_comb,          "temp\\PhPeptIntensities.tsv", sep = "\t", na = "NA")
fwrite(ph_comb_slim,     "temp\\PhPeptIntensities_slim.tsv", sep = "\t", na = "NA")
fwrite(ph_comb_long,     "temp\\PhPeptIntensities_log2FC_long.tsv", sep = "\t", na = "NA")
fwrite(ph_comb_int_long, "temp\\PhPeptIntensities_int_long.tsv", sep = "\t", na = "NA")

fwrite(list(gen_BoNT),    "temp\\Candidate_GeneNames_BoNT.tsv", col.names = FALSE)
fwrite(list(gen_CaEGTA), "temp\\Candidate_GeneNames_CaEGTA.tsv", col.names = FALSE)

})

