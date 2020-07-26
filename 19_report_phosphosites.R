
# Do:
# prepare report

# INPUT:
# "search_results\Synapt_Ph\combined\txt\Phospho (STY)Sites.txt"
# "temp\\PhPeptIntensities4.tsv"

local({
  
  library(data.table)
  library(stringr)
  
  ph <- fread("search_results\\Synapt_Ph\\combined\\txt\\Phospho (STY)Sites.txt", check.names = TRUE)
  
  df <- fread("temp\\PhPeptIntensities4.tsv")
  dim(df)
  names(df)
  
  df <- merge(df , ph[, c("id", "Fasta.headers", "PEP", "Score", "Delta.score")], by = "id")
  dim(df)
  
  take_cols <- c("id", 
                 "Site_id",
                 "Accession",
                 "Gene.name",
                 "Position",
                 "Amino.acid",
                 "Multiplicity",
                 "Localization.prob",
                 "Sequence.window_15",
                 "Protein.description",
                 "Fasta.header",
                 "Fasta.headers",
                 "PEP",
                 "Score",
                 "Delta.score",
                 "Human_Protein",
                 "Human_Position",
                 "networkin_score",
                 "Kin_mapping",
                 "Kin_gen",
                 "netphorest_group"
                 )
  
  test_cols <- c("log2FC.CaEGTA", "p.mod.CaEGTA", "q.val.CaEGTA", "Candidate.CaEGTA",
                 "log2FC.BTX", "p.mod.BTX", "q.val.BTX", "Candidate.BTX",
                 "log2FC.CaEGTA_BTX",
                 "Regulation_group_resolved"
                 )
  
  int_cols <- c("CaEGTA_01_Ca_01",
                "CaEGTA_01_Ca_02",
                "CaEGTA_01_Ca_03",
                "CaEGTA_01_EGTA_01",
                "CaEGTA_01_EGTA_02",
                "CaEGTA_01_EGTA_03",
                "CaEGTA_02_Ca_01",
                "CaEGTA_02_Ca_02",
                "CaEGTA_02_Ca_03",
                "CaEGTA_02_EGTA_01",
                "CaEGTA_02_EGTA_02",
                "CaEGTA_02_EGTA_03",
                "MockBTX_03_Mock_01",
                "MockBTX_03_Mock_02",
                "MockBTX_03_Mock_03",
                "MockBTX_03_BTX_01",
                "MockBTX_03_BTX_02",
                "MockBTX_03_BTX_03",
                "MockBTX_04_Mock_01",
                "MockBTX_04_Mock_02",
                "MockBTX_04_Mock_03",
                "MockBTX_04_BTX_01",
                "MockBTX_04_BTX_02",
                "MockBTX_04_BTX_03"
                )
  
  int_cols.norm <- c("CaEGTA_01_Ca_01.norm",
                     "CaEGTA_01_Ca_02.norm",
                     "CaEGTA_01_Ca_03.norm",
                     "CaEGTA_01_EGTA_01.norm",
                     "CaEGTA_01_EGTA_02.norm",
                     "CaEGTA_01_EGTA_03.norm",
                     "CaEGTA_02_Ca_01.norm",
                     "CaEGTA_02_Ca_02.norm",
                     "CaEGTA_02_Ca_03.norm",
                     "CaEGTA_02_EGTA_01.norm",
                     "CaEGTA_02_EGTA_02.norm",
                     "CaEGTA_02_EGTA_03.norm",
                     "MockBTX_03_Mock_01.norm",
                     "MockBTX_03_Mock_02.norm",
                     "MockBTX_03_Mock_03.norm",
                     "MockBTX_03_BTX_01.norm",
                     "MockBTX_03_BTX_02.norm",
                     "MockBTX_03_BTX_03.norm",
                     "MockBTX_04_Mock_01.norm",
                     "MockBTX_04_Mock_02.norm",
                     "MockBTX_04_Mock_03.norm",
                     "MockBTX_04_BTX_01.norm",
                     "MockBTX_04_BTX_02.norm",
                     "MockBTX_04_BTX_03.norm"
  )
  
  df_sub <- df[, c(..take_cols, ..int_cols, ..int_cols.norm, ..test_cols)]
  
  names(df_sub)[grepl("\\.BTX", names(df_sub))]
  names(df_sub) <- gsub("\\.BTX", "\\.MockBTX", names(df_sub))
  test_cols     <- gsub("\\.BTX", "\\.MockBTX", test_cols)
  
  names(df_sub)[names(df_sub) == "Regulation_group_resolved"] <- "Regulation.group"
  names(df_sub)[names(df_sub) == "Kin_gen"] <- "Kinase_gene"
  names(df_sub)[names(df_sub) == "Kin_mapping"] <- "Kinase_mapping"
  names(df_sub)[names(df_sub) == "log2FC.CaEGTA_BTX"] <- "log2FC.CaEGTA_MockBTX"  
  
  df_sub[, Site_id := paste0(Accession, "_", Amino.acid, Position, "_", Multiplicity)]
  fwrite(df_sub, "SupplData01_Phosphosite_Intensities.txt", sep = "\t")
  
  })
  