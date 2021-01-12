
# Do:
# prepare report

# INPUT:
# "search_results\Synapt_Ph\combined\txt\Phospho (STY)Sites.txt"
# "temp\\PhPeptIntensities4.tsv"

# OUTPUT:
# "SupplData\\SupplData01_Phosphosite_Intensities.txt"

local({
  
  if(!dir.exists("SupplData")) dir.create("SupplData")
  
  library(data.table)
  library(stringr)
  
  ph <- fread("search_results\\Synapt_Ph\\combined\\txt\\Phospho (STY)Sites.txt.gz", check.names = TRUE)
  
  df <- fread("temp\\PhPeptIntensities6.tsv")
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
                 "netphorest_group",
                 "occupancy_mock",
                 "PP1_substr",
                 "CN_substr"
                 )
  
  test_cols <- c("log2FC.CaEGTA", "p.mod.CaEGTA", "q.val.CaEGTA", "Candidate.CaEGTA",
                 "log2FC.BoNT", "p.mod.BoNT", "q.val.BoNT", "Candidate.BoNT",
                 "log2FC.CaEGTA_BoNT",
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
                "MockBoNT_03_Mock_01",
                "MockBoNT_03_Mock_02",
                "MockBoNT_03_Mock_03",
                "MockBoNT_03_BoNT_01",
                "MockBoNT_03_BoNT_02",
                "MockBoNT_03_BoNT_03",
                "MockBoNT_04_Mock_01",
                "MockBoNT_04_Mock_02",
                "MockBoNT_04_Mock_03",
                "MockBoNT_04_BoNT_01",
                "MockBoNT_04_BoNT_02",
                "MockBoNT_04_BoNT_03"
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
                     "MockBoNT_03_Mock_01.norm",
                     "MockBoNT_03_Mock_02.norm",
                     "MockBoNT_03_Mock_03.norm",
                     "MockBoNT_03_BoNT_01.norm",
                     "MockBoNT_03_BoNT_02.norm",
                     "MockBoNT_03_BoNT_03.norm",
                     "MockBoNT_04_Mock_01.norm",
                     "MockBoNT_04_Mock_02.norm",
                     "MockBoNT_04_Mock_03.norm",
                     "MockBoNT_04_BoNT_01.norm",
                     "MockBoNT_04_BoNT_02.norm",
                     "MockBoNT_04_BoNT_03.norm"
  )
  
  df_sub <- df[, c(..take_cols, ..int_cols, ..int_cols.norm, ..test_cols)]
  
  names(df_sub)[grepl("\\.BoNT", names(df_sub))]
  names(df_sub) <- gsub("\\.BoNT", "\\.MockBoNT", names(df_sub))
  test_cols     <- gsub("\\.BoNT", "\\.MockBoNT", test_cols)
  
  names(df_sub)[names(df_sub) == "Regulation_group_resolved"] <- "Regulation.group"
  names(df_sub)[names(df_sub) == "Kin_gen"] <- "Kinase_gene"
  names(df_sub)[names(df_sub) == "Kin_mapping"] <- "Kinase_mapping"
  names(df_sub)[names(df_sub) == "log2FC.CaEGTA_BoNT"] <- "log2FC.CaEGTA_MockBoNT"  
  
  df_sub[, Site_id := paste0(Accession, "_", Amino.acid, Position, "_", Multiplicity)]
  fwrite(df_sub, "SupplData\\SupplData01_Phosphosite_Intensities.txt", sep = "\t")
  
  })
  