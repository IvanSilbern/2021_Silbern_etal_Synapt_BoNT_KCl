# DO:
#
# plot proteins with regulated sites in CaEGTA or MockBTX experiments

# INPUT:
# "temp\\PhPeptInt_BTX.tsv"
# "temp\\PhPeptInt_CaEGTA.tsv"
# "temp\\Domains.tsv"

# OUTPUT:
# "SupplData5_Phosphoproteins_regulated_CaEGTA.pdf"
# "SupplData6_Phosphoproteins_regulated_MockBTX.pdf"

local({
  
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(svglite)
  
  lineplotSites <- function(df, df_dom, pos = NA, title = "", subtitle = "",
                            y_limit = NA, y_limit_avg = 2.0,
                            selected_line_col = "red"){
    
    # Function is designed to plot log2 fold changes of the phosphosites
    # as vertical lines across the whole length of a protein
    # required is a data.table containing "Position", "Amino.acid", "Multiplicity",
    # "Protein.description", "Gene.name", "Protein", "Siteid", "Enriched"
    # Name of the column containing log2FC values has to be specified through log2FC argument.
    # Log2FC with the highest magnitude will be selected if there are different Multiplicities for a given site.
    # If pos = NA, protein length will be assumed as highest Position + 5
    # if y_limit is specified, y scale limits will be set to +/- y_limit
    # if y_limit is NA, y limits will be set to
    # maximum of 1.1 maximum absolute log2FC in the data set or y_lim_avg
    # selected_line_col argument specifices line color of selected site
    
    # removes NAs
    df <- df[!is.na(df$log2FC)]
    
    # data table is not empty
    if(nrow(df) == 0) return(print(ggplot()))
    
    # protein annotation
    # pname <- df$Protein.description[1]
    # gn    <- df$Gene.name[1]
    # acc   <- df$Protein[1]
    
    plen <- df$Length[1]
    if(is.na(plen) || length(plen) == 0) plen <- max(df$Position) + 5
    
    if(is.na(y_limit)){
      
      y_limit <- round(max(abs(df$log2FC)), 1)*1.1
      y_limit <- max(y_limit, y_limit_avg)
      
    }
    
    # highlighting of selected site
    # df[, Selected := FALSE]
    # if(!is.na(pos) && length(pos) > 0) df[df$Position == pos[1], Selected := TRUE]
    
    # starting point to draw a line
    # (should align with the width of the rectangle
    # representing protein sequence)
    df[df$log2FC > 0, Begin := -y_limit*0.02]
    df[df$log2FC < 0, Begin :=  y_limit*0.02]
    
    # impute minimum log2FC
    df[which(df$log2FC > 0 & (df$log2FC <  y_limit*0.05)), log2FC := y_limit*0.05]
    df[which(df$log2FC < 0 & (df$log2FC > -y_limit*0.05)), log2FC := -y_limit*0.05]
    
    # Make a long table format for ggplot
    df <- melt(df, measure.vars = c("Begin", "log2FC"))
    #df[, Siteid2 := paste0(Siteid, "_", Multiplicity)]
    
    ###
    
    g <- ggplot(df, aes(x = Position, y = value, group = Siteid, label = Siteid, col = Regulation_group))
    g <- g + scale_x_continuous(limits = c(1, plen))
    g <- g + scale_y_continuous(limits = c(-y_limit, y_limit))
    g <- g + scale_color_manual(breaks = c("Ca-compensating",
                                           "Ca-dependent",
                                           "Ca-enhancing",
                                           "Cycling-dependent",
                                           "not-regulated"),
                                values = c("Ca-compensating"   = "#66c837",  # green  #abc837
                                           "Ca-dependent"      = "#ffa655",  # yellow #ffdd55
                                           "Ca-enhancing"      = "#ff7469",  # salmon #ff8c69
                                           "Cycling-dependent" = "#31a9ff",  # cyan   #aaeeff
                                           "not-regulated"     = scales::alpha("black", 0.6)),
                                labels = c("Ca-compensating",
                                           "Ca-dependent",
                                           "Ca-enhancing",
                                           "Cycling-dependent",
                                           "not-regulated"),
                                name = "")
    # g <- g + geom_segment(x = 1, xend = plen,
    #                       y = 0, yend = 0,
    #                       color = "black", size = 2)
    
    g <- g + geom_rect(xmin = 1, xmax = plen,
                       ymin = -y_limit*0.02, ymax = y_limit*0.02,
                       alpha = 0.2, fill = "grey65", color = "black")
    g <- g + geom_rect(data = df_dom,
                       aes(xmin = Start, xmax = End, fill = Description_wrapped),
                       ymin = -y_limit*0.02, ymax = y_limit*0.02, alpha = 0.8,
                       inherit.aes = FALSE)
    g <- g + geom_segment(x = 1, xend = 1,
                          y = -y_limit*0.02, yend = y_limit*0.02,
                          color = "darkgrey", size = 2, lineend = "round")
    g <- g + geom_segment(x = plen, xend = plen,
                          y = -y_limit*0.02, yend = y_limit*0.02,
                          color = "darkgrey", size = 2, lineend = "round")
    g <- g + geom_text(x = 0, y = 0, hjust = 1.4, label = "N", fontface = "bold", color = "black")
    g <- g + geom_text(x = plen, y = 0, hjust = -0.6, label = "C", fontface = "bold", color = "black")
    g <- g + geom_line(alpha = 0.8)
    #g <- g + geom_line(data = df[df$Selected], col = selected_line_col, alpha = 0.7)
    g <- g + geom_text_repel(data = df[df$variable != "Begin" & df$Enriched],
                             segment.color = "grey65",
                             size = 5 * 0.6,
                             #color = enriched_text_col,
                             fontface = "bold",
                             ylim = c(0.2, y_limit - 0.05*y_limit),
                             nudge_x = plen*0.01,
                             nudge_y = y_limit * 0.1,
                             min.segment.length = 0,
                             show.legend = FALSE
    )
    g <- g + geom_text_repel(data = df[df$variable != "Begin" & !df$Enriched],
                             segment.color = "grey65",
                             size = 5 * 0.6,
                             #color = enriched_text_col,
                             fontface = "bold",
                             ylim = c(-0.2, -(y_limit - 0.05*y_limit)),
                             nudge_x = -plen*0.01,
                             nudge_y = -y_limit * 0.1,
                             min.segment.length = 0,
                             show.legend = FALSE
    )
    g <- g + ggtitle(title, subtitle = subtitle)
    g <- g + ylab("log2 Fold Change")
    g <- g + guides(color = guide_legend(title = "Regulation:", override.aes = list(size = 1.2), order = 1),
                    fill  = guide_legend(title = "Domains and Regions:"))
    # g <- g + theme(legend.justification = c(1, 0),
    #                legend.direction     = "horizontal",
    #                legend.position      = c(0.98, 0.92),
    #                legend.background    = element_rect(fill = scales::alpha("white", 0.2)),
    #                #legend.background    = element_rect(fill = NA),
    #                legend.key           = element_rect(fill = "white"),
    #                legend.title         = element_blank(),
    #                legend.text          = element_text(face = "bold"),
    #                plot.subtitle        = element_text(face = "bold"))
    
    g
    
  }
  
  # log2FC
  df_btx <- fread("temp\\PhPeptInt_BTX.tsv")
  df_caegta <- fread("temp\\PhPeptInt_CaEGTA.tsv")
  
  # annotated domains
  dom <- fread("temp\\Domains.tsv")
  
  # sort by gene name
  df_caegta <- df_caegta[order(Gene.name)]
  
  # all accessions
  acc <- unique(df_caegta$Accession[df_caegta$Regulation_group != "not-regulated"])
  length(acc)
  
  # caegta
  pb <- txtProgressBar(min = 1, max = length(acc), char = "*", style = 3)
  
  pdf("SupplData5_Phosphoproteins_regulated_CaEGTA.pdf", width = 10, height = 6)
  
  for(i in seq_along(acc)){
    
    setTxtProgressBar(pb, i)
    
    print(lineplotSites(df = df_caegta[Accession == acc[i]],
                        df_dom = dom[Accession == acc[i]],
                        title = paste0(df_caegta$Gene.name[df_caegta$Accession == acc[i]][1]),
                        subtitle = paste0("Uniprot Accession: ",
                                          df_caegta$Accession[df_caegta$Accession == acc[i]][1],
                                          "\n",
                                          df_caegta$Protein.description[df_caegta$Accession == acc[i]][1])
                        ))
    
    
    
  }
  
  dev.off()
  
  
  
  # btx
  
  # sort by gene name
  df_btx <- df_btx[order(Gene.name)]
  
  acc <- unique(df_btx$Accession[df_btx$Regulation_group != "not-regulated"])
  length(acc)
  
  # btx
  pb <- txtProgressBar(min = 1, max = length(acc), char = "*", style = 3)
  
  pdf("SupplData6_Phosphoproteins_regulated_MockBTX.pdf", width = 10, height = 6)
  
  for(i in seq_along(acc)){
    
    setTxtProgressBar(pb, i)
    
    print(lineplotSites(df = df_btx[Accession == acc[i]],
                        df_dom = dom[Accession == acc[i]],
                        title = paste0(df_btx$Gene.name[df_btx$Accession == acc[i]][1]),
                        subtitle = paste0("Uniprot Accession: ",
                                          df_btx$Accession[df_btx$Accession == acc[i]][1],
                                          "\n",
                                          df_btx$Protein.description[df_btx$Accession == acc[i]][1])
    ))
    
    
  }
  
  dev.off()
  
})
