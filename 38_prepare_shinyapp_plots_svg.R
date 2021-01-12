# Do:
# prepare svg plots for shiny app

# INPUT:
# "temp\\PhPeptInt_BoNT.tsv"
# "temp\\PhPeptInt_CaEGTA.tsv"
# "temp\\Domains.tsv"

# OUTPUT:
# plots are stored in "ShinyApp\\Lineplots_CaEGTA" and "ShinyApp\\Lineplots_CaEGTA", respectively 
#


local({
  
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(svglite)
  
  if(!dir.exists("ShinyApp\\Lineplots_CaEGTA")) dir.create("ShinyApp\\Lineplots_CaEGTA", recursive = TRUE)
  if(!dir.exists("ShinyApp\\Lineplots_BoNT")) dir.create("ShinyApp\\Lineplots_BoNT", recursive = TRUE)
  
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
    g <- g + scale_color_manual(breaks = c("primary Ca-dependent",
                                           "SV-cycling-dependent",
                                           "not-affected"),
                                values = c("primary Ca-dependent" = "#fcb533",  # yellow #ffdd55
                                           "SV-cycling-dependent" = "#51baf4",  # cyan   #aaeeff
                                           "not-affected"     = scales::alpha("black", 0.6)),
                                labels = c(expression(paste("primary Ca"^"2+", "-dependent")),
                                           "SV-cycling-dependent",
                                           "not-affected"),
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
                             size = 5 * 0.6 * 1.2,
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
                             size = 5 * 0.6 * 1.2,
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
    g <- g + guides(color = guide_legend(title = "Regulation:", override.aes = list(size = 1.2), order = 1, label.hjust = 0),
                    fill  = guide_legend(title = "Domains and Regions:"))
    g
    
  }
  
  # log2FC
  df_bont <- fread("temp\\PhPeptInt_BoNT.tsv")
  df_caegta <- fread("temp\\PhPeptInt_CaEGTA.tsv")
  
  # annotated domains
  dom <- fread("temp\\Domains.tsv")
  
  # create an empty plot
  svglite(paste0("ShinyApp\\Lineplots_CaEGTA\\", "_Empty", ".svg"), width = 10, height = 7)
  lineplotSites(df = df_caegta[Accession == "_Empty"],
                      df_dom = dom[Accession == "_Empty"])
  dev.off()
  
  svglite(paste0("ShinyApp\\Lineplots_BoNT\\", "_Empty", ".svg"), width = 10, height = 7)
  lineplotSites(df = df_caegta[Accession == "_Empty"],
                      df_dom = dom[Accession == "_Empty"])
  dev.off()
  
  # all accessions
  acc <- unique(df_caegta$Accession)
  
  # caegta
  pb <- txtProgressBar(min = 1, max = length(acc), char = "*", style = 3)
  
  for(i in seq_along(acc)){
    
    setTxtProgressBar(pb, i)
    
    svglite(paste0("ShinyApp\\Lineplots_CaEGTA\\", acc[i], ".svg"), width = 10, height = 7)
    print(lineplotSites(df = df_caegta[Accession == acc[i]],
                        df_dom = dom[Accession == acc[i]]))
    dev.off()
    
    
  }
  
  # bont
  acc <- unique(df_bont$Accession)
  
  pb <- txtProgressBar(min = 1, max = length(acc), char = "*", style = 3)
  
  for(i in seq_along(acc)){
    
    setTxtProgressBar(pb, i)
    
    svglite(paste0("ShinyApp\\Lineplots_BoNT\\", acc[i], ".svg"), width = 10, height = 7)
    print(lineplotSites(df = df_bont[Accession == acc[i]],
                        df_dom = dom[Accession == acc[i]]))
    dev.off()
    
    
  }

})
