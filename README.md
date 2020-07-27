# README

R scripts and data are provided in support of the publication:
2020 Silbern et al

### Before start
Download and uncompress "external.zip" folder.
If you download files from GitHub using git bash, please first install git lfs (https://git-lfs.github.com/) that is aimed at handling large files (e.g. "external.zip" file). If you use web interface for downloading (“download ZIP”), you will download a placeholder for "external.zip" file. To download the actual file, find it in the GitHub repository and click on “View raw”. Save the file in your local folder. Extract files.
IMPORTANT: Files in the "external.zip" are needed for the correct script execution. The content is covered by the third-party licences and is intended for NON-COMMERCIAL use only. The data are provided only for demonstration purporses and are not subject for destribution.   

### Content
"/search_results" folder contains required MaxQuant search results produced in this study. Complete search results and raw mass spectrometry files are deposited at ProteomeXChange (PXD).  
"/prepare_external" folder contains scripts that prepare external data for subsequent analysis.  
"/plot_figures" folder contains scripts that can be used to produce figure used in the publication.
"/comparison_EK" folder contains script that are used to compare data of this study to the previous study by Engholm-Keller et al (2019).   

### Execution
Start R and set working directory to the location of "00_Main.R" file or start "Synapt_Ph" RStudio project directly. Use "00_Main.R" for the correct order of script execution.

### Output
Intermediate results are generated and stored in "/temp" folder.  
Intermediate plots are generated and stored in the "/plots" folder.  
Source data for figures provided in the paper are generated during the script execution and stored in "/Figures".
Figures that were used in the publication (omitting smaller format adjustments) are produced while executing scripts in "plot_figures" figures and stored in "/Figures" folder.
Supplementary Data are gerenrated during the execution and stored in the woking directory directly.  
Intermediate data and plots used for comparsion with Engholm-Keller data are stored in the "/comparison_EK" folder directly. 
Data used by ShinyApp for visualization of log2 fold changes in phosphorylation site intensities on proteins are stored in "/ShinyApp" folder.

### required packages:

- data.table (v. 1.12.8),
- stringr(v. 1.4.0),
- ggplot2 (v. 3.2.1),
- scales (v. 1.1.0),
- qvalue (v. 2.18.0),
- limma (v. 3.42.0),
- ggrepel (v. 0.8.1),
- igraph (v. 1.2.4.2)
- VennDiagram (v. 1.6.20)

### Session Info:
```
R version 3.6.2 (2019-12-12)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 17763)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scales_1.1.0        qvalue_2.18.0       limma_3.42.0        ggrepel_0.8.1       ggplot2_3.2.1      
 [6] igraph_1.2.4.2      VennDiagram_1.6.20  futile.logger_1.4.3 stringr_1.4.0       data.table_1.12.8  

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.3           pillar_1.4.2         compiler_3.6.2       formatR_1.7          plyr_1.8.5          
 [6] futile.options_1.0.1 tools_3.6.2          lifecycle_0.1.0      tibble_2.1.3         gtable_0.3.0        
[11] pkgconfig_2.0.3      rlang_0.4.2          rstudioapi_0.10      withr_2.1.2          dplyr_0.8.3         
[16] tidyselect_0.2.5     glue_1.3.1           R6_2.4.1             purrr_0.3.3          reshape2_1.4.3      
[21] lambda.r_1.2.4       magrittr_1.5         splines_3.6.2        assertthat_0.2.1     colorspace_1.4-1    
[26] stringi_1.4.3        lazyeval_0.2.2       munsell_0.5.0        crayon_1.3.4               
```


