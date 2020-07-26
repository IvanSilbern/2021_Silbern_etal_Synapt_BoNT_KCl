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



