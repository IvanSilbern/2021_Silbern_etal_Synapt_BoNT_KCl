
#
# If running first time or there is no alignment file external\data_BLAST\Human_vs_Rat_Networkin.txt.bz2
# prepare Rat/Human protein alignment
# use BLAST+ alignment tool (run locally, v. 2.7.1)
# to align Rat protein sequences
# to human protein sequences used by NetworKIN
# 
# 1. create a sequence database:
#     makeblastdb -in "external\data_Uniprot\201902_uniprot-rat.fasta" -dbtype prot
# 2. Alignment
#     blastp.exe -max_target_seqs 40 -num_threads 6 -word_size 6 -gapopen 12 -gapextend 1 -query external\data_Networkin\9606.protein.sequences.v9.0.fa -db external\data_Uniprot\201902_uniprot-rat.fasta -out temp\Human_vs_Rat_Networkin.asn -outfmt 11
# 3. Format output table
#     blast_formatter.exe -archive temp\Human_vs_Rat_Networkin.asn -outfmt "6 qacc sacc evalue bitscore nident pident mismatch gapopen gaps qstart qend sstart send qlen slen length qseq sseq qframe sframe" -out external\data_BLAST\Human_vs_Rat_Networkin.txt
# 4. Bz2-compress Human_vs_Rat_Networkin.txt and move to external\data_BLAST\
#

# prepare external data
source("prepare_external\\01_up_mappingtable_stringid.R")
source("prepare_external\\02_prepare_phsitedb_kin_substrate_fasta.R")
source("prepare_external\\03_prepare_phsitedb_kin_substrate_stringid.R")
source("prepare_external\\04_prepare_stringdb_graph.R")
source("prepare_external\\05_prepare_human_rat_blast_stringid.R")
source("prepare_external\\06_prepare_netphorest_kin_map.R")
source("prepare_external\\07_prepare_domain_data.R")

# analyse data
source("01_prepare_phosphosite_table.R")
source("02_extract_fasta_sequences.R")

# If run first time or there is no alignment file external\data_BLAST\PhSeq_PhsitePlus_align.txt
# prepare sequence window alignment
# use BLAST+ alignment tool (run locally, v. 2.7.1)
# 1. Use sequence windows from Kinase_Substrate_Dataset.txt
#     to create a database:
#     makeblastdb -in "Phsitedb_kin_sub.fasta" -dbtype prot
# 2. Align sequence windows
#     blastp.exe -max_target_seqs 20 -num_threads 6 -query temp\Sequence_window_15.fasta -db temp\Phsitedb_kin_sub.fasta -out temp\PhSeq_PhsitePlus.asn -outfmt 11 
# 3. Convert BLAST output into table format
#     blast_formatter.exe -archive temp\PhSeq_PhsitePlus.asn -outfmt "6 qacc sacc evalue bitscore nident pident mismatch gapopen gaps qstart qend sstart send qseq sseq qframe sframe" -out external\data_BLAST\PhSeq_PhsitePlus_align.txt
# 4. Copy PhSeq_PhsitePlus_align.txt into external\data_BLAST\
#

source("03_prepare_aligned_phsites.R")
source("04_add_kin_substrate_phsp.R")
source("05_convert_rat_human_sites.R")


# If run first time or there is no "external\\data_NetworKIN\\networking_predictions.tsv.bz2"
# NetworKIN predictions
# Use Human_Stringid_PhSites.tsv as input
# (only "Human_Protein", "Human_Position", "Human_Residue" columns)
# Submit to https://networkin.info/index_ht.shtml (high-throughput workflow)
# Save unmapped sites under "external\data_NetworKIN\Networkin_unmapped_sites.txt"
# re-run source("05_convert_rat_human_sites.R")
# re-submit data to NetworKIN
# bz2-compress and save under "external\\data_NetworKIN\\networking_predictions.tsv.bz2"

source("06_add_stringids.R")
source("07_add_networkin_predictions.R")
source("08_testing_limma.R")
source("09_netphorest_group_count_reg_sites.R")
source("10_netphorest_group_count_fisher_test_human_bcgr.R")
source("11_netphorest_group_count_fisher_test_CaEGTA_MockBoNT.R")
source("12_GOBP_netphorest_group_enrichm.R")
source("13_comparison_MockBoNT_CaEGTA.R")
source("14_restructure_tables.R")
source("15_add_functional_annotation.R")
source("16_resolve_site_regulation.R")
source("17_regulation_groups_proteins.R")
source("18_regulation_groups_kinases.R")
source("19_report_phosphosites.R")
source("20_analysis_unbound_MockBoNT_AD.R")
source("21_analysis_unbound_MockBoNT_CB.R")
source("22_analysis_supernatant_Cbotulinum.R")
source("23_analysis_test_hela_phospho.R")
source("24_prepare_shinyapp_tables.R")
source("25_prepare_shinyapp_Rdata.R")
source("26_plot_selected_genes_svg.R")
source("27_plot_selected_genes_no_legend_svg.R")

# Uncomment to run
# Plot generation might take long!
# source("28_plot_regulated_sites_pdf.R")
# source("29_prepare_shinyapp_plots_svg.R")

source("30_Annotation_Manual.R")
source("31_reactome_pathway.R")

# compare to Engholm-Keller data
source("comparison_EK//01_prepare_phosphosite_table_EK.R")
source("comparison_EK//02_prepare_phosphosite_table_IS.R")
source("comparison_EK//03_PhSiteContrasts_EK.R")
source("comparison_EK//04_comparison_log2FC.R")
source("comparison_EK//05_comparison_counts.R")

# plot figures
source("plot_figures\\Fig_2A.R")
source("plot_figures\\Fig_2D.R")
source("plot_figures\\Fig_2BC.R")
source("plot_figures\\Fig_2EF.R")
source("plot_figures\\Fig_3A.R")
source("plot_figures\\Fig_4A.R")
source("plot_figures\\Fig_4B.R")
source("plot_figures\\Fig_4C.R")
source("plot_figures\\SupplFig_2.R")
source("plot_figures\\SupplFig_5.R")
source("plot_figures\\SupplFig_15.R")
source("plot_figures\\SupplFig_16ABC.R")
source("plot_figures\\SupplFig_16DEF.R")
source("plot_figures\\SupplFig_17.R")
