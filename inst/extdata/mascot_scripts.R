## part 0 --- installation
# proteoQ package
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biobase", "Mfuzz", "limma"))

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("qzhang503/proteoQ")

# data package with fasta and Mascot examples
devtools::install_github("qzhang503/proteoQDA")


## part 1 --- setup
# fasta files to database directory
library(proteoQDA)
copy_refseq_hs("~\\proteoQ\\dbs\\fasta\\refseq")
copy_refseq_mm("~\\proteoQ\\dbs\\fasta\\refseq")

# examplary PSM data to working directory
dir.create("C:\\The\\Mascot\\Example", recursive = TRUE, showWarnings = FALSE)
dat_dir <- "C:\\The\\Mascot\\Example"
cptac_csv_1(dat_dir)

# metadata to working directory
cptac_expt_1(dat_dir)
cptac_frac_1(dat_dir)

# metadata upload
library(proteoQ)
load_expts()


## part 2 --- summarization of PSMs to peptides and proteins
# PSM tables
normPSM(
	group_psm_by = pep_seq, 
	group_pep_by = gene, 
	fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta", 
						"~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"), 
	rptr_intco = 3000,
	rm_craps = TRUE,
	rm_krts = FALSE,
	rm_outliers = FALSE, 
	annot_kinases = TRUE, 
	plot_rptr_int = TRUE, 
	plot_log2FC_cv = TRUE, 
	
	filter_peps = exprs(pep_expect <= .1, pep_score >= 15), 
	filter_by_more = exprs(pep_rank == 1),
)

# optional: purge of PSM groups under the same peptide IDs
purgePSM(max_cv = 0.5, min_n = 2)

# peptide tables
normPep(
	method_psm_pep = median, 
	method_align = MGKernel, 
	range_log2r = c(5, 95), 
	range_int = c(5, 95), 
	n_comp = 3, 
	seed = 749662, 
	maxit = 200, 
	epsilon = 1e-05, 
)

# optional: purge of peptide groups under the same protein IDs
purgePep(max_cv = .5, min_n = 2)

# peptide histograms with logFC scaling
pepHist(
	scale_log2r = TRUE,
	show_curves = TRUE, 
	show_vline = TRUE,
	xmin = -1, 
	xmax = 1,
	ncol = 10, 
)

# peptide histograms with scaling and filtration
pepHist(
	scale_log2r = TRUE,
	show_curves = TRUE, 
	show_vline = TRUE,
	xmin = -1, 
	xmax = 1,
	ncol = 10, 
	
	filter_by = exprs(pep_n_psm >= 10), 
	filename = "pepHist_npsm10.png", 
)

# peptide histograms of `BI` subset without scaling
pepHist(
 scale_log2r = FALSE, 
 col_select = BI,
 filename = Hist_BI_N.png, 
)

# peptide histograms of `BI` subset with scaling
pepHist(
 scale_log2r = TRUE, 
 col_select = BI,
 filename = Hist_BI_Z.png, 
)

# renormalization of peptide data for selected samples
normPep(
	method_psm_pep = median, 
	method_align = MGKernel, 
	range_log2r = c(5, 95), 
	range_int = c(5, 95), 
	n_comp = 3, 
	seed = 749662, 
	maxit = 200, 
	epsilon = 1e-05, 
	col_refit = Select_sub,
)

# protein tables
normPrn(
	method_pep_prn = median, 
	method_align = MGKernel, 
	range_log2r = c(20, 95), 
	range_int = c(5, 95), 
	n_comp = 2, 
	seed = 749662, 
	maxit = 200, 
	epsilon = 1e-05, 
)

# examplary protein tables with data prefiltration
normPrn(
	method_pep_prn = median, 
	method_align = MGKernel, 
	range_log2r = c(20, 95), 
	range_int = c(5, 95), 
	n_comp = 2, 
	seed = 749662, 
	maxit = 200, 
	epsilon = 1e-05, 
	
	filter_by = exprs(prot_n_psm >= 5, prot_n_pep >= 2),  
)

# protein histograms with logFC scaling
prnHist(
	scale_log2r = TRUE,
	xmin = -2,
	xmax = 2,
	ncol = 10, 
)

# protein histograms with logFC scaling and filtration
prnHist(
	scale_log2r = TRUE,
	xmin = -2,
	xmax = 2,
	ncol = 10, 
	filter_by = exprs(prot_n_psm >= 20), 
	filename = "prnHist_npsm_20.png", 	
)

# protein histograms without logFC scaling
prnHist(
	scale_log2r = FALSE,
	xmin = -2,
	xmax = 2,
	ncol = 10, 
)

# protein histograms of `BI` subset with logFC scaling
prnHist(
	scale_log2r = TRUE,
	col_select = BI, 
	xmin = -2,
	xmax = 2,
	ncol = 5, 
	filename = Hist_BI_Z.png, 
)

# protein histograms of `BI` subset without logFC scaling
prnHist(
	scale_log2r = FALSE,
	col_select = BI, 
	xmin = -2,
	xmax = 2,
	ncol = 5, 
	filename = Hist_BI_N.png, 
)

# global setting in scaling normalization
scale_log2r = TRUE


## part 3 --- basic informatics
# optional: peptide NA imputation (warning: may take a while)
pepImp(m = 2, maxit = 2)

# optional: protein NA imputation
prnImp(m = 5, maxit = 5)

# peptide MDS
pepMDS(
	show_ids = FALSE,
	width = 10,
	height = 3.75,
)

# peptide MDS with filtration
pepMDS(
	show_ids = FALSE,
	width = 10,
	height = 3.75,
	filter_by = exprs(pep_n_psm >= 20), 
	filename = "pepMDS_npsm_20.png",	
)

# peptide MDS of JHU subset
pepMDS(
	col_select = JHU,
	show_ids = FALSE,
	width = 10,
	height = 3.75,
	filename = "MDS_JHU.png",
)

# peptide MDS of JHU subset with new aesthetics
pepMDS(
  col_select = JHU,
  col_fill = Shape, # WHIMs  
  col_size = Alpha, # batches
  show_ids = FALSE,
	width = 10,
	height = 3.75,
	filename = "MDS_JHU_new_aes.png",
)

# protein MDS
prnMDS(
	show_ids = FALSE, 
	width = 8,
	height = 4,
)

# protein MDS with filtration
prnMDS(
	show_ids = FALSE, 
	width = 8,
	height = 4,
	filter_by = exprs(prot_n_pep > 5), 
	filename = "prnMDS_npep_5.png",		
)

# protein MDS of JHU subset
prnMDS(
	col_select = JHU,
	show_ids = FALSE,
	width = 8,
	height = 4,
	filename = "MDS_JHU.png",
)

# protein MDS of JHU subset with new aesthetics
prnMDS(
	col_select = JHU,
	col_color = Alpha,
	col_alpha = Size,
	show_ids = FALSE,
	width = 8,
	height = 4,
	filename = "MDS_JHU_new_aes.png",
)

# peptide PCA
pepPCA(
	show_ids = FALSE,
)

# protein PCA
prnPCA(
	show_ids = FALSE, 
)

# peptide EucDist of PNNL subset
pepEucDist(
	col_select = PNNL,
	annot_cols = c("Shape", "Alpha"),
	annot_colnames = c("WHIM", "Batch"), 

	display_numbers = TRUE, 
	number_color = "grey30", 
	number_format = "%.1f",
	
	clustering_distance_rows = "euclidean", 
	clustering_distance_cols = "euclidean", 
	
	fontsize = 16, 
	fontsize_row = 20, 
	fontsize_col = 20, 
	fontsize_number = 8, 
	
	cluster_rows = TRUE,
	show_rownames = TRUE,
	show_colnames = TRUE,
	border_color = "grey60", 
	cellwidth = 24, 
	cellheight = 24, 
	width = 14,
	height = 12, 
	filename = "EucDist_PNNL.png", 
)

# protein EucDist of PNNL subset
prnEucDist(
	col_select = PNNL,
	annot_cols = c("Shape", "Alpha"),
	annot_colnames = c("WHIM", "Batch"), 

	display_numbers = TRUE, 
	number_color = "grey30", 
	number_format = "%.1f",
	
	clustering_distance_rows = "euclidean", 
	clustering_distance_cols = "euclidean", 
	
	fontsize = 16, 
	fontsize_row = 20, 
	fontsize_col = 20, 
	fontsize_number = 8, 
	
	cluster_rows = TRUE,
	show_rownames = TRUE,
	show_colnames = TRUE,
	border_color = "grey60", 
	cellwidth = 24, 
	cellheight = 24, 
	width = 14,
	height = 12, 
	filename = "EucDist_PNNL.png", 
)

# peptide logFC correlation of PNNL subset with sample-order supervision
pepCorr_logFC(
	col_select = PNNL,
	col_order = Order, 
	filename = PNNL_ord.png
)

# peptide log10-intensity correlation of PNNL subset with supervision
pepCorr_logInt(
	col_select = PNNL,
	col_order = Order,
	filename = PNNL_ordint.png,
)

# protein logFC correlation of WHIM2 subset with supervision
prnCorr_logFC(
	col_select = W2,
	col_order = Group,
	filename = W2_ord.png,
)

# protein heat maps of human subset with row-clustering and subtrees
prnHM(
	xmin = -1, 
	xmax = 1, 
	x_margin = 0.1, 
	annot_cols = c("Group", "Color", "Alpha", "Shape"), 
	annot_colnames = c("Group", "Lab", "Batch", "WHIM"), 
	cluster_rows = TRUE, 
	cutree_rows = 10, 
	show_rownames = FALSE, 
	show_colnames = TRUE, 
	fontsize_row = 3, 
	cellwidth = 14, 
	width = 18, 
	height = 12, 
	filter_sp = exprs(species == "human"), 
)

# protein heat maps of kniases subset with row supervision by kinase classes
prnHM(
	xmin = -1, 
	xmax = 1, 
	x_margin = 0.1, 
	annot_cols = c("Group", "Color", "Alpha", "Shape"), 
	annot_colnames = c("Group", "Lab", "Batch", "WHIM"), 
	cluster_rows = FALSE, 
	annot_rows = c("kin_class"), 
	show_rownames = TRUE, 
	show_colnames = TRUE, 
	fontsize_row = 2, 
	cellheight = 2, 
	cellwidth = 14, 
	width = 22, 
	height = 22, 
	filter_kin = exprs(kin_attr),
	arrange_kin = exprs(kin_order, gene),
	filename = "kin_by_class.png", 
)

# protein significance tests
prnSig(
	impute_na = FALSE, 
	W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", "(W2.JHU.TMT2-W2.JHU.TMT1)", "(W2.PNNL.TMT2-W2.PNNL.TMT1)"], # batch effects
	W2_loc = ~ Term_2["W2.BI-W2.JHU", "W2.BI-W2.PNNL", "W2.JHU-W2.PNNL"], # location effects
)

# protein volcano plots
prnVol()

# peptide significance tests
pepSig(
	impute_na = FALSE, 
	W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", "(W2.JHU.TMT2-W2.JHU.TMT1)", "(W2.PNNL.TMT2-W2.PNNL.TMT1)"], # batch effects
	W2_loc = ~ Term_2["W2.BI-W2.JHU", "W2.BI-W2.PNNL", "W2.JHU-W2.PNNL"], # location effects
)

# peptide volcano plots
pepVol()

# protein trend analysis with filtration and sample-order supervision
anal_prnTrend(
  col_order = Order,
  n_clust = c(5:8), 
  filter_by_npep = exprs(prot_n_pep >= 2),
)

# protein trend visualization with sample-order supervision
plot_prnTrend(
  col_order = Order,
  n_clust = c(5:6), 
)

# NMF package
library(NMF)

# protein NMF: analysis with filtration
anal_prnNMF(
  impute_na = FALSE,
  col_group = Group,
  r = c(5:8),
  nrun = 200, 
  filter_by_npep = exprs(prot_n_pep >= 2),
)

# protein NMF: consensus heat maps
plot_prnNMFCon(
  impute_na = FALSE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# protein NMF: coefficients heat maps
plot_prnNMFCoef(
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# protein NMF: metagene heat maps
plot_metaNMF(
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  fontsize = 8,
  fontsize_col = 5,
)

# protein GSPA
prnGSPA(
  impute_na = FALSE, 
  pval_cutoff = 5E-2, 
  gset_nm = c("go_sets", "kegg_sets"), 
)

# GSPA under volcano plots
gspaMap(
	show_labels = TRUE, 
	pval_cutoff = 5E-3, 
	show_sig = pVal, 
	yco = 0.05, 
)

# STRING database
getStringDB(
  db_path = "~\\proteoQ\\dbs\\string",
  score_cutoff = .9,
  adjP = FALSE, 
)
