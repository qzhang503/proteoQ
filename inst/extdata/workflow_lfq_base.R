# ==============================================
# part 0 --- installation
# ==============================================
# proteoQ package
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("qzhang503/proteoQ")

# ==============================================
# part 1 --- setup
# ==============================================
library(proteoQDA)

# FASTA (all platforms)
copy_uniprot_hsmm("~/proteoQ/dbs/fasta/uniprot")
copy_crap("~/proteoQ/dbs/fasta/crap")

# PSM data
# (demo only: the same TMT from CPTAC searched with LFQ)
# (some columns in msms.txt omitted)
dat_dir <- "~/proteoQ/examples"

choose_one <- TRUE
if (choose_one) {
  ## Mascot
  copy_mascot_plfq(dat_dir)
  
  ## MaxQuant
  # copy_maxquant_plfq(dat_dir)
  
  ## MSFragger
  # copy_msfragger_plfq(dat_dir)
  
  ## Spectrum Mill
  # copy_specmill_plfq(dat_dir)
}

# metadata
copy_exptsmry_plfq(dat_dir)
copy_fracsmry_plfq(dat_dir)


# ==============================================
# part 2 --- PSMs to peptides and to proteins
# ==============================================
# metadata to workspace
library(proteoQ)
load_expts("~/proteoQ/examples")

# optional: show available MS file names in PSM data
extract_psm_raws()

normPSM(
	group_psm_by = pep_seq_mod, 
	group_pep_by = gene, 
	fasta = c("~/proteoQ/dbs/fasta/uniprot/uniprot_hsmm_2020_03.fasta", 
	          "~/proteoQ/dbs/fasta/crap/crap.fasta"),
	rm_craps = TRUE,
	rm_krts = FALSE,
	annot_kinases = TRUE, 
	plot_rptr_int = TRUE, 
)

## No applicable with LFQ
# purgePSM(ymax = 1)

# PSMs to peptides
PSM2Pep()

# peptide data merging
mergePep()

# peptide histograms (median-centered)
pepHist(
  scale_log2r = TRUE,
  cut_points = c(mean_lint = seq(6, 10, 1)),
  xmin = -4, 
  xmax = 4,
  ncol = 6, 
)

# peptide data aligned by human subsets and 
# with a minimal of three measures (out of six) in intensity
standPep(
  method_align = MGKernel, 
  range_log2r = c(5, 95), 
  range_int = c(5, 95), 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
  slice_peps_by = exprs(species == "human", count_nna >= 3),
)

# updated histograms
pepHist(
  cut_points = c(mean_lint = seq(6, 10, 1)),
  xmin = -4, 
  xmax = 4,
  ncol = 6, 
  filter_peps_by = exprs(species == "human", count_nna >= 3),
  filename = aligned_by_hs.png,
)

# optional: cleanup of peptide groups 
#   under the same protein IDs by CV
purgePep(pt_cv = .95, ymax = 4)

# peptides to proteins
Pep2Prn(use_unique_pep = TRUE)

# protein histograms (median-centered)
prnHist(
  scale_log2r = TRUE,
  cut_points = c(mean_lint = seq(6, 10, 1)),
  xmin = -5, 
  xmax = 5,
  ncol = 6, 
)

# protein data standardization by human subsets and 
# with a minimal of three measures (out of six) in intensity
standPrn(
  method_align = MGKernel, 
  range_log2r = c(5, 95), 
  range_int = c(5, 95), 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
  slice_peps_by = exprs(species == "human", count_nna > 3),
)

prnHist(
  scale_log2r = TRUE,
  cut_points = c(mean_lint = seq(6, 10, 1)),
  xmin = -5, 
  xmax = 5,
  ncol = 6, 
  filter_prns_by = exprs(species == "human", count_nna > 3),
  filename = aligned_by_hs.png,
)



# ==============================================
# part 3 --- significance tests
# ==============================================
## !!! decision on scaling normalization
scale_log2r <- TRUE

## no NA imputation
# peptides
pepSig(
  fml_1 = ~ Term["BI-JHU", 
							 "JHU-PNNL", 
							 "(BI+JHU)/2-PNNL"],
)

# proteins
prnSig()

# volcano plots
pepVol()
prnVol()

prnVol(
  fml_nms = "fml_1",
  xmax = 10, 
  filename = custom.png,
)


# ==============================================
# part 4 --- basic informatics
# ==============================================
## !!! decision on scaling normalization
scale_log2r <- TRUE

### MDS
# peptide
pepMDS(
	show_ids = FALSE,
	width = 10,
	height = 3.75,	
)

# protein
prnMDS(
	show_ids = FALSE, 
	width = 8,
	height = 4,
)

### PCA
# peptide
pepPCA(
	show_ids = FALSE,
)

# protein
prnPCA(
	show_ids = FALSE, 
)

### correlation
pepCorr_logFC()

# protein logFC
prnCorr_logFC()

### heat map
# protein, human subset with row-clustering and subtrees
prnHM(
  xmin = -2, 
  xmax = 2, 
  xmargin = 0.2, 
  annot_cols = c("Group"), 
  annot_colnames = c("Lab"), 
  cluster_rows = TRUE, 
  cutree_rows = 10, 
  show_rownames = FALSE, 
  show_colnames = TRUE, 
  fontsize_row = 3, 
  cellwidth = 14, 
  filter_sp = exprs(species == "human"), 
  filename = hu.png,
)

# protein, kniase subset with row supervision by kinase classes
prnHM(
	xmin = -2, 
	xmax = 2, 
	xmargin = 0.2, 
	annot_cols = c("Group"), 
	annot_colnames = c("Lab"), 
	cluster_rows = FALSE, 
	annot_rows = c("kin_class"), 
	show_rownames = TRUE, 
	show_colnames = TRUE, 
	fontsize_row = 2, 
	cellheight = 2, 
	cellwidth = 14, 
	filter_kin = exprs(kin_attr, species == "human"),
	arrange_kin = exprs(kin_order, gene),
	filename = hukin.png, 
)

### clustering (cemans, kmeans etc.)
# protein analysis, row filtration and sample-order supervision
anal_prnTrend(
  col_order = Order,
  n_clust = c(5:6), 
  filter_by_npep = exprs(prot_n_pep >= 2),
)

# protein visualization, sample-order supervision
plot_prnTrend(
  col_order = Order,
)

### NMF
library(NMF)

# analysis, protein
anal_prnNMF(
  impute_na = FALSE,
  col_group = Group,
  r = c(5:6),
  nrun = 50, 
)

# consensus heat maps, protein
plot_prnNMFCon(
  impute_na = FALSE,
	annot_cols = c("Group"), 
	annot_colnames = c("Lab"), 
	width = 12,
	height = 12,
)

# coefficients heat maps, protein
plot_prnNMFCoef(
  impute_na = FALSE,
  annot_cols = c("Group"), 
  annot_colnames = c("Lab"), 
  width = 12,
  height = 12,
)

# metagene heat maps, protein
plot_metaNMF(
  impute_na = FALSE,
  annot_cols = c("Group"), 
  annot_colnames = c("Lab"), 
  fontsize = 8,
  fontsize_col = 5,
)

### GSPA
# analysis, protein
prnGSPA(
  impute_na = FALSE,
	pval_cutoff = 5E-2, # protein pVal threshold
	logFC_cutoff = log2(1.2), # protein log2FC threshold
	gspval_cutoff = 5E-2, # gene-set threshold
	gset_nms = c("go_sets", "c2_msig"),
)

# volcano plot visualization, protein
gspaMap(
  impute_na = FALSE,
	gspval_cutoff = 5E-2, # gene set threshold
	gslogFC_cutoff = log2(1.2), # gene set log2FC threshold
	topn = 100, 
	show_sig = pVal, 
	yco = 0.05, 
)

### GSPA distance heat map and network 
# human subset
prnGSPAHM(
	# filter2_sp = exprs(start_with_str("hs", term)), 
	annot_cols = "ess_idx",
	annot_colnames = "Eset index",
	annot_rows = "ess_size", 
	filename = show_connectivity_at_large_dist.png,
)

### GSVA
prnGSVA(
  impute_na = FALSE,
  min.sz = 10,
  verbose = FALSE,
  parallel.sz = 0,
  mx.diff = TRUE,
	gset_nms = c("go_sets", "c2_msig"),
)

### GSEA
prnGSEA(
  var_cutoff = 0, 
  pval_cutoff = 1, 
  logFC_cutoff = log2(1), 
  # filter_by_sp = exprs(species == "human"), 
)

### STRING database
anal_prnString(
  db_nms = c(prepString(human), prepString(mouse)),
  score_cutoff = .9,
  filter_prots_by = exprs(prot_n_pep >= 2),
  filename = download_and_analysis_in_one_pot.tsv,
)

