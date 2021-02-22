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

# FASTAs (all platforms)
copy_refseq_hs("~/proteoQ/dbs/fasta/refseq")
copy_refseq_mm("~/proteoQ/dbs/fasta/refseq")

# Optional: UniProt and cRAP FASTAs
copy_uniprot_hsmm("~/proteoQ/dbs/fasta/uniprot")
copy_crap("~/proteoQ/dbs/fasta/crap")

# PSM data (choose one of the platforms)
dat_dir <- "~/proteoQ/examples"

choose_one <- TRUE
if (choose_one) {
  ## Mascot
  copy_mascot_gtmt(dat_dir)
  
  ## MaxQuant
  # copy_maxquant_gtmt(dat_dir)
  
  ## MSFragger
  # copy_msfragger_gtmt(dat_dir)
  
  ## Spectrum Mill
  # copy_specmill_gtmt(dat_dir)
}

# metadata (all platforms)
copy_exptsmry_gtmt(dat_dir)
copy_fracsmry_gtmt(dat_dir)


# ==============================================
# part 2 --- PSMs to peptides and to proteins
# ==============================================
# metadata to workspace
library(proteoQ)
load_expts("~/proteoQ/examples")

# optional: show available MS file names in PSM data
extract_psm_raws()

# PSM standardization
normPSM(
	group_psm_by = pep_seq_mod, 
	group_pep_by = gene, 
	fasta = c("~/proteoQ/dbs/fasta/refseq/refseq_hs_2013_07.fasta", 
	          "~/proteoQ/dbs/fasta/refseq/refseq_mm_2013_07.fasta"), 
	rptr_intco = 1000,
	rm_craps = TRUE,
	rm_krts = FALSE,
	rm_outliers = FALSE, 
	annot_kinases = TRUE, 
	plot_rptr_int = TRUE, 
	plot_log2FC_cv = TRUE, 
)

# optional: cleanup of PSM groups under the same peptide IDs by CV
purgePSM(ymax = 1)

# PSMs to peptides
PSM2Pep()

# peptide data merging
mergePep()

# peptide data standardization by human subset
# (column keys in Peptide.txt suitable for parameterization via varargs of `slice_`)
standPep(
  method_align = MGKernel, 
  range_log2r = c(5, 95), 
  range_int = c(5, 95), 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
  slice_peps_by = exprs(species == "human"),
)

# optional: cleanup of peptide groups under the same protein IDs by CV
purgePep(ymax = 1)

# peptide histograms with logFC scaling
pepHist(
	scale_log2r = TRUE,
	show_curves = TRUE, 
	show_vline = TRUE,
	xmin = -2, 
	xmax = 2,
	ncol = 10, 
)

# peptide histograms without logFC scaling
pepHist(
	scale_log2r = FALSE,
	show_curves = TRUE, 
	show_vline = TRUE,
	xmin = -2, 
	xmax = 2,
	ncol = 10, 
)

# peptides to proteins
Pep2Prn(use_unique_pep = TRUE)

# protein data standardization by human subset
standPrn(
  method_align = MGKernel, 
  range_log2r = c(5, 95), 
  range_int = c(5, 95), 
  n_comp = 2, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
  slice_peps_by = exprs(species == "human"),
)

# protein histograms with logFC scaling
prnHist(
	scale_log2r = TRUE,
	xmin = -2,
	xmax = 2,
	ncol = 10, 
)

# protein histograms without logFC scaling
prnHist(
	scale_log2r = FALSE,
	xmin = -2,
	xmax = 2,
	ncol = 10, 
)

# ==============================================
# part 3 --- significance tests
# ==============================================
## !!! decision on scaling normalization
scale_log2r <- TRUE

## no NA imputation
# peptides
pepSig(
  W2_bat = ~ Term["W2.BI.TMT2-W2.BI.TMT1", 
                  "W2.JHU.TMT2-W2.JHU.TMT1", 
                  "W2.PNNL.TMT2-W2.PNNL.TMT1"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

# proteins
prnSig()

# volcano plots
pepVol()
prnVol()

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  ## NA imputation (optional)
  pepImp(m = 2, maxit = 2)
  prnImp(m = 5, maxit = 5)
  
  pepSig(
    impute_na = TRUE, 
    W2_bat = ~ Term["W2.BI.TMT2-W2.BI.TMT1", 
                    "W2.JHU.TMT2-W2.JHU.TMT1", 
                    "W2.PNNL.TMT2-W2.PNNL.TMT1"],
    W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                      "W2.BI-W2.PNNL", 
                      "W2.JHU-W2.PNNL"],
    W16_vs_W2 = ~ Term_3["W16-W2"], 
  )
  
  prnSig(impute_na = TRUE)
  
  pepVol(impute_na = TRUE)
  prnVol(impute_na = TRUE)
}
## END of DO NOT RUN

# ==============================================
# part 4 --- basic informatics
# ==============================================
## !!! decision on scaling normalization
#  (resetting in case of a new R session)
scale_log2r <- TRUE

### MDS
# peptide
pepMDS(
	show_ids = FALSE,
	width = 8,
	height = 4,	
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

### LDA
# peptide
pepLDA(
  col_group = Group,
  show_ids = FALSE,
)

# protein
prnLDA(
  col_group = Group,
  show_ids = FALSE, 
)

### correlation
# peptide logFC, PNNL subset with sample-order supervision
res <- pepCorr_logFC(
	col_select = PNNL,
	col_order = Order, 
	filename = pnnl_ord.png
)

# head(res)

# protein logFC, WHIM2 subset with supervision
res <- prnCorr_logFC(
	col_select = W2,
	col_order = Group,
	filename = w2_ord.png,
)

# head(res)

### heat map
# protein, human subset with row-clustering and subtrees
prnHM(
	xmin = -1, 
	xmax = 1, 
	xmargin = 0.1, 
	annot_cols = c("Group", "Color", "Alpha", "Shape"), 
	annot_colnames = c("Group", "Lab", "Batch", "WHIM"), 
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
	xmin = -1, 
	xmax = 1, 
	xmargin = 0.1, 
	annot_cols = c("Group", "Color", "Alpha", "Shape"), 
	annot_colnames = c("Group", "Lab", "Batch", "WHIM"), 
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

# Make sure that Cytoscape is open.
# Also note that "Protein_Trend_Z_nclust5.txt" is a secondary data file.
cluego(
  df2 = Protein_Trend_Z_nclust5.txt, 
  species = c(human = "Homo Sapiens"), 
  n_clust = c(3, 5)
)

### NMF
library(NMF)

# analysis, protein
anal_prnNMF(
  impute_na = FALSE,
  col_group = Group,
  r = c(5:6),
  nrun = 20, 
)

# consensus heat maps, protein
plot_prnNMFCon(
  impute_na = FALSE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
	width = 12,
	height = 12,
)

# coefficients heat maps, protein
plot_prnNMFCoef(
  impute_na = FALSE,
	annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 12,
  height = 12,
)

# metagene heat maps, protein
plot_metaNMF(
  impute_na = FALSE,
	annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  fontsize = 8,
  fontsize_col = 5,
)

### GSPA
# analysis, protein
prnGSPA(
	pval_cutoff = 5E-2, # protein pVal threshold
	logFC_cutoff = log2(1.2), # protein log2FC threshold
	gspval_cutoff = 5E-2, # gene-set threshold
)

# volcano plot visualization, protein
gspaMap(
	gspval_cutoff = 5E-2, # gene set threshold
	gslogFC_cutoff = log2(1.2), # gene set log2FC threshold
	show_sig = pVal, 
	yco = 0.05, 
)

### GSPA distance heat map and network 
# human subset
prnGSPAHM(
	filter2_sp = exprs(start_with_str("hs", term)), 
	annot_cols = "ess_idx",
	annot_colnames = "Eset index",
	annot_rows = "ess_size", 
	filename = show_connectivity_at_large_dist.png,
)

prnGSPAHM(
	filter2_by = exprs(distance <= .33),
	filter2_sp = exprs(start_with_str("hs", term)), 
	annot_cols = "ess_idx",
	annot_colnames = "Eset index",
	annot_rows = "ess_size", 
	filename = show_human_redundancy.png,
)

### GSVA
prnGSVA(
  min.sz = 10,
  verbose = FALSE,
  parallel.sz = 0,
  mx.diff = TRUE,
)

### GSEA
prnGSEA(
  var_cutoff = 0, 
  pval_cutoff = 1, 
  logFC_cutoff = log2(1), 
  filter_by_sp = exprs(species == "human"), 
)

### STRING database
anal_prnString(
  db_nms = c(prepString(human), prepString(mouse)),
  score_cutoff = .9,
  filter_prots_by = exprs(prot_n_pep >= 2),
  filename = download_and_analysis_in_one_pot.tsv,
)

# human ppi only
anal_prnString(
  db_nms = prepString(human),
  score_cutoff = .9,
  filter_by_sp = exprs(species == "human"),
  filter_prots_by = exprs(prot_n_pep >= 2),
	filename = hu.tsv,
)

