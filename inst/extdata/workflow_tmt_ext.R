# ==============================================
# part 0 --- installation
# ==============================================
# proteoQ package
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("qzhang503/proteoQ")

# data package: FASTA and PSM examples
devtools::install_github("qzhang503/proteoQDA")


# ==============================================
# part 1 --- setup
# ==============================================
library(proteoQDA)

# FASTA (all platforms)
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

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  ## Mascot
  # columns keys in PSM files suitable for varargs of `filter_`
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
    
    filter_psms = exprs(pep_expect <= .1, pep_score >= 15), 
    filter_more_psms = exprs(pep_rank == 1),
  )
  
	# peptide sequences with different side-chain 
  #   modifications treated as the same species
	normPSM(
		group_psm_by = pep_seq, 
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
		
		filter_psms = exprs(pep_expect <= .1, pep_score >= 15), 
		filter_more_psms = exprs(pep_rank == 1),
	)
	
	## MaxQuant (e.g. `PEP` is column key)
	normPSM(
	  ...,
	  filter_psms_by = exprs(PEP <= 0.1), 
	)
	
	# MSFragger (`Expectation` etc. are column keys)
	normPSM(
	  filter_psms_at = rlang::exprs(Expectation <= 0.1, Retention >= 50), 
	  ..., 
	)
	
	# Spectrum Mill (e.g. `score` is column key)
	normPSM(
	  ...,
	  filter_psms_by = exprs(score >= 10), 
	)	
	
	# multiple accessions
	normPSM(
	  ..., 
	  fasta = c("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_09.fasta", 
	            "~/proteoQ/dbs/fasta/custom/my_accession.fasta",
	            "~/proteoQ/dbs/fasta/crap/crap.fasta"), 
	)
	
	# ok to use fasta(s) at higher organism(s)
	normPSM(
	  ..., 
	  fasta = c("~/proteoQ/dbs/fasta/uniprot/uniprot_all_species.fasta", 
	            "~/proteoQ/dbs/fasta/crap/crap.fasta"), 
	)
	
	# custom entrez db (see also ?Uni2Entrez, ?Ref2Entrez)
	normPSM(
	  ..., 
	  entrez = c("~\\proteoQ\\dbs\\entrez\\uniprot_entrez_hs.rds", 
	            "~\\proteoQ\\dbs\\entrez\\uniprot_entrez_mm.rds"), 
	)
}
## END of DO NOT RUN

# optional: cleanup of PSM groups under the same peptide IDs by CV
purgePSM()

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  purgePSM(pt_cv = 0.95)
}
## END of DO NOT RUN

# PSMs to peptides
# (reaches a set of common column keys)
PSM2Pep()

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  # need `Raw peptide match data` when exporting Mascot PSMs
  PSM2Pep(filter_ms1int = rlang::exprs(pep_tot_int >= 1E4))
}
## END of DO NOT RUN

# peptide data merging
mergePep()

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  # columns keys in TMTset1_LCMSinj1_Peptide_N.txt etc. 
  #  suitable for varargs of `filter_`
  mergePep(filter_peps_at = exprs(pep_len <= 100))
  
  ## Splines at `cut_points`
  mergePep(cut_points = c(mean_lint = seq(4, 7, .5)))
}
## END of DO NOT RUN

# peptide data standardization by human subset
# (column keys in Peptide.txt suitable for 
#   parameterization via varargs of `slice_`)
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

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  # ---------------------------------------------
  # see ?standPep for mixed-bed examples and more 
  # ---------------------------------------------
  
  # re-normalization against selected samples with human data
  standPep(
    method_align = MGKernel, 
    range_log2r = c(5, 95), 
    range_int = c(5, 95), 
    n_comp = 3, 
    seed = 749662, 
    maxit = 200, 
    epsilon = 1e-05, 
    slice_peps_by = exprs(species == "human"),
    col_select = Select_sub,	
  )
}
## END of DO NOT RUN

# optional: cleanup of peptide groups under the same protein IDs by CV
purgePep()

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  purgePep(pt_cv = 0.95)
}
## END of DO NOT RUN

# peptide histograms with logFC scaling
pepHist(
	scale_log2r = TRUE,
	show_curves = TRUE, 
	show_vline = TRUE,
	xmin = -1, 
	xmax = 1,
	ncol = 10, 
)

# peptide histograms without logFC scaling
pepHist(
	scale_log2r = FALSE,
	show_curves = TRUE, 
	show_vline = TRUE,
	xmin = -1, 
	xmax = 1,
	ncol = 10, 
)

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  # row filtration
  pepHist(
    scale_log2r = TRUE,
    show_curves = TRUE, 
    show_vline = TRUE,
    xmin = -1, 
    xmax = 1,
    ncol = 10, 
    filter_peps_by = exprs(pep_n_psm >= 10), 
    filename = npsm10.png, 
  )
  
  # `BI` subset without scaling
  pepHist(
   scale_log2r = FALSE, 
   col_select = BI,
   filename = bi_n.png, 
  )
  
  # `BI` subset with scaling
  pepHist(
   scale_log2r = TRUE, 
   col_select = BI,
   filename = bi_z.png, 
  )  
}
## END of DO NOT RUN

# peptides to proteins
Pep2Prn(use_unique_pep = TRUE)

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  # column keys in `Peptide.txt` suitable for varargs of `filter_`
  Pep2Prn(
    use_unique_pep = TRUE,
    filter_peps = exprs(prot_n_pep >= 2),
  )
  
  # phosphopeptide subset
  Pep2Prn(
    use_unique_pep = TRUE,
    filter_peps = exprs(pep_mod_sty == TRUE),
  )
  
  # Splines at `cut_points`
  Pep2Prn(cut_points = c(mean_lint = seq(4, 7, 0.5)))
}
## END of DO NOT RUN

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

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
	# ---------------------------------------------
  # see ?standPrn for mixed-bed examples and more
  # ---------------------------------------------
  
  # re-normalization against selected samples with human data
  standPrn(
    col_select = Select_sub, 
    method_align = MGKernel, 
    range_log2r = c(5, 95), 
    range_int = c(5, 95), 
    n_comp = 2, 
    seed = 749662, 
    maxit = 200, 
    epsilon = 1e-05, 
    slice_peps_by = exprs(species == "human"),
  )
}
## END of DO NOT RUN

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

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  # row filtration
  prnHist(
    scale_log2r = TRUE,
    xmin = -2,
    xmax = 2,
    ncol = 10, 
    filter_by = exprs(prot_n_psm >= 20), 
    filename = npsm20.png, 	
  )

  # `BI` subset without scaling
  prnHist(
    scale_log2r = FALSE, 
    col_select = BI,
    filename = bi_n.png, 
  )
  
  # `BI` subset with scaling
  prnHist(
    scale_log2r = TRUE, 
    col_select = BI,
    filename = bi_z.png, 
  )  
}
## END of DO NOT RUN


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

dontrun <- TRUE
if (!dontrun) {
  pepVol(
    filter_cdk1 = exprs(gene == "CDK1"), 
    fml_nms = "W16_vs_W2", 
    filename = CDK1.png,
  )
}
## END of DO NOT RUN

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
scale_log2r <- TRUE

### MDS
# peptide
pepMDS(
	show_ids = FALSE,
	width = 10,
	height = 3.75,	
)

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  # row filtration
  pepMDS(
  	show_ids = FALSE,
  	width = 10,
  	height = 3.75,
  	filter_peps_by = exprs(pep_n_psm >= 20), 
  	filename = npsm20.png,	
  )
  
  # exclude higher-abundance ones
  pepMDS(
    show_ids = FALSE,
    width = 10,
    height = 3.75,
    filter_peps_by = exprs(pep_n_psm <= 8), 
    filename = npsm_le_8.png,	
  )
  
  # JHU subset
  pepMDS(
  	col_select = JHU,
  	show_ids = FALSE,
  	width = 10,
  	height = 3.75,
  	filename = jhu.png,
  )
  
  # JHU subset with new aesthetics
  pepMDS(
    col_select = JHU,
    col_fill = Shape, # WHIMs  
    col_size = Alpha, # batches
    show_ids = FALSE,
  	width = 10,
  	height = 3.75,	
  	filename = jhu_new_aes.png,
  )  
}
## END of DO NOT RUN

# protein
prnMDS(
	show_ids = FALSE, 
	width = 8,
	height = 4,
)

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  # row filtration
  prnMDS(
  	show_ids = FALSE, 
  	width = 8,
  	height = 4,
  	filter_prots_by = exprs(prot_n_pep > 5), 
  	filename = npep5.png,		
  )
  
  # JHU subset
  prnMDS(
  	col_select = JHU,
  	show_ids = FALSE,
  	width = 8,
  	height = 4,
  	filename = jhu.png,
  )
  
  # JHU subset with new aesthetics
  prnMDS(
  	col_select = JHU,
  	col_color = Alpha,
  	col_alpha = Size,
  	show_ids = FALSE,
  	width = 8,
  	height = 4,
  	filename = jhu_new_aes.png,
  )  
}
## END of DO NOT RUN

### PCA
# peptide
pepPCA(
	show_ids = FALSE,
)

# protein
prnPCA(
	show_ids = FALSE, 
)

### Euclidean distance
# peptide, PNNL subset
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
	filename = pnnl.png, 
)

# protein, PNNL subset
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
	filename = pnnl.png, 
)

### correlation
# peptide logFC, PNNL subset with sample-order supervision
pepCorr_logFC(
	col_select = PNNL,
	col_order = Order, 
	filename = pnnl_ord.png
)

# protein logFC, WHIM2 subset with supervision
prnCorr_logFC(
	col_select = W2,
	col_order = Group,
	filename = w2_ord.png,
)

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

### NMF
library(NMF)

# analysis, protein
anal_prnNMF(
  impute_na = FALSE,
  col_group = Group,
  rank = c(5:6),
  nrun = 20, 
)

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  # row filtration
  anal_prnNMF(
    impute_na = FALSE,
    col_group = Group,
    rank = c(5:6),
    nrun = 20, 
    filter_prots = exprs(prot_n_pep >= 2),
  )
}
## END of DO NOT RUN

# consensus heat maps, protein
plot_prnNMFCon(
  impute_na = FALSE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
	width = 10,
	height = 10,
)

# coefficients heat maps, protein
plot_prnNMFCoef(
  impute_na = FALSE,
	annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
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
  impute_na = FALSE,
	pval_cutoff = 5E-2, # protein pVal threshold
	logFC_cutoff = log2(1.2), # protein log2FC threshold
	gspval_cutoff = 5E-2, # gene-set threshold
	gset_nms = c("go_sets", "kegg_sets"),
)

# volcano plot visualization, protein
gspaMap(
  impute_na = FALSE,
	gspval_cutoff = 5E-2, # gene set threshold
	gslogFC_cutoff = log2(1.2), # gene set log2FC threshold
	show_sig = pVal, 
	yco = 0.05, 
)

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
	# individualized parameters and with row filtration
	prnGSPA(
	  impute_na = FALSE,
		fml_nms = c("W2_bat", "W2_loc", "W16_vs_W2"), # formulae used in `prnSig()`
		pval_cutoff = c(5E-2, 5E-2, 1E-5), # respective protein pVal cut-offs for each formula
		logFC_cutoff = c(log2(1.1), log2(1.1), log2(1.2)), # protein log2FC cut-offs
		gspval_cutoff = c(5E-3, 5E-3, 1E-6), # gene-set pVal cut-offs 
		max_size = c(Inf, Inf, 120), # maxixmum sizes of gene sets for consideration
		
		gset_nms = c("go_sets", "c2_msig"), 
		filter_prots = exprs(prot_n_pep >= 2), 
	)
	
	gspaMap(
	  impute_na = FALSE,
		fml_nms = c("W2_bat", "W2_loc", "W16_vs_W2"),
		pval_cutoff = c(5E-2, 5E-2, 1E-10),
		logFC_cutoff = log2(1.2),
		show_sig = pVal,
		yco = 0.05,
		filter_prots = exprs(prot_n_pep >= 2),
	)
}
## END of DO NOT RUN

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
  impute_na = FALSE,
  min.sz = 10,
  verbose = FALSE,
  parallel.sz = 0,
  mx.diff = TRUE,
  gset_nms = "go_sets",
)

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  # row filtration
  prnGSVA(
    impute_na = FALSE,
    min.sz = 10,
    verbose = FALSE,
    parallel.sz = 0,
    mx.diff = TRUE,
    gset_nms = c("go_sets", "kegg_sets"),
    filter_prots = exprs(prot_n_pep >= 3), 
  )
}
## END of DO NOT RUN

### GSEA
prnGSEA(
  var_cutoff = 0, 
  pval_cutoff = 1, 
  logFC_cutoff = log2(1), 
  filter_by_sp = exprs(species == "human"), 
)

## DO NOT RUN
dontrun <- TRUE
if (!dontrun) {
  # prefiltration by variances, pVals and logFCs...
  prnGSEA(
    var_cutoff = 0.5,     
    pval_cutoff = 5E-2,
    logFC_cutoff = log2(1.2),
    filter_by_sp = exprs(species == "human"), 
  )
}
## END of DO NOT RUN

### STRING
anal_prnString(
  db_nms = c(prepString(human), prepString(mouse)),
  score_cutoff = .9,
  filter_prots_by = exprs(prot_n_pep >= 2),
  filename = download_and_analysis_in_one_pot.tsv,
)

