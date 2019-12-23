# ===================================
# Heat map
# ===================================

## !!!require the brief working example in `?load_expts`

# global option
scale_log2r <- TRUE

# row clustering
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
  width = 18,
  height = 12,
  filter_sp = exprs(species == "human", prot_n_pep >= 2),
  filename = "huprns_npep2.png",
)

# rows ordered by kinase classes then by gene names
# (error if `normPSM(annot_kinases = FALSE, ...)`)
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
  width = 22,
  height = 22,
  filter_kin = exprs(kin_attr, species == "human"),
  arrange_kin = exprs(kin_order, gene),
  filename = "hukins_rows_by_class.png",
)

# `cutree_rows` ignored at `cluster_rows = FALSE`
prnHM(
  scale_log2r = TRUE,
  annot_cols = c("Group"),
  cluster_rows = FALSE,
  clustering_distance_rows  = "maximum",
  cutree_rows = 6,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_row = 3,
  cellwidth = 14,
  width = 22,
  height = 22,
  filename = "cutree_overruled.png",
)


## significance tests for pVal columns
# (`W2_bat.pVal ((W2.BI.TMT2-W2.BI.TMT1))` is a column key in Peptide_pVals.txt)
pepSig(
  impute_na = FALSE, 
  W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", "(W2.JHU.TMT2-W2.JHU.TMT1)", "(W2.PNNL.TMT2-W2.PNNL.TMT1)"], # batch effects
  W2_loc = ~ Term_2["W2.BI-W2.JHU", "W2.BI-W2.PNNL", "W2.JHU-W2.PNNL"], # location effects
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = FALSE)

# now can have row filtration by pVals
prnHM(
  xmin = -1,
  xmax = 1,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = TRUE,
  cutree_rows = 10,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 3,
  cellwidth = 14,
  filter_sp = exprs(species == "human", prot_n_pep >= 2),
  filter_by = exprs(`W2_bat.pVal ((W2.BI.TMT2-W2.BI.TMT1))` <= 1e-6), 
  filename = "pval_cutoff_at_1e6.png", 
)


## impute NA
# peptide
pepImp(m = 2, maxit = 2)

# protein
prnImp(m = 5, maxit = 5)

# with NA imputation
prnHM(
  impute_na = TRUE, 
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
  cellheight = 4,
  width = 18,
  height = 12,
  filter_sp = exprs(species == "human", prot_n_pep >= 2),
  filename = "huprns_npep2_impna.png",
)


# peptides under selected protein(s)
pepHM(
  xmin = -2,
  xmax = 2,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = TRUE,
  annot_rows = c("gene"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  cellwidth = 12,
  cellheight = 12,
  width = 18,
  height = 12,
  filter_by = exprs(gene %in% c("NCL", "Ncl")),
  filename = "ncl_all.png",
)

# rows ordered by gene 
pepHM(
  xmin = -2,
  xmax = 2,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = FALSE,
  annot_rows = c("gene"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  cellwidth = 12,
  cellheight = 12,
  width = 18,
  height = 12,
  
  filter_by = exprs(gene %in% c("NCL", "Ncl")),
  arrange_peps_by = exprs(gene),
  filename = "ncl_rows_by_gene.png",
)

# ordered by sequence 
# (may try `exprs(pep_seq)` if `pep_seq_mod` not a column key in `Peptide.txt`)
pepHM(
  xmin = -2,
  xmax = 2,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = FALSE,
  annot_rows = c("gene"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  cellwidth = 12,
  cellheight = 12,
  width = 18,
  height = 12,
  filter_by = exprs(gene %in% c("NCL", "Ncl")),
  arrange_peps_by = exprs(pep_seq_mod),
  filename = "ncl_rows_by_seq.png",
)

# sample selection
pepHM(
  col_select = BI_1, 
  xmin = -2,
  xmax = 2,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = TRUE,
  annot_rows = c("gene"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  cellwidth = 12,
  cellheight = 12,
  width = 18,
  height = 12,
  filter_by = exprs(gene %in% c("NCL", "Ncl")),
  arrange_peps_by = exprs(gene),  
  filename = "ncl_bi1.png",
)

