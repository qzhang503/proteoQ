# ===================================
# NMF analysis
# ===================================

## !!!require the brief working example in `?load_expts`

# global option
scale_log2r <- TRUE

# required
library(NMF)

## analysis
# peptides, at two different r(ank)s
anal_pepNMF(
  col_group = Group,
  r = c(6, 8),
  nrun = 200,
  filter_peps_by_npsm = exprs(pep_n_psm >= 2),
)

# proteins, over a range of ranks
anal_prnNMF(
  impute_na = FALSE,
  col_group = Group,
  r = c(5:8),
  nrun = 200, 
  filter_prots_by_npep = exprs(prot_n_pep >= 2),
)

## consensus heat maps
# peptides, at specific ranks
plot_pepNMFCon(
  r = c(5, 6),
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# peptides, at all available ranks
plot_pepNMFCon(
  impute_na = FALSE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# proteins, at specific ranks
plot_prnNMFCon(
  r = c(7:8),
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# proteins, at all available ranks
plot_prnNMFCon(
  impute_na = FALSE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10
)

## coefficient heat maps
# peptides, at all ranks
plot_pepNMFCoef(
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10
)

# proteins, at all ranks
plot_prnNMFCoef(
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10
)

## metagene heat maps 
# at all available ranks
# additional arguments for `pheatmap`
plot_metaNMF(
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  fontsize = 8,
  fontsize_col = 5
)

