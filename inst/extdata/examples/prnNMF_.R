# ===================================
# NMF
# ===================================

## !!!require the brief working example in `?load_expts`

## global option
scale_log2r <- TRUE

library(NMF)

# ===================================
# Analysis
# ===================================
## base (proteins)
anal_prnNMF(
  impute_na = FALSE,
  col_group = Group,
  r = c(3:4),
  nrun = 20, 
)

## row filtration (proteins)
anal_prnNMF(
  impute_na = FALSE,
  col_group = Group,
  r = c(3:4),
  nrun = 20, 
  filter_prots = exprs(prot_n_pep >= 3),
)

## additional row filtration by pVals (proteins, impute_na = FALSE)
# if not yet, run prerequisitive significance tests at `impute_na = FALSE`
pepSig(
  impute_na = FALSE, 
  W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", 
                  "(W2.JHU.TMT2-W2.JHU.TMT1)", 
                  "(W2.PNNL.TMT2-W2.PNNL.TMT1)"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = FALSE)

# (`W2_bat.pVal ((W2.BI.TMT2-W2.BI.TMT1))` is now a column key)
anal_prnNMF(
  impute_na = FALSE,
  col_group = Group,
  r = c(3:4),
  nrun = 20, 
  filter_prots_by_npep = exprs(prot_n_pep >= 3), 
  filter_prots_by_pval = exprs(`W2_bat.pVal ((W2.BI.TMT2-W2.BI.TMT1))` <= 1e-6), 
)

## additional row filtration by pVals (impute_na = TRUE)
# if not yet, run prerequisitive NA imputation
pepImp(m = 2, maxit = 2)
prnImp(m = 5, maxit = 5)

# if not yet, run prerequisitive significance tests at `impute_na = TRUE`
pepSig(
  impute_na = TRUE, 
  W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", 
                  "(W2.JHU.TMT2-W2.JHU.TMT1)", 
                  "(W2.PNNL.TMT2-W2.PNNL.TMT1)"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = TRUE)

anal_prnNMF(
  impute_na = TRUE,
  col_group = Group,
  r = c(3:4),
  nrun = 20, 
  filter_prots_by_npep = exprs(prot_n_pep >= 3), 
  filter_prots_by_pval = exprs(`W2_bat.pVal ((W2.BI.TMT2-W2.BI.TMT1))` <= 1e-6), 
)

## analogous peptides
anal_pepNMF(
  impute_na = TRUE,
  col_group = Group,
  r = c(3:4),
  nrun = 20, 
  filter_peps = exprs(pep_n_psm >= 2),
  filter_prots_by_npep = exprs(prot_n_pep >= 3), 
  filter_prots_by_pval = exprs(`W2_bat.pVal ((W2.BI.TMT2-W2.BI.TMT1))` <= 1e-6), 
)

# ===================================
# consensus heat maps
# ===================================
## no NA imputation 
# proteins, all available ranks
plot_prnNMFCon(
  impute_na = FALSE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# proteins, specific rank(s)
plot_prnNMFCon(
  impute_na = FALSE,
  r = 3,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# analogous peptides
plot_pepNMFCon(
  impute_na = FALSE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

## NA imputation 
# proteins, all available ranks
plot_prnNMFCon(
  impute_na = TRUE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# proteins, specific rank(s)
plot_prnNMFCon(
  impute_na = TRUE,
  r = 3,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# analogous peptides
plot_pepNMFCon(
  impute_na = TRUE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# ===================================
# coefficient heat maps
# ===================================
## no NA imputation 
# proteins, all available ranks
plot_prnNMFCoef(
  impute_na = FALSE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# proteins, specific rank(s)
plot_prnNMFCoef(
  impute_na = FALSE,
  r = 3,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# analogous peptides
plot_pepNMFCoef(
  impute_na = FALSE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

## NA imputation 
# proteins, all available ranks
plot_prnNMFCoef(
  impute_na = TRUE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# proteins, specific rank(s)
plot_prnNMFCoef(
  impute_na = TRUE,
  r = 3,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# analogous peptides
plot_pepNMFCoef(
  impute_na = TRUE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
)

# ===================================
# Metagene heat maps
# ===================================
## no NA imputation 
# proteins, all available ranks
plot_metaNMF(
  impute_na = FALSE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  
  # additional arguments for `pheatmap`
  fontsize = 8,
  fontsize_col = 5,
)

# proteins, specific rank(s)
plot_metaNMF(
  impute_na = FALSE,
  r = 3,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  fontsize = 8,
  fontsize_col = 5,
)

## NA imputation 
# proteins, all available ranks
plot_metaNMF(
  impute_na = TRUE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  fontsize = 8,
  fontsize_col = 5,
)

# proteins, specific rank(s)
plot_metaNMF(
  impute_na = TRUE,
  r = 3,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  fontsize = 8,
  fontsize_col = 5,
)
