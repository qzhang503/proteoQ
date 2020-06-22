\donttest{
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
library(NMF)

anal_prnNMF(
  impute_na = FALSE,
  col_group = Group,
  r = c(3:4),
  nrun = 20, 
)

# passing a different `method`
anal_prnNMF(
  impute_na = FALSE,
  col_group = Group,
  method = "lee",
  r = c(3:4),
  nrun = 20, 
  filename = lee.txt,
)

## row filtration and selected samples (proteins)
anal_prnNMF(
  impute_na = FALSE,
  col_select = BI,
  col_group = Group,
  r = c(3:4),
  nrun = 20, 
  filter_prots = exprs(prot_n_pep >= 3),
  filename = bi_npep3.txt,
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

# (`W16_vs_W2.pVal (W16-W2)` now a column key)
anal_prnNMF(
  impute_na = FALSE,
  col_group = Group,
  r = c(3:4),
  nrun = 20, 
  filter_prots_by_npep = exprs(prot_n_pep >= 3), 
  filter_prots_by_pval = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
  filename = pval.txt,
)

## additional row filtration by pVals (impute_na = TRUE)
# if not yet, run prerequisitive NA imputation and corresponding 
# significance tests at `impute_na = TRUE`
pepImp(m = 2, maxit = 2)
prnImp(m = 5, maxit = 5)

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
  filter_prots_by_pval = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
  filename = pval2.txt,
)

## analogous peptides
anal_pepNMF(
  impute_na = TRUE,
  col_group = Group,
  r = c(3:4),
  nrun = 20, 
  filter_prots_by_npep = exprs(prot_n_pep >= 3), 
  filter_prots_by_pval = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
)

anal_pepNMF(
  impute_na = FALSE,
  col_group = Group,
  r = c(3:4),
  nrun = 20, 
  filter_prots_by_npep = exprs(prot_n_pep >= 3), 
  filter_prots_by_pval = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
)


# ===================================
# consensus heat maps
# ===================================
## no NA imputation 
# proteins, all available ranks
library(NMF)

plot_prnNMFCon(
  impute_na = FALSE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 14,
  height = 14,
)

# analogous peptides
plot_pepNMFCon(
  impute_na = FALSE,
  col_select = BI,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Spectral"))(50), 
  width = 10,
  height = 10,
  filename = bi.pdf,
)

# manual selection of input data file(s)
# may be used for optimizing individual plots
plot_prnNMFCon(
  df2 = c("Protein_NMF_Z_rank3_consensus.txt", "Protein_NMF_Z_rank4_consensus.txt"),
  impute_na = FALSE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 14,
  height = 14,
)

## NA imputation 
# proteins, all available ranks
plot_prnNMFCon(
  impute_na = TRUE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 14,
  height = 14,
)

# analogous peptides
plot_pepNMFCon(
  impute_na = TRUE,
  col_select = BI,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 10,
  height = 10,
  filename = bi_con.png,
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
  width = 12,
  height = 12,
)

# manual selection of input data file(s)
# may be used for optimizing individual plots
plot_prnNMFCoef(
  df2 = c("Protein_NMF_Z_rank3_coef.txt"),  
  impute_na = FALSE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  width = 12,
  height = 12,
)

# analogous peptides
plot_pepNMFCoef(
  impute_na = FALSE,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  color = colorRampPalette(brewer.pal(n = 7, name = "Spectral"))(50), 
  width = 12,
  height = 12,
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

# proteins, selected sample(s)
plot_metaNMF(
  impute_na = FALSE,
  col_select = BI_1,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  fontsize = 8,
  fontsize_col = 5,
  cellwidth = 6, 
  filename = bi1.png,
)

# proteins, selected sample(s) and row ordering
plot_metaNMF(
  impute_na = FALSE,
  col_select = BI_1,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  fontsize = 8,
  fontsize_col = 5,
  cellwidth = 6, 
  cluster_rows = FALSE,
  arrange_prots_by = exprs(gene),
  filename = bi1_row_by_genes.png,
)

# manual selection of input .rda file(s)
# may be used for optimizing individual plots
plot_metaNMF(
  df2 = c("Protein_NMF_Z_rank3.rda"),  
  impute_na = FALSE,
  col_select = BI_1,
  annot_cols = c("Color", "Alpha", "Shape"),
  annot_colnames = c("Lab", "Batch", "WHIM"),
  fontsize = 8,
  fontsize_col = 5,
  cellwidth = 6, 
  cluster_rows = FALSE,
  arrange_prots_by = exprs(gene),
  filename = bi1_row_by_genes.png,
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
}
