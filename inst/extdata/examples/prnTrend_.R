# ===================================
# Trend analysis
# ===================================

## !!!require the brief working example in `?load_expts`

## global option
scale_log2r <- TRUE

# ===================================
# Analysis
# ===================================
## base (proteins, with sample order supervision)
anal_prnTrend(
  impute_na = FALSE,
  col_order = Order,
  n_clust = c(5:6), 
)

## against selected samples
anal_prnTrend(
  impute_na = FALSE,
  col_order = Order,
  n_clust = c(5:6), 
  col_select = BI,
  filename = sel.txt,
)

## row filtration (proteins)
anal_prnTrend(
  impute_na = FALSE,
  col_order = Order,
  n_clust = c(5:6), 
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

# (`W16_vs_W2.pVal (W16-W2)` now a column key)
anal_prnTrend(
  impute_na = FALSE,
  col_order = Order,
  n_clust = c(5:6), 
  filter_prots_by_npep = exprs(prot_n_pep >= 3), 
  filter_prots_by_pval = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
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

anal_prnTrend(
  impute_na = TRUE,
  col_order = Order,
  n_clust = c(5:6), 
  filter_prots_by_npep = exprs(prot_n_pep >= 3), 
  filter_prots_by_pval = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
)


# ===================================
# Visualization
# ===================================
## no NA imputation 
plot_prnTrend(
  col_order = Order, 
)

# at specific cluster ID(s)
# (`cluster` is a column key in `Protein_Trend_[...].txt`)
plot_prnTrend(
  impute_na = FALSE, 
  col_order = Order,
  filter_by_clusters = exprs(cluster == 5),
  width = 8, 
  height = 10,
  filename = cl5.png,
)

## NA imputation
# also save as pdf
plot_prnTrend(
  impute_na = TRUE,
  col_order = Order,
  filename = my.pdf,
)

## against selected samples
plot_prnTrend(
  col_order = Order, 
  col_select = BI,
  filename = bi.png,
)

