# ===================================
# Trend analysis
# ===================================

## !!!require the brief working example in `?load_expts`

# global option
scale_log2r <- TRUE

# analysis 
# (with sample order supervision)
anal_prnTrend(
  col_order = Order,
  n_clust = c(5:8), 
  filter_by_npep = exprs(prot_n_pep >= 2),
)

# visualization
plot_prnTrend(
  col_order = Order,
  n_clust = c(5:6), 
)

# subsets visualization
# (proteins with four or more identifying peptides)
plot_prnTrend(
  col_order = Order,
  filter_by_npep = exprs(prot_n_pep >= 4),
)

