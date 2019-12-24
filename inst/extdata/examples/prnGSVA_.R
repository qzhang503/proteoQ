# ===================================
# GSVA
# ===================================

## !!!require the brief working example in `?load_expts`

# global option
scale_log2r <- TRUE

# base
prnGSVA(
  min.sz = 10,
  verbose = FALSE,
  parallel.sz = 0,
  mx.diff = TRUE,
  gset_nms = c("go_sets", "kegg_sets"),
)

# row filtration
prnGSVA(
  min.sz = 10,
  verbose = FALSE,
  parallel.sz = 0,
  mx.diff = TRUE,
  gset_nms = "go_sets",
  filter_prots = exprs(prot_n_pep >= 3), 
)

