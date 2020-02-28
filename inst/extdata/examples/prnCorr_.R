\donttest{
# ===================================
# Correlation
# ===================================

## !!!require the brief working example in `?load_expts`

## global option
scale_log2r <- TRUE

# peptide log2FC with sample ID ordering
#(no more than 40 samples for visualization)
pepCorr_logFC(
  col_select = BI,
  col_order = Order, 
  width = 25,
  height = 25,
  filter_peps_by = exprs(pep_n_psm >= 3),
  filename = bi_npsm3.png,
)

# protein log2FC
prnCorr_logFC(
  col_select = W2,
  col_order = Order, 
  width = 40,
  height = 40,
  filter_prots = exprs(prot_n_pep >= 2),
  filename = w2_npep2.png,
)


## Not run: 
# at most 40 samples
pepCorr_logFC(
  col_order = Order, 
  width = 40,
  height = 40,
  filter_peps_by = exprs(pep_n_psm >= 3),
  filename = too_many_cols.png,
)

# interplex comparison of peptide intensity
# (modest correlation in interplex reporter-ion intensity at data-dependant acquistion)
pepCorr_logInt(
  width = 10,
  height = 10,
  filter_peps_by = exprs(pep_n_psm >= 3),
  filename = pepcorr_int_npsm3.png,
)

# protein intensity
prnCorr_logInt(
  width = 10,
  height = 10,
  filter_prots_by = exprs(prot_n_pep >= 5),
  filename = prncorr_int_npep5.png,
)
## End(Not run)
}

