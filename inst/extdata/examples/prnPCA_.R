# ===================================
# PCA
# ===================================

## !!!require the brief working example in `?load_expts`

# global option
scale_log2r <- TRUE

pepPCA(
  col_select = Select, 
  filter_peps_by = exprs(pep_n_psm >= 10),
  show_ids = FALSE, 
  filename = "peps_rowfil.png",
)

pepPCA(
  col_select = BI, 
  col_shape = Shape,   
  col_color = Alpha, 
  filter_peps_by = exprs(pep_n_psm >= 10),
  show_ids = FALSE, 
  filename = "peps_rowfil_colsel.png",
)

prnPCA(
  col_color = Color,
  col_shape = Shape,
  show_ids = FALSE,
  filter_peps_by = exprs(prot_n_pep >= 5),
  filename = "prns_rowfil.png",
)

# custom palette
prnPCA(
  col_shape = Shape,
  color_brewer = Set1,
  show_ids = FALSE,
  filename = "my_palette.png",
)

# by features
prnPCA(
  type = feats,
  scale_log2r = TRUE,
  filename = "by_feats.png",
)

## Not run: 
prnPCA(
  col_color = "column_key_not_existed",
  col_shape = "another_missing_column_key"
)

## End(Not run)


