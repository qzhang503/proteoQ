# ===================================
# MDS
# ===================================

## !!!require the reduced working example in `?load_expts`

# global option
scale_log2r <- TRUE

# peptide
pepMDS(
  scale_log2r = TRUE,
  col_select = Select, 
  filter_by_npsm = exprs(pep_n_psm >= 10),
  filename = "pepMDS_filtered.png",
)

# protein
prnMDS(
  scale_log2r = TRUE,
  col_color = Color,
  col_shape = Shape,
  show_ids = TRUE,
  filter_by_npep = exprs(prot_n_pep >= 5),
  filename = "prnMDS_filtered.png",
)

# custom palette
prnMDS(
  scale_log2r = TRUE,
  col_shape = Shape,
  color_brewer = Set1,
  show_ids = TRUE,
  filename = "my_palette.png",
)

## Not run: 
prnMDS(
  col_color = "column_key_not_existed",
  col_shape = "another_missing_column_key"
)

## End(Not run)