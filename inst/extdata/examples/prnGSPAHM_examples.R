# ===================================
# Heat maps of enriched gene sets
# ===================================
# distance heat map and network of GSPA terms
# a `term` is a subset of an `ess_term` if the distance is zero
# `ess_idx` is a column key in `essmap_.*.csv`
# `ess_size` is a column key in metadata file `essmeta_.*.csv`
prnGSPAHM(
  filter_by = exprs(distance <= .6),
  annot_cols = "ess_idx",
  annot_colnames = "Eset index",
  annot_rows = "ess_size",
  filename = show_some_redundancy.png,
)

# human terms only
prnGSPAHM(
  filter_num = exprs(distance <= .95),
  filter_sp = exprs(start_with_str("hs", term)),
  annot_cols = "ess_idx",
  annot_colnames = "Eset index",
  filename = show_more_connectivity.png,
)

# custom color palette
prnGSPAHM(
  annot_cols = c("ess_idx", "ess_size"),
  annot_colnames = c("Eset index", "Size"),
  filter_by = exprs(distance <= .95),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  filename = "custom_colors.png"
)



