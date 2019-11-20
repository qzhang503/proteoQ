# ===================================
# String DB
# ===================================
getStringDB(
  db_path = "~\\proteoQ\\dbs\\string",
  score_cutoff = .9,
  adjP = FALSE,
  filter_by_sp = exprs(species == "human"),
  filter_by_npep = exprs(prot_n_pep >= 2),
)

