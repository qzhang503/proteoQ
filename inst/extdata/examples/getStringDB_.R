# ===================================
# String DB
# ===================================

## !!!require the brief working example in `?load_expts`

string_dir <- "~\\proteoQ\\dbs\\string"
dir.create(string_dir, recursive = TRUE, showWarnings = FALSE)

# download
dl_stringdbs(
  species = human,
  db_path = string_dir,
)

# protein-protein interaction network
getStringDB(
  db_path = string_dir,
  score_cutoff = .9,
  adjP = FALSE,
  filter_by_sp = exprs(species == "human"),
  filter_prots_by = exprs(prot_n_pep >= 2),
)

