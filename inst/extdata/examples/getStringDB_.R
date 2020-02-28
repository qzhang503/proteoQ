\donttest{
# ===================================
# String DB
# ===================================

## !!!require the brief working example in `?load_expts`

string_dir <- "~\\proteoQ\\dbs\\string"
dir.create(string_dir, recursive = TRUE, showWarnings = FALSE)

## download
dl_stringdbs(
  species = human,
  db_path = string_dir,
)

## ppi and expression
# human and mouse and probably 'unknown' species
anal_prnString(
  db_path = string_dir,
  score_cutoff = .9,
  filter_prots_by = exprs(prot_n_pep >= 2),
)

# cleaner outputs by specifying species
# ('unknown' will be removed)
anal_prnString(
  db_path = string_dir,
  score_cutoff = .9,
  filter_by_sp = exprs(species == "human"),
  filter_prots_by = exprs(prot_n_pep >= 2),
)

anal_prnString(
  db_path = string_dir,
  score_cutoff = .9,
  filter_by_sp = exprs(species == "mouse"),
  filter_prots_by = exprs(prot_n_pep >= 2),
)

# custom filename with outputs in: 
# prot_npep2 W16_vs_W2.pVal (W16-W2)_human_expr.tsv
# prot_npep2 W16_vs_W2.pVal (W16-W2)_human_ppi.tsv
# prot_npep2 W16_vs_W2.pVal (W16-W2)_mouse_expr.tsv
# prot_npep2 W16_vs_W2.pVal (W16-W2)_mouse_ppi.tsv
anal_prnString(
  db_path = string_dir,
  score_cutoff = .9,
  complete_cases = TRUE, 
  filter_by = exprs(species %in% c("human", "mouse"), 
                    prot_n_pep >= 2, 
                    `W16_vs_W2.pVal (W16-W2)` <= 1e-6),
  filename = "prot_npep2 W16_vs_W2.pVal (W16-W2).txt",
)
}
