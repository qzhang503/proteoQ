\donttest{
# ===================================
# String DB
# ===================================

## !!!require the brief working example in `?load_expts`

library(proteoQ)

# `human` and `mouse` STRING using default urls;
prepString(human)
prepString(mouse)

# custom `human` and `mouse` STRING
prepString(
  species = does_not_matter_at_custom_urls,
  links_url = "https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz",
  aliases_url = "https://stringdb-static.org/download/protein.aliases.v11.0/9606.protein.aliases.v11.0.txt.gz",
  info_url = "https://stringdb-static.org/download/protein.info.v11.0/9606.protein.info.v11.0.txt.gz", 
  filename = my_hs.rds,
)

prepString(
  # species = this_mouse,
  links_url = "https://stringdb-static.org/download/protein.links.full.v11.0/10090.protein.links.full.v11.0.txt.gz",
  aliases_url = "https://stringdb-static.org/download/protein.aliases.v11.0/10090.protein.aliases.v11.0.txt.gz",
  info_url = "https://stringdb-static.org/download/protein.info.v11.0/10090.protein.info.v11.0.txt.gz", 
  filename = my_mm.rds,
)
}

\dontrun{
identical(
  readRDS(file.path("~\\proteoQ\\dbs\\string\\string_hs.rds")), 
  readRDS(file.path("~\\proteoQ\\dbs\\string\\my_hs.rds"))
)
}

\donttest{
# analysis: both `human` and `mouse`
anal_prnString(
  db_nms = c("~\\proteoQ\\dbs\\string\\string_hs.rds",
             "~\\proteoQ\\dbs\\string\\string_mm.rds"),
  score_cutoff = .9,
  filter_prots_by = exprs(prot_n_pep >= 2),
)


# `human` only ('unknown' species will be removed)
# OK to include both `string_hs.rds` and `string_mm.rds`
anal_prnString(
  db_nms = c("~\\proteoQ\\dbs\\string\\string_hs.rds",
             "~\\proteoQ\\dbs\\string\\string_mm.rds"),
  score_cutoff = .9,
  filter_by_sp = exprs(species == "human"),
  filter_prots_by = exprs(prot_n_pep >= 2),
  filename = human.tsv,
)

# `mouse` only
anal_prnString(
  db_nms = c("~\\proteoQ\\dbs\\string\\string_hs.rds",
             "~\\proteoQ\\dbs\\string\\string_mm.rds"),
  score_cutoff = .9,
  filter_by_sp = exprs(species == "mouse"),
  filter_prots_by = exprs(prot_n_pep >= 2),
  filename = mouse.tsv,
)


# additional filtration by `pVals` and `log2FC`; 
# `W16_vs_W2.pVal (W16-W2)` is a column key in `Protein_pVals.txt`
anal_prnString(
  db_nms = "~\\proteoQ\\dbs\\string\\string_hs.rds",
  score_cutoff = .9,
  filter_by_sp = exprs(species == "human", 
                       `W16_vs_W2.pVal (W16-W2)` <= 1E-6,
                       abs(`W16_vs_W2.log2Ratio (W16-W2)`) >= 1.2),
  filter_prots_by = exprs(prot_n_pep >= 2),
  filename = human_sigs.tsv,
)

anal_prnString(
  db_nms = "~\\proteoQ\\dbs\\string\\string_mm.rds",
  score_cutoff = .9,
  filter_by_sp = exprs(species == "mouse", 
                       `W16_vs_W2.pVal (W16-W2)` <= 1E-6,
                       abs(`W16_vs_W2.log2Ratio (W16-W2)`) >= 1.2),
  filter_prots_by = exprs(prot_n_pep >= 2),
  filename = mouse_sigs.tsv,
)

# can incorporate `prepString` into `anal_prnString`
anal_prnString(
  db_nms = c(prepString(human),
             prepString(mouse)),
  score_cutoff = .9,
  filter_prots_by = exprs(prot_n_pep >= 2),
  filename = one_pot.tsv,
)
}