\donttest{
# ===================================
# PhosphoSitePlus (PSP)
# ===================================

## !!!require the brief working example in `?load_expts`

library(proteoQ)

# expression data for kinases and their substrates
# (need to first download the ".txt" from PSP)
anal_KinSub("~/proteoQ/dbs/psp/Kinase_Substrate_Dataset.txt")

# `human` only ('unknown' species will be removed)
anal_KinSub(
  filter_by_sp = exprs(species == "human"),
  filter_prots_by = exprs(prot_n_pep >= 2),
  filename = human_2peps.txt,
)

# proteins
anal_KinSub(type = protein)
}

\dontrun{
  # peptides: CDK1 substrates
  anal_KinSub(
    filter_by_gene = exprs(species == "human", gene == "CDK1"),
    filename = hu_CDK1.txt,
  )
  
  # proteins: CDK1 substrates
  anal_KinSub(
    type = protein,
    filter_by_gene = exprs(species == "human", gene == "CDK1"),
    filename = hu_CDK1.txt,
  )
}
