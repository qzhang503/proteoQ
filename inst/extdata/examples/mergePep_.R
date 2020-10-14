\donttest{
# ===================================
# Merge peptide data
# ===================================

## !!!require the brief working example in `?load_expts`

# everything included
mergePep()

# row filtrations against column keys in `TMTset1_LCMSinj1_Peptide_N.txt`...
mergePep(
  filter_peps_by_sp = exprs(species == "human", pep_len <= 50),
)

# alignment of data by segments
mergePep(cut_points = c(mean_lint = seq(4, 7, .5)))

# alignment of data by empirical protein abundance
# `10^prot_icover - 1` comparable to emPAI
mergePep(cut_points = c(prot_icover = seq(0, 1, .25)))
}
