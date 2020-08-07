\donttest{
# ===================================
# Merge peptide data
# ===================================

## !!!require the brief working example in `?load_expts`

# everything included
mergePep()

# row filtrations against column keys in `TMTset1_LCMSinj1_Peptide_N.txt`...
# (with Mascot column keys)
mergePep(
  filter_peps_by_sp = exprs(species == "human", pep_len <= 50),
)

# (with MaxQuant column keys)
mergePep(
  filter_peps_by_sp = exprs(species == "human", `Missed cleavages` <= 3),
)

# (with Spectrum Mill column keys)
mergePep(
  filter_peps_by_sp = exprs(species == "human", pep_miss <= 3),
)


# alignment of data by segments
mergePep(cut_points = seq(4, 7, .5))
}
