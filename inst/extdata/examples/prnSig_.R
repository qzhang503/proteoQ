# ===================================
# Significance tests
# ===================================

## !!!require the brief working example in `?load_expts`

# global option
scale_log2r <- TRUE

## significance tests and volcano plot visualization
# peptide
# (`Term`, `Term_2`, `Term_3` are user-defined, column keys in `expt_smry.xlsx`)
pepSig(
  impute_na = FALSE, 
  W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", 
                  "(W2.JHU.TMT2-W2.JHU.TMT1)", 
                  "(W2.PNNL.TMT2-W2.PNNL.TMT1)"], # batch effects
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"], # location effects
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

pepVol()

# protein
# (formulae matched to `pepSig` by default)
prnSig(impute_na = FALSE)
prnVol()

# note the incongruity in peptide and protein fold changes
# (no measures for peptides but for proteins)

#  sequence | sample_1 | sample_2 |    FC  
# --------------------------------------------
# prnX_pep1 |   2.15   |   NA     |    NA
# prnX_pep2 |    NA    |  1.05    |    NA
 
#  protein  | sample_1 | sample_2 |    FC  
# --------------------------------------------
#   prnX    |   2.15   |   1.05   |   2.05

## averaged batch effect
# (suggest always both `pepSig` and `prnSig` for consistent formulae)
pepSig(
  impute_na = FALSE,
  W2_loc_2 = ~ Term["(W2.BI.TMT2+W2.BI.TMT1)/2 - (W2.JHU.TMT2+W2.JHU.TMT1)/2"], 
)

prnSig(impute_na = FALSE)

pepVol()
prnVol()


## impute NA
# (suggested for models with random effects)
pepImp(m = 2, maxit = 2)
prnImp(m = 5, maxit = 5)

## random effects
# single
pepSig(
  impute_na = TRUE,
  W2_vs_W16_fix = ~ Term_3["W16-W2"], # fixed
  W2_vs_W16_mix = ~ Term_3["W16-W2"] + (1|TMT_Set), # one fixed and one random
)

prnSig(impute_na = TRUE)

pepVol()
prnVol()

# one to multiple
# (method `lm` for multiple random)
pepSig(
  impute_na = TRUE,
  method = lm,
  W2_vs_W16_fix = ~ Term_3["W16-W2"], # one fixed
  W2_vs_W16_mix = ~ Term_3["W16-W2"] + (1|TMT_Set), # one fixed and one random
  W2_vs_W16_mix_2 = ~ Term_3["W16-W2"] + (1|TMT_Set) + (1|Color), # one fixed and two randoms
)

prnSig(
  impute_na = TRUE,
  method = lm,
)

pepVol()
prnVol()

