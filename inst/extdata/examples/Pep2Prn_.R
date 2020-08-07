\donttest{
# ===================================
# Peptides to proteins
# ===================================

## !!!require the brief working example in `?load_expts`

# use unique peptides
Pep2Prn()

# include shared peptides
Pep2Prn(use_unique_pep = FALSE)  

# alignment of data by segments
Pep2Prn(cut_points = seq(4, 7, 0.5))
}
