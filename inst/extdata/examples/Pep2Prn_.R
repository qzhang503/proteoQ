\donttest{
# ===================================
# Peptides to proteins
# ===================================

## !!!Require the brief working example in `?load_expts`

# Use unique peptides
Pep2Prn()

# Include shared peptides
Pep2Prn(use_unique_pep = FALSE)  

# Alignment of data by segments
Pep2Prn(cut_points = seq(4, 7, 0.5))
}

