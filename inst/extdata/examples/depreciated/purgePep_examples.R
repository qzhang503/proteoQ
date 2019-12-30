# ===================================
# Peptide cleanup by CV
# ===================================
# by percent CV
purgePep(pt_cv = .95)

# by max CV
purgePep(max_cv = .5)

# by `max_cv` then by `pt_cv`
purgePep(max_cv = .5)
purgePep(pt_cv = .95)

# actually 90% CV 
purgePep(pt_cv = .95)
purgePep(pt_cv = .95)
