# ===================================
# Peptide cleanup by CV
# ===================================
# percent CV
purgePep(pt_cv = .95)

# max CV
purgePep(max_cv = .5)

# `max_cv` then `pt_cv`
purgePep(max_cv = .5)
purgePep(pt_cv = .95)

# actually 90% CV 
purgePep(pt_cv = .95)
purgePep(pt_cv = .95)
