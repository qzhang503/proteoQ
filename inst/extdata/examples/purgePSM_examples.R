# ===================================
# PSM cleanup by CV
# ===================================
# percent CV
purgePSM(pt_cv = .95)

# max CV
purgePSM(max_cv = .5)

# `max_cv` then `pt_cv`
purgePSM(max_cv = .5)
purgePSM(pt_cv = .95)

# actually 90% CV 
purgePSM(pt_cv = .95)
purgePSM(pt_cv = .95)
