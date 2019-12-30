# ===================================
# PSM cleanup by CV
# ===================================
# visualization only
purgePSM()

# by percent CV
purgePSM(pt_cv = .95)

# by max CV
purgePSM(max_cv = .5)

# by `max_cv` then by `pt_cv`
purgePSM(max_cv = .5)
purgePSM(pt_cv = .95)

# actually 90% CV 
purgePSM(pt_cv = .95)
purgePSM(pt_cv = .95)

