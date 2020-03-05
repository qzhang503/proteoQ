\donttest{
# ===================================
# PSM cleanup by CV
# ===================================

## !!!require the brief working example in `?load_expts`

## additional examples
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

## custom theme
purgePSM(dpi = 300, theme = ggplot2::theme_dark())  
}
