\donttest{
# ===================================
# Peptide cleanup by CV
# ===================================

## !!!require the brief working example in `?load_expts`

## additional examples
# visualization only
purgePep()

# visualization of non-void entries under column `BI_1` in `expt_smry.xlsx`
purgePep(
  col_select = BI_1,
  width = 8,
  height = 8,
  filename = bi_1.png,
)

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

## custom theme
purgePep(dpi = 300, theme = ggplot2::theme_get(), filename = my_theme.png)
}

