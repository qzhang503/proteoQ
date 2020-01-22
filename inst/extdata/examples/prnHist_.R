# ===================================
# Histogram
# ===================================

## !!!require the brief working example in `?load_expts`

## examplary `MGKernel` alignment
standPep(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
)

standPrn(
  method_align = MGKernel, 
  n_comp = 2, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
)

## (1) effects of data scaling
# peptide without log2FC scaling
pepHist(scale_log2r = FALSE)

# with scaling
pepHist(scale_log2r = TRUE)

## (2) sample column selection
# sample IDs indicated under column `Select` in `expt_smry.xlsx`
pepHist(col_select = Select, filename = colsel.png)

# protein data for samples under column `W2` in `expt_smry.xlsx`
prnHist(col_select = W2, filename = w2.png)

## (3) row filtration of data
# exclude oxidized methione or deamidated asparagine
pepHist(
  # filter_by = exprs(!grepl("[mn]", pep_seq_mod)),
  filter_by = exprs(not_contain_chars_in("mn", pep_seq_mod)),
  filename = "no_mn.png",
)

# phosphopeptide subset (error message if no matches)
pepHist(
  filter_peps = exprs(contain_chars_in("sty", pep_seq_mod)), 
  scale_y = FALSE, 
  filename = phospho.png,
)

# or use `grepl` directly
pepHist(
  filter_by = exprs(grepl("[sty]", pep_seq_mod)),
  filename = same_phospho.png,
)

## (4) between lead and lag
# leading profiles
pepHist(
  filename = lead.png,
)

# lagging profiles at
#   (1) n_psm >= 10
#   (2) and no methionine oxidation or asparagine deamidation
pepHist(
  filter_peps_by_npsm = exprs(pep_n_psm >= 10),
  filter_peps_by_mn = exprs(not_contain_chars_in("mn", pep_seq_mod)),
  filename = lag.png,
)

## custom theme
my_histo_theme <- theme_bw() + theme(
  axis.text.x  = element_text(angle=0, vjust=0.5, size=18),
  axis.ticks.x  = element_blank(), # x-axis ticks
  axis.text.y  = element_text(angle=0, vjust=0.5, size=18),
  axis.title.x = element_text(colour="black", size=24),
  axis.title.y = element_text(colour="black", size=24),
  plot.title = element_text(colour="black", size=24, hjust=.5, vjust=.5),
  
  strip.text.x = element_text(size = 18, colour = "black", angle = 0),
  strip.text.y = element_text(size = 18, colour = "black", angle = 90),
  
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  
  legend.key = element_rect(colour = NA, fill = 'transparent'),
  legend.background = element_rect(colour = NA,  fill = "transparent"),
  legend.title = element_blank(),
  legend.text = element_text(colour="black", size=18),
  legend.text.align = 0,
  legend.box = NULL
)

pepHist(
  theme = my_histo_theme,
  filename = my_theme.png,
)

## Not run: 
# sample selection
pepHist(
  col_select = "a_column_key_not_in_`expt_smry.xlsx`",
)

# data filtration
pepHist(
  filter_by = exprs(!grepl("[m]", a_column_key_not_in_data_table)),
)

prnHist(
  lhs_not_start_with_filter_ = exprs(n_psm >= 5),
)

## End(Not run)