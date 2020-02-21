# ===================================
# Volcano plots
# ===================================

## !!!require the brief working example in `?load_expts`

## global option
scale_log2r <- TRUE

## for all peptides or proteins
# all peptides
pepVol()

# all proteins
prnVol(
  show_labels = TRUE,
  xco = 1.2,
  yco = 0.01,
)

# kinases and prot_n_pep >= 2
prnVol(
  show_labels = TRUE,
  xco = 1.2,
  yco = 0.01,
  filter_prots_by_kin = exprs(kin_attr, prot_n_pep >= 2),
  filename = "kin_npep2.png"
)

# custom theme
my_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 24),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, size = 24),
    axis.title.x = element_text(colour = "black", size = 24),
    axis.title.y = element_text(colour="black", size = 24),
    plot.title = element_text(face = "bold", colour = "black", size = 14, hjust = .5, vjust = .5),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    strip.text.x = element_text(size = 16, colour = "black", angle = 0),
    strip.text.y = element_text(size = 16, colour = "black", angle = 90),
    
    legend.key = element_rect(colour = NA, fill = 'transparent'),
    legend.background = element_rect(colour = NA,  fill = "transparent"),
    legend.position = "none",
    legend.title = element_text(colour="black", size = 18),
    legend.text = element_text(colour="black", size = 18),
    legend.text.align = 0,
    legend.box = NULL
  )

prnVol(theme = my_theme, filename = my_theme.png)

## protein subgroups by gene sets
# prerequisite analysis of GSPA
prnGSPA(
  impute_na = FALSE,
  pval_cutoff = 5E-2,
  logFC_cutoff = log2(1.2),
  gspval_cutoff = 5E-2,
  gset_nms = c("go_sets"),
)

# filtered by proteins with two or more identifying peptides for visualization
gspaMap(
  gspval_cutoff = 5E-3,
  gslogFC_cutoff = log2(1.2),
  gset_nms = c("go_sets"),
  topn = 100, 
  show_sig = pVal,
  show_labels = TRUE,
  yco = 0.01,
  filter_prots_by_npep = exprs(prot_n_pep >= 2),
  # `filename`(s) will be automated, i.e., by gene-set names
)

# customized thresholds for the corresponding formulae in `pepSig` or `prnSig()`
# (may be suitable with the examplary differences in `W16_vs_W2` being much greater than `W2_bat`...)
gspaMap(
  fml_nms = c("W2_bat", "W2_loc", "W16_vs_W2"),
  gspval_cutoff = c(5E-2, 5E-2, 1E-10),
  gslogFC_cutoff = log2(1.2),
  topn = 100, 
  show_sig = pVal,
  show_labels = TRUE,
  yco = 0.05,
  filter_prots = exprs(prot_n_pep >= 2),
)

## gspaMap(...) maps secondary results of `[...]Protein_GSPA_{NZ}[_impNA].txt` 
#  from prnGSPA(...) onto a primary `df` of `Protein[_impNA]_pVal.txt` 
#  
#  see also ?prnGSPA for additional examples involving both `df` and `df2`, 
#  as well as `filter_` and `filter2_`


