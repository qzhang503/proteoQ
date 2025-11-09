\donttest{
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
  xco = 1.2,
  yco = 0.01,
)

# hide `xco` and/or `yco` lines
# (xco = 0 -> log2(xco) = - Inf)
prnVol(
  xco = 0,
  yco = Inf,
  filename = no_xylines.png,
)

# shows vertical center line at log2(1)
# (xco = 1 -> log2(xco) = 0)
prnVol(
  xco = 1,
  yco = Inf,
  filename = no_xylines.png,
)

# kinases and prot_n_pep >= 2
prnVol(
  xco = 1.2,
  yco = 0.01,
  filter_prots_by_kin = exprs(kin_attr, prot_n_pep >= 2),
  filename = "kin_npep2.png"
)

# selected formula and/or customization
prnVol(
  fml_nms = "W2_bat",
  xmin = -5,
  xmax = 5, 
  ymin = 0, 
  ymax = 30,
  x_label = "Ratio ("*log[2]*")",
  y_label = "pVal ("*-log[10]*")", 
  height = 6,
  width = 6*2.7,
  filename = custom.png,
)

# custom theme
library(ggplot2)
my_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 24),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, size = 24),
    axis.title.x = element_text(colour = "black", size = 24),
    axis.title.y = element_text(colour="black", size = 24),
    plot.title = element_text(face = "bold", colour = "black", size = 14, 
                              hjust = .5, vjust = .5),
    
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

# custom plot
# ("W2_bat" etc. are contrast names in `pepSig`)
prnVol(fml_nms = c("W2_bat", "W2_loc"), filename = foo.png)

res <- readRDS(file.path(dat_dir, "Protein/Volcano/W2_bat/foo.rds"))
# names(res)

p <- ggplot() +
  geom_point(data = res$data, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
             size = 3, colour = "#f0f0f0", shape = 20, alpha = .5) +
  geom_point(data = res$greater, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
             size = 3, color = res$palette[2], shape = 20, alpha = .8) +
  geom_point(data = res$less, mapping = aes(x = log2Ratio, y = -log10(pVal)), 
             size = 3, color = res$palette[1], shape = 20, alpha = .8) +
  geom_hline(yintercept = -log10(res$yco), linetype = "longdash", linewidth = .5) +
  geom_vline(xintercept = -log2(res$xco), linetype = "longdash", linewidth = .5) +
  geom_vline(xintercept = log2(res$xco), linetype = "longdash", linewidth = .5) +
  scale_x_continuous(limits = c(res$xmin, res$xmax)) +
  scale_y_continuous(limits = c(res$ymin, res$ymax)) +
  labs(title = res$title, x = res$x_label, y = res$y_label) +
  res$theme

p <- p + geom_text(data = res$topns, 
                   mapping = aes(x = log2Ratio, 
                                 y = -log10(pVal), 
                                 label = Index, 
                                 color = Index),
                   size = 3, 
                   alpha = .5, 
                   hjust = 0, 
                   nudge_x = 0.05, 
                   vjust = 0, 
                   nudge_y = 0.05, 
                   na.rm = TRUE)

p <- p + facet_wrap(~ Contrast, nrow = 1, labeller = label_value)

p <- p + geom_table(data = res$topn_labels, aes(table = gene), 
                    x = -res$xmax*.85, y = res$ymax/2)

# Highlight
prnVol(
  highlights = rlang::exprs(gene %in% c("ACTB", "GAPDH")), 
  filename = highlights.png
)


## protein subgroups by gene sets
# prerequisite analysis of GSPA
prnGSPA(
  impute_na = FALSE,
  pval_cutoff = 1E-2, # protein pVal threshold
  logFC_cutoff = log2(1.1), # protein log2FC threshold
  gspval_cutoff = 5E-2, # gene-set pVal threshold
  gslogFC_cutoff = log2(1.2), # gene-set log2FC threshold
  gset_nms = c("go_sets"),
)

# mapping gene sets to volcano-plot visualization
# (1) forced lines of `pval_cutoff` and `logFC_cutoff`  
#   according to the corresponding `prnGSPA` in red; 
# (2) optional lines of `xco` and `yco` in grey
gspaMap(
  impute_na = FALSE,
  topn_gsets = 20, 
  show_sig = pVal, 
)

# disable the lines of `xco` and `yco`, 
gspaMap(
  impute_na = FALSE,
  topn_gsets = 20, 
  show_sig = pVal, 
  xco = 0, 
  yco = Inf, 
)

# customized thresholds for visualization
gspaMap(
  fml_nms = c("W2_bat", "W2_loc", "W16_vs_W2"),
  gspval_cutoff = c(5E-2, 5E-2, 1E-10),
  gslogFC_cutoff = log2(1.2),
  topn_gsets = 20, 
  topn_labels = 0,
  show_sig = pVal,
  xco = 0, 
  yco = Inf, 
)

## gspaMap(...) maps secondary results of `[...]Protein_GSPA_{NZ}[_impNA].txt` 
#  from prnGSPA(...) onto a primary `df` of `Protein[_impNA]_pVal.txt` 
#  
#  see also ?prnGSPA for additional examples involving both `df` and `df2`, 
#  as well as `filter_` and `filter2_`
}
