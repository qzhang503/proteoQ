\donttest{
# ===================================
# Trend analysis
# ===================================

## !!!require the brief working example in `?load_expts`

## global option
scale_log2r <- TRUE


# ===================================
# Analysis
# ===================================
## base (proteins, with sample order supervision)
anal_prnTrend(
  impute_na = FALSE,
  col_order = Order,
  n_clust = c(5:6), 
)

## against selected samples
anal_prnTrend(
  col_select = BI,
  impute_na = FALSE,
  col_order = Order,
  n_clust = c(5:6), 
  filename = sel.txt,
)

## row filtration (proteins)
anal_prnTrend(
  impute_na = FALSE,
  col_order = Order,
  n_clust = c(5:6), 
  filter_prots_by = exprs(prot_n_pep >= 3),
)

## manual m degree of fuzziness (proteins)
anal_prnTrend(
  impute_na = FALSE,
  col_order = Order,
  n_clust = c(5:6), 
  filter_prots = exprs(prot_n_pep >= 3),
  m = 1.5,
  filename = my_m.txt,
)

## additional row filtration by pVals (proteins, impute_na = FALSE)
# if not yet, run prerequisitive significance tests at `impute_na = FALSE`
pepSig(
  impute_na = FALSE, 
  W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", 
                  "(W2.JHU.TMT2-W2.JHU.TMT1)", 
                  "(W2.PNNL.TMT2-W2.PNNL.TMT1)"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = FALSE)

# (`W16_vs_W2.pVal (W16-W2)` now a column key)
anal_prnTrend(
  impute_na = FALSE,
  col_order = Order,
  n_clust = c(5:6), 
  filter_prots_by_npep = exprs(prot_n_pep >= 3), 
  filter_prots_by_pval = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
)

## additional row filtration by pVals (impute_na = TRUE)
# if not yet, run prerequisitive NA imputation and corresponding 
# significance tests at `impute_na = TRUE`
pepImp(m = 2, maxit = 2)
prnImp(m = 5, maxit = 5)

pepSig(
  impute_na = TRUE, 
  W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", 
                  "(W2.JHU.TMT2-W2.JHU.TMT1)", 
                  "(W2.PNNL.TMT2-W2.PNNL.TMT1)"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = TRUE)

anal_prnTrend(
  impute_na = TRUE,
  col_order = Order,
  n_clust = c(5:6), 
  filter_prots_by_npep = exprs(prot_n_pep >= 3), 
  filter_prots_by_pval = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
)


# ===================================
# Visualization
# ===================================
## base (proteins, no NA imputation) 
plot_prnTrend(
  col_order = Order, 
)

# at specific cluster ID(s)
# (`cluster` is a column key in `Protein_Trend_[...].txt`)
plot_prnTrend(
  impute_na = FALSE, 
  col_order = Order,
  filter2_by_clusters = exprs(cluster == 5),
  width = 8, 
  height = 10,
  filename = cl5.png,
)

# manual selection of secondary input data file(s)
# may be used for optimizing individual plots
plot_prnTrend(
  df2 = c("Protein_Trend_Z_nclust5.txt"),
  col_order = Order, 
  filename = n5.png,
)

# manual secondary input(s) at specific rank(s)
plot_prnTrend(
  df2 = c("Protein_Trend_Z_nclust5.txt"),
  impute_na = FALSE, 
  col_order = Order,
  filter2_by_clusters = exprs(cluster == 5),
  width = 8, 
  height = 10,
  filename = n5_cl5.png,
)

## NA imputation
# also save as pdf
plot_prnTrend(
  impute_na = TRUE,
  col_order = Order,
  filename = my.pdf,
)

## against selected samples
plot_prnTrend(
  col_order = Order, 
  col_select = BI,
  filename = bi.png,
)

## custom theme
library(ggplot2)
my_trend_theme <- theme_bw() + theme(
  axis.text.x  = element_text(angle=60, vjust=0.5, size=24),
  axis.ticks.x  = element_blank(), 
  axis.text.y  = element_text(angle=0, vjust=0.5, size=24),
  axis.title.x = element_text(colour="black", size=24),
  axis.title.y = element_text(colour="black", size=24),
  plot.title = element_text(face="bold", colour="black",
                            size=20, hjust=.5, vjust=.5),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.background = element_rect(fill = '#0868ac', colour = 'red'),
  
  strip.text.x = element_text(size = 24, colour = "black", angle = 0),
  strip.text.y = element_text(size = 24, colour = "black", angle = 90),
  
  plot.margin = unit(c(5.5, 55, 5.5, 5.5), "points"), 
  
  legend.key = element_rect(colour = NA, fill = 'transparent'),
  legend.background = element_rect(colour = NA,  fill = "transparent"),
  legend.position = "none",
  legend.title = element_text(colour="black", size=18),
  legend.text = element_text(colour="black", size=18),
  legend.text.align = 0,
  legend.box = NULL
)

plot_prnTrend(
  col_order = Order, 
  col_select = BI,
  theme = my_trend_theme,
  filename = my_theme.png,
)

## Cytoscape visualization
# (Make sure that Cytoscape is open.)
# Human
cluego(
  df2 = Protein_Trend_Z_nclust5.txt, 
  species = c(human = "Homo Sapiens"), 
  n_clust = c(3, 5)
)

# Mouse
cluego(
  df2 = Protein_Trend_Z_nclust5.txt, 
  species = c(mouse = "Mus Musculus"), 
  n_clust = c(3:4)
)

}
