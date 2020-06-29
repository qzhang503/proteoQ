\donttest{
# ===================================
# PCA
# ===================================

## !!!require the brief working example in `?load_expts`

## global option
scale_log2r <- TRUE

# peptides, all samples
pepPCA(
  col_select = Select, 
  filter_peps_by = exprs(pep_n_psm >= 3),
  show_ids = FALSE, 
  filename = "peps_rowfil.png",
)

# peptides, samples under column `BI`
pepPCA(
  col_select = BI, 
  col_shape = Shape,   
  col_color = Alpha, 
  filter_peps_by = exprs(pep_n_psm >= 10),
  show_ids = FALSE, 
  filename = "peps_rowfil_colsel.png",
)

# proteins
prnPCA(
  col_color = Color,
  col_shape = Shape,
  show_ids = FALSE,
  filter_peps_by = exprs(prot_n_pep >= 5),
  filename = "prns_rowfil.png",
)

# subset by mean deviation values
# deviations to means may not be symmetric;
prnPCA(
  col_select = Select, 
  filter_peps_by = exprs(prot_mean_z >= -.25, prot_mean_z <= .3),
  show_ids = FALSE, 
  filename = "subset_by_mean_dev.png",
)

# proteins, custom palette
prnPCA(
  col_shape = Shape,
  color_brewer = Set1,
  show_ids = FALSE,
  filename = "my_palette.png",
)

# proteins, by features
prnPCA(
  type = feats,
  scale_log2r = TRUE,
  filename = "by_feats.png",
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
prnPCA(
  col_color = Color,
  col_shape = Shape,
  show_ids = FALSE,
  filter_peps_by = exprs(prot_n_pep >= 5),
  filter_by = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
  filename = pvalcutoff.png, 
)

# analogous peptides
prnPCA(
  col_color = Color,
  col_shape = Shape,
  show_ids = FALSE,
  filter_peps_by = exprs(prot_n_pep >= 5),
  filter_by = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
  filename = pvalcutoff.png, 
)

## additional row filtration by pVals (proteins, impute_na = TRUE)
# if not yet, run prerequisitive NA imputation
pepImp(m = 2, maxit = 2)
prnImp(m = 5, maxit = 5)

# if not yet, run prerequisitive significance tests at `impute_na = TRUE`
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

prnPCA(
  impute_na = TRUE,
  col_color = Color,
  col_shape = Shape,
  show_ids = FALSE,
  filter_peps_by = exprs(prot_n_pep >= 5),
  filter_by = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
  filename = filpvals_impna.png, 
)

# analogous peptides
pepPCA(
  impute_na = TRUE,
  col_color = Color,
  col_shape = Shape,
  show_ids = FALSE,
  filter_peps_by = exprs(prot_n_pep >= 5),
  filter_by = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
  filename = filpvals_impna.png,
)

## a higher dimension
pepPCA(
  show_ids = FALSE,
  rank. = 5, 
  dimension = 3,
  filename = d3.pdf,
)

prnPCA(
  show_ids = TRUE,
  rank. = 4, 
  dimension = 3,
  filename = d3.png,
)

prnPCA(
  type = feats,
  rank. = 4, 
  dimension = 3,
  filename = feat_d3.png,
)

# show ellipses
# (column `expt_smry.xlsx::Color` codes `labs`.)
prnPCA(
  show_ids = FALSE,
  show_ellipses = TRUE,
  col_group = Color, 
  rank. = 4, 
  dimension = 3,
  filename = d3_labs.png,
)

# (column `expt_smry.xlsx::Shape` codes `WHIMs`.)
prnPCA(
  show_ids = FALSE,
  show_ellipses = TRUE,
  col_group = Shape, 
  rank. = 4, 
  dimension = 3,
  filename = d3_whims.png,
)

## custom theme
library(ggplot2)
my_theme <- theme_bw() + theme(
  axis.text.x  = element_text(angle=0, vjust=0.5, size=20),
  axis.text.y  = element_text(angle=0, vjust=0.5, size=20),
  axis.title.x = element_text(colour="black", size=20),
  axis.title.y = element_text(colour="black", size=20),
  plot.title = element_text(face="bold", colour="black", size=20, hjust=0.5, vjust=0.5),
  
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  
  legend.key = element_rect(colour = NA, fill = 'transparent'),
  legend.background = element_rect(colour = NA,  fill = "transparent"),
  legend.title = element_blank(),
  legend.text = element_text(colour="black", size=14),
  legend.text.align = 0,
  legend.box = NULL
)

pepPCA(
  impute_na = TRUE,
  col_color = Color,
  col_shape = Shape,
  show_ids = FALSE,
  filter_peps_by = exprs(prot_n_pep >= 5),
  filter_by = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
  theme = my_theme, 
  filename = my_theme.png,
)

## direct uses of ggplot2
library(ggplot2)
res <- prnPCA(filename = foo.png)

# names(res)

p_fil <- ggplot(res$pca, aes(PC1, PC2)) +
  geom_point(aes(colour = Color, shape = Shape, alpha = Alpha), size = 4, stroke = 0.02) + 
  scale_alpha_manual(values = c(.5, .9)) + 
  stat_ellipse(aes(fill = Shape), geom = "polygon", alpha = .4) + 
  guides(fill = FALSE) + 
  labs(title = "", 
       x = paste0("PC1 (", res$prop_var[1], ")"), 
       y = paste0("PC2 (", res$prop_var[2], ")")) +
  coord_fixed() 

ggsave(file.path(dat_dir, "Protein/PCA/my_ggplot2_fil.png"))

\dontrun{
prnPCA(
  col_color = "column_key_not_existed",
  col_shape = "another_missing_column_key"
)  
}
}

