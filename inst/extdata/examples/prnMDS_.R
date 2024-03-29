\donttest{
# ===================================
# MDS
# ===================================

## !!!require the brief working example in `?load_expts`

# global option
scale_log2r <- TRUE

## peptides
# all samples
pepMDS(
  col_select = Select, 
  filter_peps_by = exprs(pep_n_psm >= 10),
  show_ids = FALSE, 
  filename = "peps_rowfil.png",
)

# selected samples
pepMDS(
  col_select = BI, 
  col_shape = Shape,   
  col_color = Alpha, 
  filter_peps_by = exprs(pep_n_psm >= 10),
  show_ids = FALSE, 
  filename = "peps_rowfil_colsel.png",
)

# column `Alpha` will be used at the default of
# `col_alpha = NULL`;
# To bypass the aesthetics under column `Alpha`, 
# use `col_alpha = NA`
# (the same applies to other aesthetics, and PCA and LDA)
pepMDS(
  col_select = Select, 
  col_alpha = NA, 
  filter_peps_by = exprs(pep_n_psm >= 10),
  show_ids = FALSE, 
  filename = "peps_rowfil_no_alpha.png",
)


## proteins
prnMDS(
  col_color = Color,
  col_shape = Shape,
  show_ids = FALSE,
  filter_peps_by = exprs(prot_n_pep >= 5),
  filename = "prns_rowfil.png",
)

# custom palette
prnMDS(
  col_shape = Shape,
  color_brewer = Set1,
  show_ids = FALSE,
  filename = "my_palette.png",
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
prnMDS(
  col_color = Color,
  col_shape = Shape,
  show_ids = FALSE,
  filter_peps_by = exprs(prot_n_pep >= 5),
  filter_by = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
  filename = pvalcutoff.png, 
)

# analogous peptides
pepMDS(
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

prnMDS(
  impute_na = TRUE,
  col_color = Color,
  col_shape = Shape,
  show_ids = FALSE,
  filter_peps_by = exprs(prot_n_pep >= 5),
  filter_by = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
  filename = filpvals_impna.png, 
)

# analogous peptides
pepMDS(
  impute_na = TRUE,
  col_color = Color,
  col_shape = Shape,
  show_ids = FALSE,
  filter_peps_by = exprs(prot_n_pep >= 5),
  filter_by = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
  filename = filpvals_impna.png,
)

## show ellipses
prnMDS(
  show_ellipses = TRUE,
  col_group = Shape, 
  show_ids = FALSE,
  filename = ellipses_by_whims.png,
)

prnMDS(
  show_ellipses = TRUE,
  col_group = Color, 
  show_ids = FALSE,
  filename = ellipses_by_labs.png,
)

## a higher dimension
pepMDS(
  show_ids = FALSE,
  k = 5, 
  dimension = 3,
  filename = d3.pdf,
)

prnMDS(
  show_ids = TRUE,
  k = 4, 
  dimension = 3,
  filename = d3.png,
)

# show ellipses
# (column `expt_smry.xlsx::Color` codes `labs`.)
prnMDS(
  show_ids = FALSE,
  show_ellipses = TRUE,
  col_group = Color, 
  k = 4, 
  dimension = 3,
  filename = d3_labs.png,
)

# (column `expt_smry.xlsx::Shape` codes `WHIMs`.)
prnMDS(
  show_ids = FALSE,
  show_ellipses = TRUE,
  col_group = Shape, 
  k = 4, 
  dimension = 3,
  filename = d3_whims.png,
)


# toy example of finding samples(s) that are 
# most different in large fold changes;
prnMDS(
  show_ids = TRUE, 
  dist_co = log2(4),
  filename = where_are_the_large_diffs.png,
)


## custom theme
library(ggplot2)
my_mds_theme <- theme_bw() + theme(
  axis.text.x  = element_text(angle=0, vjust=0.5, size=16),
  axis.text.y  = element_text(angle=0, vjust=0.5, size=16),
  axis.title.x = element_text(colour="black", size=18),
  axis.title.y = element_text(colour="black", size=18),
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

pepMDS(
  impute_na = FALSE,
  col_color = Color,
  col_shape = Shape,
  show_ids = FALSE,
  filter_peps_by = exprs(prot_n_pep >= 5),
  filter_by = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
  theme = my_mds_theme,
  filename = my_theme.png,
)

## direct uses of ggplot2
library(ggplot2)
res <- prnMDS(filename = foo.png)

p_fil <- ggplot(res, aes(x = Coordinate.1, y = Coordinate.2)) +
  geom_point(aes(colour = Color, shape = Shape, alpha = Alpha), size = 4, stroke = 0.02) + 
  scale_alpha_manual(values = c(.5, .9)) + 
  stat_ellipse(aes(fill = Shape), geom = "polygon", alpha = .4) + 
  guides(fill = FALSE) + 
  labs(title = "", x = "Coordinate 1", y = "Coordinate 2") +
  coord_fixed() 

ggsave(file.path(dat_dir, "Protein/MDS/my_ggplot2_fil.png"))

\dontrun{
prnMDS(
  col_color = "column_key_not_existed",
  col_shape = "another_missing_column_key"
)  
}
}
