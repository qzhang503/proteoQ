\donttest{
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
library(ggplot2)
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

pepHist(
  col_select = BI_1,
  theme = theme_dark(),
  filename = bi1_dark.png,
)


## direct uses of ggplot2
library(ggplot2)
res <- pepHist(filename = default.png)

# names(res)

p <- ggplot() +
  geom_histogram(data = res$raw, aes(x = value, y = ..count.., fill = Int_index),
                 color = "white", alpha = .8, binwidth = .05, size = .1) +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  labs(title = "", x = expression("Ratio (" * log[2] * ")"), y = expression("Frequency")) +
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1),
                     labels = as.character(seq(-2, 2, by = 1))) +
  scale_y_continuous(limits = NULL) + 
  facet_wrap(~ Sample_ID, ncol = 5, scales = "fixed") # + 
  # my_histo_theme

p <- p + 
  geom_line(data = res$fitted, mapping = aes(x = x, y = value, colour = variable), size = .2) +
  scale_colour_manual(values = c("gray", "gray", "gray", "black"), name = "Gaussian",
                      breaks = c(c("G1", "G2", "G3"), paste(c("G1", "G2", "G3"), collapse = " + ")),
                      labels = c("G1", "G2", "G3", "G1 + G2 + G3"))

p <- p + geom_vline(xintercept = 0, size = .25, linetype = "dashed")

ggsave(file.path(dat_dir, "Peptide\\Histogram\\my_ggplot2.png"), 
       width = 22, height = 48, limitsize = FALSE)

\dontrun{
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
}
}
