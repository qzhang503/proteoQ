\donttest{
# ===================================
# Heat map
# ===================================

## !!!require the brief working example in `?load_expts`

## global option
scale_log2r <- TRUE

## proteins
# row clustering
prnHM(
  xmin = -1,
  xmax = 1,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = TRUE,
  cutree_rows = 10,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_row = 3,
  cellwidth = 14,
  width = 18,
  height = 12,
  filter_sp = exprs(species == "human", prot_n_pep >= 2),
  filename = "huprns_npep2.png",
)

# rows ordered by kinase classes then by gene names
# (error if `normPSM(annot_kinases = FALSE, ...)`)
prnHM(
  xmin = -1,
  xmax = 1,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = FALSE,
  annot_rows = c("kin_class"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 2,
  cellheight = 2,
  cellwidth = 14,
  width = 22,
  height = 22,
  filter_kin = exprs(kin_attr, species == "human"),
  arrange_kin = exprs(kin_order, gene),
  filename = "hukins_rows_by_class.png",
)

# `cutree_rows` ignored at `cluster_rows = FALSE`
prnHM(
  scale_log2r = TRUE,
  annot_cols = c("Group"),
  cluster_rows = FALSE,
  clustering_distance_rows  = "maximum",
  cutree_rows = 6,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_row = 3,
  cellwidth = 14,
  width = 22,
  height = 22,
  filename = "cutree_overruled.png",
)

# `minkowski` distance and `ward.D2` clustering
prnHM(
  xmin = -1,
  xmax = 1,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = TRUE,
  cutree_rows = 10,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_row = 3,
  cellwidth = 14,
  width = 18,
  height = 12,
  filter_sp = exprs(species == "human", prot_n_pep >= 2),
  hc_method_rows = "ward.D2", 
  hc_method_cols = "ward.D2", 
  clustering_distance_rows = "minkowski", 
  clustering_distance_cols = "minkowski", 
  p_dist_rows = 2,
  p_dist_cols = 2,
  clustering_distance_cols = "manhattan", 
  filename = "rowminko2_colman_clustward.D2.png",
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
prnHM(
  xmin = -1,
  xmax = 1,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = TRUE,
  cutree_rows = 10,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 3,
  cellwidth = 14,
  filter_sp = exprs(species == "human", prot_n_pep >= 2),
  filter_by = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
  filename = "pval_cutoff_at_1e6.png", 
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

prnHM(
  impute_na = TRUE, 
  xmin = -1,
  xmax = 1,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = TRUE,
  cutree_rows = 10,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_row = 3,
  cellwidth = 14,
  width = 18,
  height = 12,
  filter_prots_by_sp_npep = exprs(species == "human", prot_n_pep >= 3),
  filter_prots_by_pvals = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
  filename = "huprns_fil_impna.png",
)

## peptides
# under selected protein(s)
pepHM(
  xmin = -2,
  xmax = 2,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = TRUE,
  annot_rows = c("gene"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  cellwidth = 12,
  cellheight = 12,
  width = 18,
  height = 12,
  filter_by = exprs(gene %in% c("NCL", "Ncl")),
  filename = "ncl_all.png",
)

# rows ordered by gene 
pepHM(
  xmin = -2,
  xmax = 2,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = FALSE,
  annot_rows = c("gene"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  cellwidth = 12,
  cellheight = 12,
  width = 18,
  height = 12,
  
  filter_by = exprs(gene %in% c("NCL", "Ncl")),
  arrange_peps_by = exprs(gene),
  filename = "ncl_rows_by_gene.png",
)

# rows ordered by sequence 
# (may try alternatively `exprs(pep_seq)` if `pep_seq_mod` not a column key in `Peptide.txt`)
pepHM(
  xmin = -2,
  xmax = 2,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = FALSE,
  annot_rows = c("gene"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  cellwidth = 12,
  cellheight = 12,
  width = 18,
  height = 12,
  filter_by = exprs(gene %in% c("NCL", "Ncl")),
  arrange_peps_by = exprs(pep_seq_mod),
  filename = "ncl_rows_by_seq.png",
)

# more options
pepHM(
  xmin = -2,
  xmax = 2,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = FALSE,
  annot_rows = c("gene", "W16_vs_W2.pVal (W16-W2)"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  cellwidth = 12,
  cellheight = 12,
  width = 18,
  height = 12,
  filter_by = exprs(gene %in% c("NCL", "Ncl")),
  filter_prots_by_pvals = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-5), 
  arrange_by = exprs(gene, -`W16_vs_W2.pVal (W16-W2)`), 
  filename = "ncl_more.png",
)

# selected samples
pepHM(
  col_select = BI_1, 
  xmin = -2,
  xmax = 2,
  xmargin = 0.1,
  annot_cols = c("Group", "Color", "Alpha", "Shape"),
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
  cluster_rows = TRUE,
  annot_rows = c("gene"),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  cellwidth = 12,
  cellheight = 12,
  width = 18,
  height = 12,
  filter_by = exprs(gene %in% c("NCL", "Ncl")),
  arrange_peps_by = exprs(gene),  
  filename = "ncl_bi1.png",
)

## multiple genes
genes <- c("NCL", "Ncl")

lapply(genes, function (gene) {
  gn <- gene
  
  pepHM(
    xmin = -2, 
    xmax = 2, 
    xmargin = 0.1, 
    annot_cols = c("Group", "Color", "Alpha", "Shape"),
    annot_colnames = c("Group", "Lab", "Batch", "WHIM"),
    cluster_rows = FALSE, 
    show_rownames = TRUE, 
    show_colnames = TRUE, 
    fontsize_row = 10, 
    cellwidth = 12, 
    cellheight = 12, 
    width = 18,
    height = 12,
    arrange_pep = exprs(pep_start, pep_end), 
    filter_sp = exprs(gene == !!gn), 
    filename = !!paste0(gene, ".png"),
  )
})

## Customer annotation colors 
annot_colors_group <- colorRampPalette(brewer.pal(n = 9, "Set1"))(12)
names(annot_colors_group) <- c("W16.BI.TMT1", "W16.BI.TMT2", 
                               "W16.JHU.TMT1", "W16.JHU.TMT2", 
                               "W16.PNNL.TMT1", "W16.PNNL.TMT2",
                               "W2.BI.TMT1", "W2.BI.TMT2", 
                               "W2.JHU.TMT1", "W2.JHU.TMT2", 
                               "W2.PNNL.TMT1", "W2.PNNL.TMT2")

annot_colors_lab <- brewer.pal(n = 3, "Set2")
names(annot_colors_lab) <- c("BI", "JHU", "PNNL")

annot_colors_batch <- brewer.pal(n = 4, "Set3")[1:2]
names(annot_colors_batch) <- c("TMT1", "TMT2")

annot_colors_whim <- brewer.pal(n = 4, "Set3")[3:4]
names(annot_colors_whim) <- c("W16", "W2")

annot_colors <- list(Group = annot_colors_group, 
                     Lab = annot_colors_lab, 
                     Batch = annot_colors_batch, 
                     WHIM = annot_colors_whim)

prnHM(
  xmin = -1, 
  xmax = 1, 
  xmargin = 0.1, 
  annot_cols = c("Group", "Color", "Alpha", "Shape"), 
  annot_colnames = c("Group", "Lab", "Batch", "WHIM"), 
  annotation_colors = annot_colors, 
  cluster_rows = TRUE, 
  cutree_rows = 10, 
  show_rownames = FALSE, 
  show_colnames = TRUE, 
  fontsize_row = 3, 
  cellwidth = 14, 
  filter_sp = exprs(species == "human"), 
  filename = custom.png,
)
}

