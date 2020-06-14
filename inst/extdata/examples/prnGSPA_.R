\donttest{
# ===================================
# GSPA tests
# ===================================

## !!!require the brief working example in `?load_expts`

## global option
scale_log2r <- TRUE

## base
# prerequisites in significance tests (impute_na = FALSE)
# as `pVals` are required
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

# GSPA analysis
prnGSPA(
  impute_na = FALSE,
  pval_cutoff = 5E-2, # protein pVal threshold
  logFC_cutoff = log2(1.2), # protein log2FC threshold
  gspval_cutoff = 5E-2, # gene-set pVal threshold
  gslogFC_cutoff = log2(1.2), # gene-set log2FC threshold
  gset_nms = c("go_sets", "c2_msig"), 
)

# GSPA visualization under volcano plots 
# forced lines of `pval_cutoff` and `logFC_cutoff`  
#   from the corresponding `prnGSPA`;
# more examples of aesthetics in `?gspaMap`
gspaMap(
  impute_na = FALSE,
  show_labels = TRUE, 
  topn = 50, 
  show_sig = pVal, 
  xco = 1.2, # optional x-axis cut-off lines
  yco = 0.05, # a optional y-axis cut-off line 
)

## row filtration
prnGSPA(
  impute_na = FALSE,
  gset_nms = c("go_sets", "c2_msig"),
  filter_prots_by_npep = exprs(prot_n_pep >= 2),
  filter_prots_by_species = exprs(species == "human"), 
  filename = row_filters.txt,
)

# by default, volcano-plot visualization for all GSPA result files, 
#   which are "Protein_GSPA_Z.txt" and "row_filters_Protein_GSPA_Z.txt" 
#   up to this point
# will match the primary `fliter_` varargs to 
#   those in the corresponding `prnGSPA(...)`
gspaMap(
  impute_na = FALSE,
  show_labels = FALSE, 
  topn = 50, 
  show_sig = pVal, 
)

# or manually supply specific secondary file(s) with `df2`; 
gspaMap(
  df2 = "row_filters_Protein_GSPA_Z.txt", 
  impute_na = FALSE,
  topn = 10, 
  show_labels = TRUE, 
  show_sig = pVal, 
)

# may apply additionally secondary `fliter2_` varargs against `df2`
gspaMap(
  df2 = "row_filters_Protein_GSPA_Z.txt", 
  filter2_ess_size = exprs(ess_size >= 1),   
  impute_na = FALSE,
  topn = 10, 
  show_labels = FALSE, 
  show_sig = pVal, 
)

# can also visualize results for specific formula(e) and/or 
#   specific `df2`
gspaMap(
  fml_nms = "W16_vs_W2",
  df2 = "row_filters_Protein_GSPA_Z.txt", 
  filter2_ess_size = exprs(ess_size >= 1),   
  impute_na = FALSE,
  topn = 20, 
  show_labels = TRUE, 
  show_sig = pVal, 
)

## customized thresholds
# (formula names defined a priori in `pepSig(...)`)
# (higher thresholds for "W16_vs_W2"; also set `max_size` to 
#  avoid trapping with large terms, such as "cell parts"...)
prnGSPA(
  impute_na = FALSE, 
  fml_nms = c("W2_bat", "W2_loc", "W16_vs_W2"),
  pval_cutoff = c(5E-2, 5E-2, 1E-6), 
  logFC_cutoff = c(log2(1.2)), 
  gslogFC_cutoff = c(log2(1.2)), 
  gspval_cutoff = c(5E-2, 5E-2, 1E-4), 
  max_size = c(Inf, Inf, 120), 
  gset_nms = c("go_sets"), 
  filter_by_npep = exprs(prot_n_pep >= 2), 
  filter_prots_by_species = exprs(species == "human"), 
  filename = thresholds.txt,
)

gspaMap(
  df2 = "thresholds_Protein_GSPA_Z.txt", 
  impute_na = FALSE,
  topn = 10, 
  show_labels = FALSE, 
  show_sig = pVal, 
)

# "limma" for "W16_vs_W2"
prnGSPA(
  impute_na = FALSE, 
  fml_nms = c("W2_bat", "W2_loc", "W16_vs_W2"),
  method = c("mean", "mean", "limma"),
  max_size = c(Inf, Inf, 120), 
  pval_cutoff = c(5E-2), 
  logFC_cutoff = c(log2(1.2)), 
  gspval_cutoff = c(5E-2), 
  gset_nms = c("go_sets", "c2_msig"), 
  filter_by = exprs(prot_n_pep >= 2, species == "human"), 
  filename = methods.txt,
)

gspaMap(
  df2 = "methods_Protein_GSPA_Z.txt", 
  impute_na = FALSE,
  topn = 25, 
  show_labels = FALSE, 
  show_sig = pVal, 
)

## NA imputation
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

prnGSPA(
  impute_na = TRUE,
  gset_nms = c("go_sets"),
  filter_prots_by_npep = exprs(prot_n_pep >= 2),
)

# ===================================
# Distance heat maps of gene sets
# (also interactive networks)
# ===================================
# a `term` is a subset of an `ess_term` if the distance is zero
prnGSPAHM(
  filter2_by = exprs(distance <= .7),
  annot_cols = "ess_idx",
  annot_colnames = "Eset index",
  annot_rows = "ess_size",
  filename = show_some_redundancy.png,
)

# human terms only
prnGSPAHM(
  filter2_num = exprs(distance <= .95),
  filter2_sp = exprs(start_with_str("hs", term)),
  annot_cols = "ess_idx",
  annot_colnames = "Eset index",
  filename = show_more_connectivity.png,
)

# custom color palette
prnGSPAHM(
  annot_cols = c("ess_idx", "ess_size"),
  annot_colnames = c("Eset index", "Size"),
  filter2_by = exprs(distance <= .95),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  filename = custom_colors.png,
)

# can also visualize results for specific `df2` and formula(e)
prnGSPAHM(
  df2 = "thresholds_Protein_GSPA_Z_essmap.txt", 
  fml_nms = "W16_vs_W2", 
  annot_cols = c("ess_idx", "ess_size"),
  annot_colnames = c("Eset index", "Size"),
  filter2_by = exprs(distance <= .95),
  filename = w16_vs_w2.png,
)


# ===================================
# custom GO databases
# ===================================
# start over
unlink(file.path(dat_dir, "Protein/GSPA"), recursive = TRUE, force = TRUE)

# prepare GO databases (see also ?prepGO)
prepGO(
  species = human,
  db_path = "~/proteoQ/dbs/go",
  gaf_url = "http://current.geneontology.org/annotations/goa_human.gaf.gz",
  obo_url = "http://purl.obolibrary.org/obo/go/go-basic.obo",
  filename = go_hs.rds,
)

prepGO(
  species = mouse,
  db_path = "~/proteoQ/dbs/go",
  gaf_url = "http://current.geneontology.org/annotations/mgi.gaf.gz",
  obo_url = "http://purl.obolibrary.org/obo/go/go-basic.obo",
  filename = go_mm.rds,
)

head(readRDS(file.path("~/proteoQ/dbs/go", "go_hs.rds")))
head(readRDS(file.path("~/proteoQ/dbs/go", "go_mm.rds")))

# analysis
prnGSPA(
  impute_na = FALSE,
  gset_nms = c("~/proteoQ/dbs/go/go_hs.rds",
               "~/proteoQ/dbs/go/go_mm.rds"),
)

# It is possible that the system "go_sets" may be partially different to the
# custom GO. Be explicit on `gset_nms`; otherwise, the default of `gset_nms =
# c("go_sets", "c2_msig")` will be applied.
gspaMap(
  gset_nms = c("~/proteoQ/dbs/go/go_hs.rds",
               "~/proteoQ/dbs/go/go_mm.rds"),
  impute_na = FALSE,
  topn = 20, 
  show_sig = pVal, 
)

# ===================================
# custom MSig databases
# ===================================
# start over
unlink(file.path(dat_dir, "Protein/GSPA"), recursive = TRUE, force = TRUE)

# prepare MSig databases (see also ?prepMSig)
prepMSig(species = human)
prepMSig(species = mouse)

head(readRDS(file.path("~/proteoQ/dbs/msig", "msig_hs.rds")))
head(readRDS(file.path("~/proteoQ/dbs/msig", "msig_mm.rds")))

# analysis
prnGSPA(
  impute_na = FALSE,
  gset_nms = c("~/proteoQ/dbs/msig/msig_hs.rds",
               "~/proteoQ/dbs/msig/msig_mm.rds"),
)

# It is possible that the system "c2_msig" may be partially different to the
# custom databases. Be explicit on `gset_nms`; otherwise, the default of
# `gset_nms = c("go_sets", "c2_msig")` will be applied.
gspaMap(
  gset_nms = c("~/proteoQ/dbs/msig/msig_hs.rds",
               "~/proteoQ/dbs/msig/msig_mm.rds"),
  impute_na = FALSE,
  topn = 20, 
  show_sig = pVal, 
)

## put custom `prepGO`, `prepMSig` preparation into `prnGSPA` analysis
# start over with fresh database preparation
unlink(file.path("~/proteoQ/dbs/go"), recursive = TRUE, force = TRUE)
unlink(file.path("~/proteoQ/dbs/msig"), recursive = TRUE, force = TRUE)

library(magrittr)

prnGSPA(
  impute_na = FALSE,
  pval_cutoff = 5E-2,
  logFC_cutoff = log2(1.2),
  gspval_cutoff = 5E-2,
  gset_nms = c(prepGO(human), prepGO(mouse), prepMSig(human), prepMSig(mouse)),
  filename = cmbn.txt,
) %>% 
  gspaMap(
    # gset_nms = ., 
    df2 = "cmbn_Protein_GSPA_Z.txt", 
    impute_na = FALSE,
    topn = 20, 
    show_labels = FALSE, 
    show_sig = pVal, 
  )
}


