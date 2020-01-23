# ===================================
# GSVA
# ===================================

## !!!require the brief working example in `?load_expts`

## global option
scale_log2r <- TRUE

## base
prnGSVA(
  impute_na = FALSE,
  min.sz = 10,
  verbose = FALSE,
  parallel.sz = 0,
  mx.diff = TRUE,
  gset_nms = "go_sets",
)

## row filtration
prnGSVA(
  impute_na = FALSE,
  min.sz = 10,
  verbose = FALSE,
  parallel.sz = 0,
  mx.diff = TRUE,
  gset_nms = c("go_sets", "kegg_sets"),
  filter_prots = exprs(prot_n_pep >= 3), 
  filename = fil.txt,
)

## additional row filtration by pVals (impute_na = FALSE)
# if not yet, run prerequisitive significance tests at `impute_na = FALSE`
pepSig(
  impute_na = FALSE, 
  W2_bat = ~ Term["W2.BI.TMT2-W2.BI.TMT1", 
                  "W2.JHU.TMT2-W2.JHU.TMT1", 
                  "W2.PNNL.TMT2-W2.PNNL.TMT1"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = FALSE)

# (`W16_vs_W2.pVal (W16-W2)` is now a column key)
prnGSVA(
  min.sz = 10,
  verbose = FALSE,
  parallel.sz = 0,
  mx.diff = TRUE,
  gset_nms = "go_sets",
  filter_prots_by_npep = exprs(prot_n_pep >= 3), 
  filter_prots_by_pval = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
)

## additional row filtration by pVals (impute_na = TRUE)
# if not yet, run prerequisitive NA imputation
pepImp(m = 2, maxit = 2)
prnImp(m = 5, maxit = 5)

# if not yet, run prerequisitive significance tests at `impute_na = TRUE`
pepSig(
  impute_na = TRUE, 
  W2_bat = ~ Term["W2.BI.TMT2-W2.BI.TMT1", 
                  "W2.JHU.TMT2-W2.JHU.TMT1", 
                  "W2.PNNL.TMT2-W2.PNNL.TMT1"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = TRUE)

prnGSVA(
  impute_na = TRUE, 
  min.sz = 10,
  verbose = FALSE,
  parallel.sz = 0,
  mx.diff = TRUE,
  gset_nms = "go_sets",
  filter_prots_by_npep = exprs(prot_n_pep >= 3), 
  filter_prots_by_pval = exprs(`W16_vs_W2.pVal (W16-W2)` <= 1e-6), 
)
