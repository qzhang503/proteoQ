# ===================================
# Significance test
# ===================================
scale_log2r <- TRUE

## RUN `Mascot or Maxquant but not both`
dontrun <- TRUE
if (!dontrun) {
  ## Mascot
  # prior significance tests in protein abundance changes
	# (a) no NA imputation; output: Protein_pVals.txt
  prnSig(
    impute_na = FALSE,
    W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", "(W2.JHU.TMT2-W2.JHU.TMT1)", "(W2.PNNL.TMT2-W2.PNNL.TMT1)"],
    W2_loc = ~ Term_2["W2.BI-W2.JHU", "W2.BI-W2.PNNL", "W2.JHU-W2.PNNL"],
    W16_vs_W2 = ~ Term_3["W16-W2"],
  )

  # (b) or NA imputation; output: Protein_impNA_pvals.txt
  prnSig(
    impute_na = TRUE,
    W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", "(W2.JHU.TMT2-W2.JHU.TMT1)", "(W2.PNNL.TMT2-W2.PNNL.TMT1)"],
    W2_loc = ~ Term_2["W2.BI-W2.JHU", "W2.BI-W2.PNNL", "W2.JHU-W2.PNNL"],
    W16_vs_W2 = ~ Term_3["W16-W2"],
  )

  # gene-set enrichment tests
	# Protein_pVals.txt will be loaded at `impute_na = FALSE`
  prnGSPA(
    pval_cutoff = 5E-2,
    logFC_cutoff = log2(1.2),
    gspval_cutoff = 5E-2,
    gset_nms = c("go_sets", "kegg_sets"),
    impute_na = FALSE,
  )

  # for proteins with two or more identifying peptides
  prnGSPA(
    pval_cutoff = 5E-2,
    logFC_cutoff = log2(1.2),
    gspval_cutoff = 5E-2,
    gset_nms = c("go_sets", "kegg_sets"),
    filter_by_npep = exprs(prot_n_pep >= 2),
    impute_na = FALSE,
  )

  # customized thresholds for the corresponding formulae in `prnSig()`
  prnGSPA(
    fml_nms = c("W2_bat", "W2_loc", "W16_vs_W2"), # formulae used in `prnSig()`
    pval_cutoff = c(5E-2, 5E-2, 1E-10), # respective protein pVal cut-offs for each formula
    logFC_cutoff = log2(1.2), # i.e., the same protein log2FC cut-off for all formulae
    gspval_cutoff = c(5E-2, 5E-2, 1E-5), # gene-set pVal cut-offs 
    max_size = c(Inf, Inf, 120), # maxixmum sizes of gene sets for consideration
    
    gset_nms = c("go_sets", "kegg_sets"), # the gene sets 
    filter_by_npep = exprs(prot_n_pep >= 2), # only consider proteins with two or more identifying peptides
    impute_na = FALSE, 
  )


  ## MaxQuant
  # prior significance tests in protein abundance changes
	# (a) no NA imputation; output: Protein_pVals.txt	
  prnSig(
    impute_na = FALSE, 
    W16_vs_W2_fine = ~ Term["W16.BI-W2.BI", "W16.JHU-W2.JHU", "W16.PNNL-W2.PNNL"],
    W16_vs_W2_course = ~ Term_2["W16-W2"], 
  )
	
  # (b) or NA imputation; output: Protein_impNA_pvals.txt
  prnSig(
    impute_na = TRUE, 
    W16_vs_W2_fine = ~ Term["W16.BI-W2.BI", "W16.JHU-W2.JHU", "W16.PNNL-W2.PNNL"],
    W16_vs_W2_course = ~ Term_2["W16-W2"], 
  )
	
  # gene-set enrichment tests
	# Protein_pVals.txt will be loaded at `impute_na = FALSE`	
  prnGSPA(
    pval_cutoff = 5E-2,
    logFC_cutoff = log2(1.2),
    gspval_cutoff = 5E-2,
    gset_nms = c("go_sets", "kegg_sets"),
    impute_na = FALSE,
  )
	
  # customized thresholds for the corresponding formulae in `prnSig()`
  prnGSPA(
    fml_nms = c("W16_vs_W2_fine", "W16_vs_W2_course"), # formulae used in `prnSig()`
    pval_cutoff = c(5E-2, 5E-4), # respective protein pVal cut-offs for each formula
    logFC_cutoff = log2(1.2), # i.e., the same protein log2FC cut-off for all formulae
    gspval_cutoff = c(5E-2, 1E-6), # gene-set pVal cut-offs 
    max_size = c(Inf, 120), # maxixmum sizes of gene sets for consideration
    
    gset_nms = c("go_sets", "kegg_sets"), # the gene sets 
    filter_by_npep = exprs(prot_n_pep >= 2), # only consider proteins with two or more identifying peptides
    impute_na = FALSE,
  )
}
## END of RUN `Mascot or Maxquant but not both`


