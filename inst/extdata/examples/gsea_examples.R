# ===================================
# Significance test
# ===================================
scale_log2r <- TRUE

## RUN `Mascot or Maxquant but not both`
dontrun <- TRUE
if (!dontrun) {
  ## Mascot
	
  # suggested prior significance tests of peptide log2FC
  # (1a) no NA imputation; output: Peptide_pVals.txt
  pepSig(
    impute_na = FALSE, 
    cmbn_contrs = ~ Term["(W16.BI.TMT1+W16.BI.TMT2)/2-(W2.JHU.TMT1+W2.JHU.TMT2)/2"],
    W16_vs_W2_fine = ~ Term["W16.BI.TMT1-W2.BI.TMT1", "W16.JHU.TMT1-W2.JHU.TMT1", "W16.PNNL.TMT1-W2.PNNL.TMT1"],
  )
	
  # (1b) or NA imputation; output: Peptide_impNA_pvals.txt
  pepSig(
    impute_na = TRUE, 
    cmbn_contrs = ~ Term["(W16.BI.TMT1+W16.BI.TMT2)/2-(W2.JHU.TMT1+W2.JHU.TMT2)/2"],
    W16_vs_W2_fine = ~ Term["W16.BI.TMT1-W2.BI.TMT1", "W16.JHU.TMT1-W2.JHU.TMT1", "W16.PNNL.TMT1-W2.PNNL.TMT1"],
  )
	
  # required prior significance tests of protein log2FC
  # all formulae in the latest `pepSig()` will be assessed by default 
  # (2a) no NA imputation; output: Protein_pVals.txt
  prnSig(impute_na = FALSE)

  # (2b) or NA imputation; output: Protein_impNA_pvals.txt
  prnSig(impute_na = TRUE)

	
  ## MaxQuant
  # suggested prior significance tests of peptide log2FC
  # (1a) no NA imputation; output: Peptide_pVals.txt
  pepSig(
    impute_na = FALSE, 
    cmbn_contrs = ~ Term["(W16.BI+W16.JHU+W16.PNNL)/3-(W2.BI+W2.JHU+W2.PNNL)/3"],
    W16_vs_W2_fine = ~ Term["W16.BI-W2.BI", "W16.JHU-W2.JHU", "W16.PNNL-W2.PNNL"],
  )
	
  # (1b) or NA imputation; output: Peptide_impNA_pvals.txt
  pepSig(
    impute_na = TRUE, 
    cmbn_contrs = ~ Term["(W16.BI+W16.JHU+W16.PNNL)/3-(W2.BI+W2.JHU+W2.PNNL)/3"],
    W16_vs_W2_fine = ~ Term["W16.BI-W2.BI", "W16.JHU-W2.JHU", "W16.PNNL-W2.PNNL"],
  )
  
  # see the above Mascot example for remaining steps
  
}
## END of RUN `Mascot or Maxquant but not both`


## examples
# (1a) all human proteins; input: Protein_pVals.txt
prnGSEA(
  impute_na = FALSE, 
  var_cutoff = 0, 
  pval_cutoff = 1, 
  logFC_cutoff = log2(1), 
  filter_by_sp = exprs(species == "human"), 
)

# (1b) all human proteins; input: Protein_impNA_pvals.txt
prnGSEA(
  impute_na = TRUE, 
  var_cutoff = 0, 
  pval_cutoff = 1, 
  logFC_cutoff = log2(1), 
  filter_by_sp = exprs(species == "human"), 
)

# (2) prefiltration by variances, pVals and logFCs
prnGSEA(
  impute_na = FALSE, 
  var_cutoff = 0.5,     
  pval_cutoff = 5E-2,
  logFC_cutoff = log2(1.2),
  filter_by_sp = exprs(species == "human"), 
)

# (3) cases that are complete with no missing values
prnGSEA(
  impute_na = FALSE, 
  var_cutoff = 0.5,     
  pval_cutoff = 5E-2,
  logFC_cutoff = log2(1.2),
  complete_cases = TRUE, 
  filter_by_sp = exprs(species == "human"), 
)


