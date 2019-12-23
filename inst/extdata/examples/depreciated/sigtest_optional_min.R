# ===================================
# Optional significance tests
# (for data filtration by pVals...)
# ===================================
scale_log2r <- TRUE

## RUN `Mascot or Maxquant but not both`
dontrun <- TRUE
if (!dontrun) {
  ## Mascot
  # peptide significance tests
	# (1a) no NA imputation; output: Peptide_pVals.txt
  pepSig(
    impute_na = FALSE, 
    W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", "(W2.JHU.TMT2-W2.JHU.TMT1)", "(W2.PNNL.TMT2-W2.PNNL.TMT1)"], # batch effects
    W2_loc = ~ Term_2["W2.BI-W2.JHU", "W2.BI-W2.PNNL", "W2.JHU-W2.PNNL"], # location effects
    W16_vs_W2 = ~ Term_3["W16-W2"], 
  )
  
  # (1b) or NA imputation; output: Peptide_impNA_pvals.txt
  pepSig(
    impute_na = TRUE, 
    W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", "(W2.JHU.TMT2-W2.JHU.TMT1)", "(W2.PNNL.TMT2-W2.PNNL.TMT1)"], # batch effects
    W2_loc = ~ Term_2["W2.BI-W2.JHU", "W2.BI-W2.PNNL", "W2.JHU-W2.PNNL"], # location effects
    W16_vs_W2 = ~ Term_3["W16-W2"], 
  )
	
  # protein significance tests
  # all formulae in the latest `pepSig()` will be assessed by default 
  # (2a) no NA imputation; output: Protein_pVals.txt
  prnSig(impute_na = FALSE)
  
  # (2b) or NA imputation; output: Protein_impNA_pvals.txt
  prnSig(impute_na = TRUE)
  
	
  ## MaxQuant
  # peptide significance tests
	# (1a) no NA imputation; output: Peptide_pVals.txt
  pepSig(
    impute_na = FALSE, 
    W16_vs_W2_fine = ~ Term["W16.BI-W2.BI", "W16.JHU-W2.JHU", "W16.PNNL-W2.PNNL"],
    W16_vs_W2_course = ~ Term_2["W16-W2"], 
  )
	
  # (1b) or NA imputation; output: Peptide_impNA_pvals.txt
  pepSig(
    impute_na = TRUE, 
    W16_vs_W2_fine = ~ Term["W16.BI-W2.BI", "W16.JHU-W2.JHU", "W16.PNNL-W2.PNNL"],
    W16_vs_W2_course = ~ Term_2["W16-W2"], 
  )
  
  # protein significance tests
	# formula(e) in the latest `pepSig()` will be taken by default
  # (2a) no NA imputation; output: Protein_pVals.txt
  prnSig(impute_na = FALSE)
  
  # (2b) or NA imputation; output: Protein_impNA_pvals.txt
  prnSig(impute_na = TRUE)
}
## END of RUN `Mascot or Maxquant but not both`


