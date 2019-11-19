# ===================================
# Significance tests
# ===================================
scale_log2r <- TRUE

## RUN `Mascot or Maxquant but not both`
dontrun <- TRUE
if (!dontrun) {
  # Mascot
	# peptide significance tests
	pepSig(
		impute_na = FALSE, 
		W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", "(W2.JHU.TMT2-W2.JHU.TMT1)", "(W2.PNNL.TMT2-W2.PNNL.TMT1)"], # batch effects
		W2_loc = ~ Term_2["W2.BI-W2.JHU", "W2.BI-W2.PNNL", "W2.JHU-W2.PNNL"], # location effects
		W16_vs_W2 = ~ Term_3["W16-W2"], 
	)

	# protein significance tests
	prnSig(
		impute_na = FALSE, 
		W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", "(W2.JHU.TMT2-W2.JHU.TMT1)", "(W2.PNNL.TMT2-W2.PNNL.TMT1)"], # batch effects
		W2_loc = ~ Term_2["W2.BI-W2.JHU", "W2.BI-W2.PNNL", "W2.JHU-W2.PNNL"], # location effects
		W16_vs_W2 = ~ Term_3["W16-W2"], # between two WHIMs
	)

	# more examples
	prnSig(
		impute_na = FALSE,
		W2_loc_2 = ~ Term["(W2.BI.TMT2+W2.BI.TMT1)/2 - (W2.JHU.TMT2+W2.JHU.TMT1)/2"], # locations with batch average
	)	

	prnSig(
		impute_na = TRUE,
		W2_vs_W16_fix = ~ Term_3["W16-W2"], # fixed effect
		W2_vs_W16_mix = ~ Term_3["W16-W2"] + (1|TMT_Set), # one fixed and one random effect
	)

	prnSig(
		impute_na = TRUE,
		method = lm,
		W2_vs_W16_fix = ~ Term_3["W16-W2"], # one fixed effect
		W2_vs_W16_mix = ~ Term_3["W16-W2"] + (1|TMT_Set), # one fixed and one random effect
		W2_vs_W16_mix_2 = ~ Term_3["W16-W2"] + (1|TMT_Set) + (1|Color), # one fixed and two random effects
	)
	
  # MaxQuant
	# peptide significance tests
	pepSig(
		impute_na = FALSE, 
		W16_vs_W2_fine = ~ Term["W16.BI-W2.BI", "W16.JHU-W2.JHU", "W16.PNNL-W2.PNNL"],
		W16_vs_W2_course = ~ Term_2["W16-W2"], 
	)

	# protein significance tests
	prnSig(
		impute_na = FALSE, 
		W16_vs_W2_fine = ~ Term["W16.BI-W2.BI", "W16.JHU-W2.JHU", "W16.PNNL-W2.PNNL"],
		W16_vs_W2_course = ~ Term_2["W16-W2"], 
	)
}
## END of RUN `Mascot or Maxquant but not both`







