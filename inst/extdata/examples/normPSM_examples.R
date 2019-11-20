# ===================================
# PSM normalization
# ===================================
library(proteoQ)

## RUN `Mascot or Maxquant but not both`
dontrun <- TRUE
if (!dontrun) {
  # Mascot
  dat_dir <- "C:\\The\\Mascot\\Example"
  dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)
  load_expts()
	
  normPSM(
    group_psm_by = pep_seq_mod,
    group_pep_by = prot_acc,
    fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
              "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"),
    rptr_intco = 3000,
    rm_craps = TRUE,
    rm_krts = FALSE,
    rm_outliers = FALSE,
    annot_kinases = TRUE,
    plot_rptr_int = TRUE,
    plot_log2FC_cv = TRUE,
    
    filter_peps = exprs(pep_expect <= .1),
    filter_by_more = exprs(pep_rank == 1, pep_exp_z > 1),
  )

  # MaxQuant
  dat_dir <- c("C:\\The\\MQ\\Example")
  dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)
  load_expts()
	
  normPSM(
    group_psm_by = pep_seq_mod,
    group_pep_by = prot_acc,
    fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
              "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"),
    corrected_int = TRUE,
    rm_reverses = TRUE,
    rptr_intco = 3000,
    rm_craps = TRUE,
    rm_krts = FALSE,
    rm_outliers = FALSE,
    annot_kinases = TRUE,
    plot_rptr_int = TRUE,
    plot_log2FC_cv = TRUE,
    
    filter_peps = exprs(PEP <= 0.1),
  )
}
## END of RUN `Mascot or Maxquant but not both`



