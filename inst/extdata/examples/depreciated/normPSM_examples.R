# ===================================
# PSM normalization
# ===================================
library(proteoQ)
dat_dir <- "~\\proteoQ\\examples"
load_expts()

# PSM (choose one of the platforms)
choose_one <- TRUE
if (!choose_one) {
  # Mascot
  normPSM(
    group_psm_by = pep_seq_mod,
    group_pep_by = prot_acc,
    fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
              "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"),

    # variable argument statement(s)
    filter_psms_at = exprs(pep_expect <= .1),
    filter_by_more = exprs(pep_rank == 1, pep_exp_z > 1),
  )
  
  # or MaxQuant
  normPSM(
    group_psm_by = pep_seq_mod,
    group_pep_by = prot_acc,
    fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
              "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"),
    corrected_int = TRUE,
    rm_reverses = TRUE,

    # vararg statement(s)
    filter_psms_at = exprs(PEP <= 0.1),
  )
  
  # or Spectrum Mill
  normPSM(
    group_psm_by = pep_seq_mod,
    group_pep_by = prot_acc,
    fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
              "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"),

    # vararg statement(s)
    filter_psms = exprs(score >= 10),
  )
}

