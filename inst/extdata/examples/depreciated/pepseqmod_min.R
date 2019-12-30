# ===================================
# Prerequisite PSM normalization
# ===================================
library(proteoQ)
dat_dir <- "~\\proteoQ\\examples"
load_expts()

normPSM(
  group_psm_by = pep_seq_mod, 
  fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta", 
            "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"), 
)

