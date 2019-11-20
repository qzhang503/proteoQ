# ===================================
# Prerequisite PSM normalization
# ===================================
library(proteoQ)

## RUN `Mascot or Maxquant but not both`
dontrun <- TRUE
if (!dontrun) {
  # Mascot
  dat_dir <- "C:\\The\\Mascot\\Example"

  # Maxquant
  dat_dir <- c("C:\\The\\MQ\\Example")
}
## END of RUN `Mascot or Maxquant but not both`

dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)

load_expts()

normPSM(
  group_psm_by = pep_seq_mod, 
  fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta", 
            "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"), 
)


