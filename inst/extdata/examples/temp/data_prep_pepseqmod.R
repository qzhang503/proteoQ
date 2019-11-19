# ===================================
# Fasta and PSM files
# ===================================
library(proteoQDA)
fasta_dir <- "~\\proteoQ\\dbs\\fasta\\refseq"
dir.create(fasta_dir, recursive = TRUE, showWarnings = FALSE)
copy_refseq_hs(fasta_dir)
copy_refseq_mm(fasta_dir)

## RUN `Mascot or Maxquant but not both`
dontrun <- TRUE
if (!dontrun) {
  # Mascot
  dat_dir <- "C:\\The\\Mascot\\Example"
  dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)
  cptac_csv_1(dat_dir)
  cptac_expt_1(dat_dir)
  cptac_frac_1(dat_dir)
  
  # Maxquant
  #   to install package `proteoQDB` containing MaxQuant examples: 
  #    (1) clone https://github.com/qiangzhang503/proteoQDB.git through `Github Desktop`
  #    (2) local installation of proteoQDB, e.g., devtools::install("~\\my_dir\\proteoQDB")
  library(proteoQDB)
  dat_dir <- c("C:\\The\\MQ\\Example")
  dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)
  cptac_mqpsm_txt(dat_dir)
  cptac_mqpsm_expt(dat_dir)
  cptac_mqpsm_frac(dat_dir)	
}
## END of RUN `Mascot or Maxquant but not both`


# ===================================
# PSM, peptide and protein processing
# ===================================
library(proteoQ)
load_expts()

normPSM(
	group_psm_by = pep_seq_mod, 
	fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta", 
            "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"), 
)

normPep(
  method_align = MGKernel, 
  n_comp = 3, 
  seed = 749662, 
)

normPrn(
  method_align = MGKernel, 
  n_comp = 2, 
  seed = 749662, 
)


