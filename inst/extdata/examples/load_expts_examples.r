# An examplary "expt_smry.xlsx"
system.file("extdata", "expt_smry.xlsx", package = "proteoQ")

# An examplary "frac_smry.xlsx"
system.file("extdata", "frac_smry.xlsx", package = "proteoQ")


# ===================================
# Fasta and PSM files
# ===================================
# fasta
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
  
  # PSM
  cptac_csv_1(dat_dir)
  
  # "expt_smry.xlsx" and "frac_smry.xlsx"
  cptac_expt_1(dat_dir)
  cptac_frac_1(dat_dir)
  
  # or Maxquant
  #   to install package `proteoQDB` containing MaxQuant examples: 
  #    (1) clone https://github.com/qiangzhang503/proteoQDB.git through `Github Desktop`
  #    (2) local installation of proteoQDB, e.g., devtools::install("~\\my_dir\\proteoQDB")
  library(proteoQDB)
  dat_dir <- c("C:\\The\\MQ\\Example")
  dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)
	
  # PSM
  cptac_mqpsm_txt(dat_dir)
  
  # "expt_smry.xlsx" and "frac_smry.xlsx"
  cptac_mqpsm_expt(dat_dir)
  cptac_mqpsm_frac(dat_dir)	
}
## END of RUN `Mascot or Maxquant but not both`


# ===================================
# PSM, peptide and protein processing
# ===================================
# load experiments
library(proteoQ)
load_expts()

# PSM processing with in-function filtration of data by `filter_`
normPSM(
  group_psm_by = pep_seq_mod,
  group_pep_by = gene,
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

# exemplary PSM purging
purgePSM(pt_cv = .95)

# peptides results with exemplary `filter_...`
normPep(
  method_psm_pep = median,
  method_align = MGKernel,
  range_log2r = c(5, 95),
  range_int = c(5, 95),
  n_comp = 3,
  seed = 749662,
  maxit = 200,
  epsilon = 1e-05,
  # filter_by_sp = exprs(species == "human"),
)

# exemplary peptide purging
purgePep(pt_cv = .95)

# proteins results with examplary `filter_...`
normPrn(
  method_pep_prn = median,
  method_align = MGKernel,
  range_log2r = c(20, 95),
  range_int = c(5, 95),
  n_comp = 2,
  seed = 749662,
  maxit = 200,
  epsilon = 1e-05,
  filter_prns = exprs(prot_n_pep >= 2)
)


