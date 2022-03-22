\donttest{
# ***********************************
# ************    TMT    ************
# ***********************************
  
# ===================================
# Fasta and PSM files
# ===================================
# fasta (all platforms)
library(proteoQDA)
fasta_dir <- "~/proteoQ/dbs/fasta/refseq"
copy_refseq_hs(fasta_dir)
copy_refseq_mm(fasta_dir)

# working directory (all platforms)
dat_dir <- "~/proteoQ/examples"

# metadata (all platforms)
copy_exptsmry_gtmt(dat_dir)
copy_fracsmry_gtmt(dat_dir)

# PSM (choose one of the platforms)
choose_one <- TRUE
if (!choose_one) {
  ## Mascot
  copy_mascot_gtmt(dat_dir)
  
  ## or MaxQuant
  # copy_maxquant_gtmt(dat_dir)
  
  ## or MSFragger
  # copy_msfragger_gtmt(dat_dir)
  
  ## or proteoM
  # copy_proteom_gtmt(dat_dir)
  
  ## or Spectrum Mill
  # (temporarily unavailable)
}

# ===================================
# PSM, peptide and protein processing
# ===================================
library(proteoQ)
load_expts("~/proteoQ/examples")

# PSM data standardization
normPSM(
  group_psm_by = pep_seq_mod, 
  group_pep_by = gene, 
  annot_kinases = TRUE, 
  
  # no default and required
  fasta = c("~/proteoQ/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
            "~/proteoQ/dbs/fasta/refseq/refseq_mm_2013_07.fasta"),
)

# optional PSM purging
purgePSM()

# PSMs to peptides
PSM2Pep()

# peptide data merging
mergePep()

# peptide data standardization
standPep()

# peptide data histograms
pepHist()

# optional peptide purging
purgePep()

# peptides to proteins
Pep2Prn(use_unique_pep = TRUE)

# protein data standardization
standPrn()

# protein data histograms
prnHist()

# ===================================
# Optional significance tests
# (no NA imputation)
# ===================================
pepSig(
  W2_bat = ~ Term["W2.BI.TMT2-W2.BI.TMT1", 
                  "W2.JHU.TMT2-W2.JHU.TMT1", 
                  "W2.PNNL.TMT2-W2.PNNL.TMT1"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig()

# ===================================
# optional NA imputation
# ===================================
pepImp(m = 2, maxit = 2)
prnImp(m = 5, maxit = 5)

# ===================================
# Optional significance tests
# (with NA imputation)
# ===================================
pepSig(
  impute_na = TRUE, 
  W2_bat = ~ Term["W2.BI.TMT2-W2.BI.TMT1", 
                  "W2.JHU.TMT2-W2.JHU.TMT1", 
                  "W2.PNNL.TMT2-W2.PNNL.TMT1"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = TRUE)



# ***********************************
# ************    LFQ    ************
# ***********************************

# ===================================
# Fasta and PSM files
# ===================================
# fasta (all platforms)
library(proteoQDA)
fasta_dir <- "~/proteoQ/dbs/fasta/uniprot"
copy_uniprot_hsmm(fasta_dir)

# working directory (all platforms)
dat_dir <- "~/proteoQ/examples"

# metadata (all platforms)
copy_exptsmry_plfq(dat_dir)
copy_fracsmry_plfq(dat_dir)

# PSM (choose one of the platforms)
choose_one <- TRUE
if (!choose_one) {
  ## Mascot
  copy_mascot_plfq(dat_dir)
  
  ## or MaxQuant
  # copy_maxquant_plfq(dat_dir)
  
  ## or MSFragger
  # copy_msfragger_plfq(dat_dir)
  
  ## or proteoM
  # copy_proteom_plfq(dat_dir)
  
  ## or Spectrum Mill
  # (temporarily unavailable)
}


# ===================================
# PSM, peptide and protein processing
# ===================================
library(proteoQ)
load_expts("~/proteoQ/examples")

# PSM data standardization
normPSM(
  group_psm_by = pep_seq_mod, 
  group_pep_by = gene, 
  annot_kinases = TRUE, 
  fasta = c("~/proteoQ/dbs/fasta/uniprot/uniprot_hsmm_2020_03.fasta"),
)

# PSM purging not applicable with LFQ
# purgePSM()

# PSMs to peptides
PSM2Pep()

# peptide data merging
mergePep()

# peptide data standardization
standPep()

# peptide data histograms
pepHist()

# optional peptide purging
purgePep()

# peptides to proteins
Pep2Prn(use_unique_pep = TRUE)

# protein data standardization
standPrn()

# protein data histograms
prnHist()

# ===================================
# Optional significance tests
# (no NA imputation)
# ===================================
pepSig(
  fml_1 = ~ Term["BI-JHU", 
                 "JHU-PNNL", 
                 "(BI+JHU)/2-PNNL"],
)

prnSig()

# ===================================
# optional NA imputation
# ===================================
pepImp(m = 2, maxit = 2)
prnImp(m = 5, maxit = 5)

# ===================================
# Optional significance tests
# (with NA imputation)
# ===================================
pepSig(
  impute_na = TRUE, 
  fml_1 = ~ Term["BI-JHU", 
                 "JHU-PNNL", 
                 "(BI+JHU)/2-PNNL"],
)

prnSig(impute_na = TRUE)


# ***********************************
# ***********    SILAC    ***********
# ***********************************

# Database searches
library(proteoM)

matchMS(
  silac_mix = list(base = NULL, heavy = c("K8 (K)", "R10 (R)")),
  ...
)

# The remaining is the same as LFQ
# ...

}


\dontrun{
load_expts(dat_dir = "~/proteoQ/examples", expt_smry = "expt_smry.xlsx")

# not working; `expt_smry = my_expt` is an expression
my_expt <- "expt_smry.xlsx"
load_expts(dat_dir = "~/proteoQ/examples", expt_smry = my_expt)

# need unquoting; 
# see also: https://dplyr.tidyverse.org/articles/programming.html
load_expts(dat_dir = "~/proteoQ/examples", expt_smry = !!my_expt)
}



