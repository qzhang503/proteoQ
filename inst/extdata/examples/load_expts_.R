\donttest{
# ===================================
# Fasta and PSM files
# ===================================
# fasta (all platforms)
library(proteoQDA)
fasta_dir <- "~\\proteoQ\\dbs\\fasta\\refseq"
dir.create(fasta_dir, recursive = TRUE, showWarnings = FALSE)
copy_refseq_hs(fasta_dir)
copy_refseq_mm(fasta_dir)

# working directory (all platforms)
dat_dir <- "~\\proteoQ\\examples"
dir.create(dat_dir, recursive = TRUE, showWarnings = FALSE)

# metadata (all platforms)
copy_global_exptsmry(dat_dir)
copy_global_fracsmry(dat_dir)

# PSM (choose one of the platforms)
choose_one <- TRUE
if (!choose_one) {
  # Mascot
  copy_global_mascot(dat_dir)
  
  # or MaxQuant
  copy_global_maxquant(dat_dir)
  
  # or Spectrum Mill
  copy_global_sm(dat_dir)
}


# ===================================
# PSM, peptide and protein processing
# ===================================
library(proteoQ)
load_expts("~\\proteoQ\\examples")

# PSM data standardization
normPSM(
  # alternative to default
  group_psm_by = pep_seq_mod, 
  group_pep_by = gene, 
  annot_kinases = TRUE, 
  
  # no default and required
  fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
            "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"),
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
  impute_na = FALSE, 
  W2_bat = ~ Term["W2.BI.TMT2-W2.BI.TMT1", 
                  "W2.JHU.TMT2-W2.JHU.TMT1", 
                  "W2.PNNL.TMT2-W2.PNNL.TMT1"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = FALSE)

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

}
