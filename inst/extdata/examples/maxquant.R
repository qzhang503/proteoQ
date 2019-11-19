## =========================================================================================
## or MaxQuant workflow
## =========================================================================================

# special guide for the installation of data package proteoQDB containing MaxQuant examples: 
#  (1) clone https://github.com/qiangzhang503/proteoQDB.git through `Github Desktop`
#  (2) local installation of proteoQDB, e.g., devtools::install("~\\my_dir\\proteoQDB")

## part 1 --- setup
# fasta files to database directory
library(proteoQDA)
dir.create("~\\proteoQ\\dbs\\fasta\\refseq", recursive = TRUE, showWarnings = FALSE)
copy_refseq_hs("~\\proteoQ\\dbs\\fasta\\refseq")
copy_refseq_mm("~\\proteoQ\\dbs\\fasta\\refseq")

# examplary PSM data to working directory
library(proteoQDB)
dir.create("C:\\The\\MQ\\Example", recursive = TRUE, showWarnings = FALSE)
dat_dir <- c("C:\\The\\MQ\\Example")
cptac_mqpsm_txt(dat_dir)

# metadata to working directory
cptac_mqpsm_expt(dat_dir)
cptac_mqpsm_frac(dat_dir)

# metadata to workspace
library(proteoQ)
load_expts()


## part 2 --- summarization of PSMs to peptides and proteins
# PSM tables
normPSM(
	group_psm_by = pep_seq, 
	group_pep_by = gene, 
	fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta", 
					"~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"), 
	rptr_intco = 3000,					
	corrected_int = TRUE,
	rm_reverses = TRUE,
	rm_craps = TRUE,
	rm_krts = FALSE,
	rm_outliers = FALSE, 
	annot_kinases = TRUE,	
	plot_rptr_int = TRUE, 
	plot_log2FC_cv = TRUE, 
	
	filter_psms_by = exprs(PEP <= 0.1), 
)

# peptide tables
normPep(
	method_psm_pep = median, 
	method_align = MGKernel, 
	range_log2r = c(5, 95), 
	range_int = c(5, 95), 
	n_comp = 3, 
	seed = 749662, 
	maxit = 200, 
	epsilon = 1e-05, 
)

# protein tables
normPrn(
	use_unique_pep = TRUE, 
	method_pep_prn = median, 
	method_align = MGKernel, 
	range_log2r = c(20, 95), 
	range_int = c(5, 95), 
	n_comp = 2, 
	seed = 749662, 
	maxit = 200, 
	epsilon = 1e-05, 
	filter_prots_by = exprs(prot_n_pep >= 2),
)

## global setting of scaling normalization
scale_log2r = TRUE


# optional: peptide NA imputation (warning: may take a while)
pepImp(m = 2, maxit = 2)

# optional: protein NA imputation
prnImp(m = 5, maxit = 5)


