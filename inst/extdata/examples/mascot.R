## =========================================================================================
## Mascot workflow
## =========================================================================================

## part 1 --- setup
# FASTA files to database directory
library(proteoQDA)
dir.create("~\\proteoQ\\dbs\\fasta\\refseq", recursive = TRUE, showWarnings = FALSE)
copy_refseq_hs("~\\proteoQ\\dbs\\fasta\\refseq")
copy_refseq_mm("~\\proteoQ\\dbs\\fasta\\refseq")

# examplary Mascot PSM data to working directory
dir.create("C:\\The\\Mascot\\Example", recursive = TRUE, showWarnings = FALSE)
dat_dir <- "C:\\The\\Mascot\\Example"
cptac_csv_1(dat_dir)

# metadata to working directory
cptac_expt_1(dat_dir)
cptac_frac_1(dat_dir)

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
	rm_craps = TRUE,
	rm_krts = FALSE,
	rm_outliers = FALSE, 
	annot_kinases = TRUE, 
	plot_rptr_int = TRUE, 
	plot_log2FC_cv = TRUE, 
	
	filter_psms = exprs(pep_expect <= .1, pep_score >= 15), 
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
	method_pep_prn = median, 
	method_align = MGKernel, 
	range_log2r = c(20, 95), 
	range_int = c(5, 95), 
	n_comp = 2, 
	seed = 749662, 
	maxit = 200, 
	epsilon = 1e-05, 
)

## global setting of scaling normalization
scale_log2r = TRUE

# optional: peptide NA imputation (warning: may take a while)
pepImp(m = 2, maxit = 2)

# optional: protein NA imputation
prnImp(m = 5, maxit = 5)


