# ===================================
# PSM normalization
# ===================================

## !!!require the brief working example in `?load_expts`

## additional examples
# Mascot
normPSM(
  group_psm_by = pep_seq_mod,
  group_pep_by = prot_acc,
  fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
            "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"),
  
  # variable argument statement(s)
  filter_psms_at = exprs(pep_expect <= .1),
  filter_psms_more = exprs(pep_rank == 1, pep_exp_z > 1),
)

# MaxQuant
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

# Spectrum Mill
normPSM(
  group_psm_by = pep_seq_mod,
  group_pep_by = prot_acc,
  fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
            "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"),
  
  # vararg statement(s)
  filter_psms = exprs(score >= 10),
)


## Not run: 
# wrong fasta 
normPSM(
  fasta = "~\\proteoQ\\dbs\\fasta\\wrong.fasta",
)

# no mouse entry annotation
normPSM(
  fasta = "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
)

# bad vararg statement
normPSM(
  fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
            "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"),
  filter_psms_at = exprs(column_key_not_in_psm_tables <= .1),
)
## End(Not run)
