\donttest{
## Update `human` and `mouse` entrez lookups
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")

library(org.Hs.eg.db)
library(org.Mm.eg.db)

library(proteoQ)
prepEntrez(human)
prepEntrez(mouse)

# my_entrez_hs <- readRDS(file.path("~\\proteoQ\\dbs\\entrez\\uniprot_entrez_hs.rds"))
# my_entrez_mm <- readRDS(file.path("~\\proteoQ\\dbs\\entrez\\uniprot_entrez_mm.rds"))

## Not run: 
# replace the default `entrez` lookups with custom databases in PSM normalization
normPSM(
  group_psm_by = pep_seq_mod, 
  group_pep_by = gene, 
  fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
            "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"),
  entrez = c("~\\proteoQ\\dbs\\entrez\\uniprot_entrez_hs.rds", 
             "~\\proteoQ\\dbs\\entrez\\uniprot_entrez_mm.rds"),
)
## End(Not run)


## custom `species` name: 
prepEntrez(species = this_human, abbr_species = Hs, filename = my_human.rds)
prepEntrez(species = this_mouse, abbr_species = Mm, filename = my_mouse.rds)

## Not run: 
#  PSM and subsequent outputs under the column `species` 
#    will be shown as "this_human" or "this_mouse"
normPSM(
  group_psm_by = pep_seq_mod, 
  group_pep_by = gene, 
  fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
            "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"),
  entrez = c("~\\proteoQ\\dbs\\entrez\\my_human.rds", 
             "~\\proteoQ\\dbs\\entrez\\my_mouse.rds"),
)
## End(Not run)


## Custom databases required for species other than `human`, `mouse` and `rat`
BiocManager::install("org.Hs.eg.db")
library(org.Ce.eg.db)

library(proteoQ)
prepEntrez(species = "worm", abbr_species = "Ce", filename = uniprot_entrez_ce.rds)

## Not run: 
normPSM(
  fasta = "~\\proteoQ\\dbs\\fasta\\specify_your_worm.fasta",
  entrez = c("~\\proteoQ\\dbs\\entrez\\uniprot_entrez_ce.rds"),
)
## End(Not run)


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



}
