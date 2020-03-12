\donttest{
# ===============================================
# Apply custom `human` and `mouse` Entrez lookups
# ===============================================
## A RefSeq example
# Prep I: fetch up-to-date `org.Xx.eg.db`
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")

# Prep II: set up RefSeq Fasta(s) if not yet
library(proteoQDA)
db_path <- "~\\proteoQ\\dbs\\fasta\\refseq"
copy_refseq_hs(db_path)
copy_refseq_mm(db_path)

# Prep III: copy metadata and PSMs if not yet
dat_dir <- "~\\proteoQ\\custom_refseq_lookups"
copy_global_exptsmry(dat_dir)
copy_global_fracsmry(dat_dir)

# (for MaxQuant, use `copy_global_maxquant()`)
# (for Spectrum Mill, use `copy_global_sm()`)
copy_global_mascot(dat_dir)


# --- workflow begins ---
library(proteoQ)
dat_dir <- "~\\proteoQ\\custom_refseq_lookups"
load_expts()

# prepare RefSeq-to-Entrez lookups
Ref2Entrez(species = human)
Ref2Entrez(species = mouse)

# head(readRDS(file.path("~\\proteoQ\\dbs\\entrez\\refseq_entrez_hs.rds")))
# head(readRDS(file.path("~\\proteoQ\\dbs\\entrez\\refseq_entrez_mm.rds")))

# overrule the default `Entrez` lookups with the custom databases
normPSM(
  group_psm_by = pep_seq_mod, 
  group_pep_by = gene, 
  fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
            "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"),
  entrez = c("~\\proteoQ\\dbs\\entrez\\refseq_entrez_hs.rds", 
             "~\\proteoQ\\dbs\\entrez\\refseq_entrez_mm.rds"),
)


## A UniProt example
# Prep I: set up UniProt Fasta(s) if not yet
library(proteoQDA)
db_path <- "~\\proteoQ\\dbs\\fasta\\uniprot"
copy_uniprot_hs(db_path)
copy_uniprot_mm(db_path)

# Prep II: copy metadata and PSMs if not yet
dat_dir <- "~\\proteoQ\\custom_uniprot_lookups"
copy_global_exptsmry(dat_dir)
copy_global_fracsmry(dat_dir)

# (for Mascot, use `copy_global_mascot()`)
# (for Spectrum Mill, use `copy_global_sm()`)
copy_global_maxquant(dat_dir)

# Prep III: simulate UniProt data from RefSeq PSMs
# (for Mascot, use `simulUniprotPSM(Mascot)`)
library(proteoQ)
simulUniprotPSM(MaxQuant)


# --- workflow begins ---
library(proteoQ)
dat_dir <- "~\\proteoQ\\custom_uniprot_lookups"
load_expts()

# prepare UniProt-to-Entrez lookups
Uni2Entrez(species = human)
Uni2Entrez(species = mouse)

# head(readRDS(file.path("~\\proteoQ\\dbs\\entrez\\uniprot_entrez_hs.rds")))
# head(readRDS(file.path("~\\proteoQ\\dbs\\entrez\\uniprot_entrez_mm.rds")))

normPSM(
  group_psm_by = pep_seq_mod, 
  group_pep_by = gene, 
  fasta = c("~\\proteoQ\\dbs\\fasta\\uniprot\\uniprot_hs_2014_07.fasta",
            "~\\proteoQ\\dbs\\fasta\\uniprot\\uniprot_mm_2014_07.fasta"),
  entrez = c("~\\proteoQ\\dbs\\entrez\\uniprot_entrez_hs.rds", 
             "~\\proteoQ\\dbs\\entrez\\uniprot_entrez_mm.rds"),
)
}


\dontrun{
# name your `species`
Uni2Entrez(species = this_human, abbr_species = Hs, filename = my_human.rds)
Uni2Entrez(species = this_mouse, abbr_species = Mm, filename = my_mouse.rds)

# head(readRDS(file.path("~\\proteoQ\\dbs\\entrez\\my_human.rds")))
# head(readRDS(file.path("~\\proteoQ\\dbs\\entrez\\my_mouse.rds")))

# in PSM and subsequent outputs, values under column `species` 
#  will be shown as "this_human" or "this_mouse"
normPSM(
  group_psm_by = pep_seq_mod, 
  group_pep_by = gene, 
  fasta = c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta",
            "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta"),
  entrez = c("~\\proteoQ\\dbs\\entrez\\my_human.rds", 
             "~\\proteoQ\\dbs\\entrez\\my_mouse.rds"),
)  
}


\dontrun{
## Custom database(s) are required for workflows 
#  with species other than `human`, `mouse` and `rat`
BiocManager::install("org.Ce.eg.db")
library(org.Ce.eg.db)

library(proteoQ)
Uni2Entrez(species = "worm", abbr_species = "Ce", filename = uniprot_entrez_ce.rds)

# --- PAUSE: prepare Fasta file(s) before proceeding to `normPSM`. ---

normPSM(
  fasta = "~\\proteoQ\\dbs\\fasta\\specify_your_worm.fasta",
  entrez = c("~\\proteoQ\\dbs\\entrez\\uniprot_entrez_ce.rds"),
)
}


\dontrun{
# wrong `abbr_species` provided `species` other than "human", "mouse" and "rat"
Uni2Entrez(species = "my human", abbr_species = Hu, filename = my_human.rds, overwrite = TRUE)
  
# the value of `abbr_species` ignored at `species` among "human", "mouse" and "rat"
Uni2Entrez(species = human, abbr_species = ok_not_Hs, filename = my_human.rds, overwrite = TRUE)
}
