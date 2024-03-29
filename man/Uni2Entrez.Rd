% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/go.R
\name{Uni2Entrez}
\alias{Uni2Entrez}
\alias{Ref2Entrez}
\title{Map UniProt accessions to Entrez IDs}
\usage{
Uni2Entrez(
  species = "human",
  abbr_species = NULL,
  filename = NULL,
  db_path = "~/proteoQ/dbs/entrez",
  overwrite = FALSE
)

Ref2Entrez(
  species = "human",
  abbr_species = NULL,
  filename = NULL,
  db_path = "~/proteoQ/dbs/entrez",
  overwrite = FALSE
)
}
\arguments{
\item{species}{Character string; the name of a species.}

\item{abbr_species}{Two-letter character string; the abbreviated name of
species used with
\href{https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html}{org.Xx.eg.db}.
The value of \code{abbr_species} will be determined automatically if the
species is in one of \code{c("human", "mouse", "rat")}. Otherwise, for
example, users need to provide \code{abbr_species = Ce} for fetching the
\code{org.Ce.eg.db} package in the name space of proteoQ. 

For analysis against
\href{http://current.geneontology.org/products/pages/downloads.html}{gene
ontology} and \href{https://www.gsea-msigdb.org/gsea/index.jsp}{Molecular
Signatures}, the argument is further applied to differentiate the same
biological terms under different species; e.g., \code{GO~0072686 mitotic
spindle} becomes \code{hs_GO~0072686 mitotic spindle} for human and
\code{mm_GO~0072686 mitotic spindle} for mouse.}

\item{filename}{Character string; An output file name. At the \code{NULL}
default, the name will be determined automatically at a given
\code{species}; i.e., \code{msig_hs.rds} for \code{human} data. The file is
saved as a \code{.rds} object for uses with \code{\link{prnGSPA}}.}

\item{db_path}{Character string; the local path for database(s). The default
is \code{"~/proteoQ/dbs/entrez"}.}

\item{overwrite}{Logical; if TRUE, overwrite the downloaded database(s). The
default is FALSE.}
}
\description{
\code{Uni2Entrez} prepares lookup tables between UniProt accessions and
Entrez IDs for uses with \link{normPSM} and downstream gene-set analysis such
as \link{prnGSPA}. The utility is optional for \code{human}, \code{mouse} and
\code{rat} data. It is \strong{required} for other species with \link{prnGSPA}
in users' workflows. It can also be used to update and overrule the lookups
for \code{human}, \code{mouse} and \code{rat} that are defaulted by
\code{proteoQ}.

\code{Ref2Entrez} prepares lookup tables between RefSeq accessions and
Entrez IDs and gene names for uses with \link{normPSM} and downstream gene-set
analysis such as \link{prnGSPA}. The utility is optional for \code{human} and
\code{mouse} data. It is \strong{required} for other species with
\link{prnGSPA} in users' workflows. It can also be used to update and overrule
the lookups for \code{human} and \code{mouse} that are defaulted by
\code{proteoQ}.
}
\examples{
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

# Prep II: make available RefSeq Fasta(s) if not yet
library(proteoQDA)
db_path <- "~/proteoQ/dbs/fasta/refseq"
copy_refseq_hs(db_path)
copy_refseq_mm(db_path)

# Prep III: copy metadata and PSMs if not yet
dat_dir <- "~/proteoQ/custom_refseq_lookups"
copy_global_exptsmry(dat_dir)
copy_global_fracsmry(dat_dir)

# (for MaxQuant, use `copy_global_maxquant()`)
# (for Spectrum Mill, use `copy_global_sm()`)
copy_global_mascot(dat_dir)


# --- workflow begins ---
library(proteoQ)
load_expts("~/proteoQ/custom_refseq_lookups")

# prepare RefSeq-to-Entrez lookups
Ref2Entrez(species = human)
Ref2Entrez(species = mouse)

# head(readRDS(file.path("~/proteoQ/dbs/entrez/refseq_entrez_hs.rds")))
# head(readRDS(file.path("~/proteoQ/dbs/entrez/refseq_entrez_mm.rds")))

# overrule the default `Entrez` lookups with the custom databases
normPSM(
  group_psm_by = pep_seq_mod, 
  group_pep_by = gene, 
  fasta = c("~/proteoQ/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
            "~/proteoQ/dbs/fasta/refseq/refseq_mm_2013_07.fasta"),
  entrez = c("~/proteoQ/dbs/entrez/refseq_entrez_hs.rds", 
             "~/proteoQ/dbs/entrez/refseq_entrez_mm.rds"),
)


## A UniProt example
# Prep I: set up UniProt Fasta(s) if not yet
library(proteoQDA)
db_path <- "~/proteoQ/dbs/fasta/uniprot"
copy_uniprot_hs(db_path)
copy_uniprot_mm(db_path)

# Prep II: copy metadata and PSMs if not yet
dat_dir <- "~/proteoQ/custom_uniprot_lookups"
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
dat_dir <- "~/proteoQ/custom_uniprot_lookups"
load_expts()

# prepare UniProt-to-Entrez lookups
Uni2Entrez(species = human)
Uni2Entrez(species = mouse)

# head(readRDS(file.path("~/proteoQ/dbs/entrez/uniprot_entrez_hs.rds")))
# head(readRDS(file.path("~/proteoQ/dbs/entrez/uniprot_entrez_mm.rds")))

normPSM(
  group_psm_by = pep_seq_mod, 
  group_pep_by = gene, 
  fasta = c("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2014_07.fasta",
            "~/proteoQ/dbs/fasta/uniprot/uniprot_mm_2014_07.fasta"),
  entrez = c("~/proteoQ/dbs/entrez/uniprot_entrez_hs.rds", 
             "~/proteoQ/dbs/entrez/uniprot_entrez_mm.rds"),
)
}


\dontrun{
# name your `species`
Uni2Entrez(species = this_human, abbr_species = Hs, filename = my_human.rds)
Uni2Entrez(species = this_mouse, abbr_species = Mm, filename = my_mouse.rds)

# head(readRDS(file.path("~/proteoQ/dbs/entrez/my_human.rds")))
# head(readRDS(file.path("~/proteoQ/dbs/entrez/my_mouse.rds")))

# in PSM and subsequent outputs, values under column `species` 
#  will be shown as "this_human" or "this_mouse"
normPSM(
  group_psm_by = pep_seq_mod, 
  group_pep_by = gene, 
  fasta = c("~/proteoQ/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
            "~/proteoQ/dbs/fasta/refseq/refseq_mm_2013_07.fasta"),
  entrez = c("~/proteoQ/dbs/entrez/my_human.rds", 
             "~/proteoQ/dbs/entrez/my_mouse.rds"),
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
  fasta = "~/proteoQ/dbs/fasta/specify_your_worm.fasta",
  entrez = c("~/proteoQ/dbs/entrez/uniprot_entrez_ce.rds"),
)
}


\dontrun{
# wrong `abbr_species` provided `species` other than "human", "mouse" and "rat"
Uni2Entrez(species = "my human", abbr_species = Hu, filename = my_human.rds, overwrite = TRUE)
  
# the value of `abbr_species` ignored at `species` among "human", "mouse" and "rat"
Uni2Entrez(species = human, abbr_species = ok_not_Hs, filename = my_human.rds, overwrite = TRUE)
}
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

# Prep II: make available RefSeq Fasta(s) if not yet
library(proteoQDA)
db_path <- "~/proteoQ/dbs/fasta/refseq"
copy_refseq_hs(db_path)
copy_refseq_mm(db_path)

# Prep III: copy metadata and PSMs if not yet
dat_dir <- "~/proteoQ/custom_refseq_lookups"
copy_global_exptsmry(dat_dir)
copy_global_fracsmry(dat_dir)

# (for MaxQuant, use `copy_global_maxquant()`)
# (for Spectrum Mill, use `copy_global_sm()`)
copy_global_mascot(dat_dir)


# --- workflow begins ---
library(proteoQ)
load_expts("~/proteoQ/custom_refseq_lookups")

# prepare RefSeq-to-Entrez lookups
Ref2Entrez(species = human)
Ref2Entrez(species = mouse)

# head(readRDS(file.path("~/proteoQ/dbs/entrez/refseq_entrez_hs.rds")))
# head(readRDS(file.path("~/proteoQ/dbs/entrez/refseq_entrez_mm.rds")))

# overrule the default `Entrez` lookups with the custom databases
normPSM(
  group_psm_by = pep_seq_mod, 
  group_pep_by = gene, 
  fasta = c("~/proteoQ/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
            "~/proteoQ/dbs/fasta/refseq/refseq_mm_2013_07.fasta"),
  entrez = c("~/proteoQ/dbs/entrez/refseq_entrez_hs.rds", 
             "~/proteoQ/dbs/entrez/refseq_entrez_mm.rds"),
)


## A UniProt example
# Prep I: set up UniProt Fasta(s) if not yet
library(proteoQDA)
db_path <- "~/proteoQ/dbs/fasta/uniprot"
copy_uniprot_hs(db_path)
copy_uniprot_mm(db_path)

# Prep II: copy metadata and PSMs if not yet
dat_dir <- "~/proteoQ/custom_uniprot_lookups"
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
dat_dir <- "~/proteoQ/custom_uniprot_lookups"
load_expts()

# prepare UniProt-to-Entrez lookups
Uni2Entrez(species = human)
Uni2Entrez(species = mouse)

# head(readRDS(file.path("~/proteoQ/dbs/entrez/uniprot_entrez_hs.rds")))
# head(readRDS(file.path("~/proteoQ/dbs/entrez/uniprot_entrez_mm.rds")))

normPSM(
  group_psm_by = pep_seq_mod, 
  group_pep_by = gene, 
  fasta = c("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2014_07.fasta",
            "~/proteoQ/dbs/fasta/uniprot/uniprot_mm_2014_07.fasta"),
  entrez = c("~/proteoQ/dbs/entrez/uniprot_entrez_hs.rds", 
             "~/proteoQ/dbs/entrez/uniprot_entrez_mm.rds"),
)
}


\dontrun{
# name your `species`
Uni2Entrez(species = this_human, abbr_species = Hs, filename = my_human.rds)
Uni2Entrez(species = this_mouse, abbr_species = Mm, filename = my_mouse.rds)

# head(readRDS(file.path("~/proteoQ/dbs/entrez/my_human.rds")))
# head(readRDS(file.path("~/proteoQ/dbs/entrez/my_mouse.rds")))

# in PSM and subsequent outputs, values under column `species` 
#  will be shown as "this_human" or "this_mouse"
normPSM(
  group_psm_by = pep_seq_mod, 
  group_pep_by = gene, 
  fasta = c("~/proteoQ/dbs/fasta/refseq/refseq_hs_2013_07.fasta",
            "~/proteoQ/dbs/fasta/refseq/refseq_mm_2013_07.fasta"),
  entrez = c("~/proteoQ/dbs/entrez/my_human.rds", 
             "~/proteoQ/dbs/entrez/my_mouse.rds"),
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
  fasta = "~/proteoQ/dbs/fasta/specify_your_worm.fasta",
  entrez = c("~/proteoQ/dbs/entrez/uniprot_entrez_ce.rds"),
)
}


\dontrun{
# wrong `abbr_species` provided `species` other than "human", "mouse" and "rat"
Uni2Entrez(species = "my human", abbr_species = Hu, filename = my_human.rds, overwrite = TRUE)
  
# the value of `abbr_species` ignored at `species` among "human", "mouse" and "rat"
Uni2Entrez(species = human, abbr_species = ok_not_Hs, filename = my_human.rds, overwrite = TRUE)
}
}
