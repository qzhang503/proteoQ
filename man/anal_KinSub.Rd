% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kin.R
\name{anal_KinSub}
\alias{anal_KinSub}
\title{PSP outputs of kinase-substrate interactions}
\usage{
anal_KinSub(
  db_nms = "~/proteoQ/dbs/psp/Kinase_Substrate_Dataset.txt",
  type = c("peptide", "protein"),
  match_orgs = TRUE,
  scale_log2r = TRUE,
  complete_cases = FALSE,
  impute_na = FALSE,
  df = NULL,
  filepath = NULL,
  filename = NULL,
  ...
)
}
\arguments{
\item{db_nms}{Character string(s) to the name(s) of PSP database(s) with
prepended directory path(s). Users need to download the kinase-substrate
table, e.g. \code{Kinase_Substrate_Dataset.txt} directly from the PSP
website and supply the corresponding file path(s) and name(s). Currently
assume single database file.}

\item{type}{The type of data for annotation. The default is \code{peptide}.}

\item{match_orgs}{Logical; if TRUE, matches the organism between kinases and
their acting substrates. The default is TRUE.}

\item{scale_log2r}{Not currently used. Values before and after scaling will
be both reported.}

\item{complete_cases}{Logical; if TRUE, only cases that are complete with no
missing values will be used. The default is FALSE.}

\item{impute_na}{Logical; if TRUE, data with the imputation of missing values
will be used. The default is FALSE.}

\item{df}{The name of a primary data file. By default, it will be determined
automatically after matching the types of data and analysis with an
\code{id} among \code{c("pep_seq", "pep_seq_mod", "prot_acc", "gene")}. A
primary file contains normalized peptide or protein data and is among
\code{c("Peptide.txt", "Peptide_pVal.txt", "Peptide_impNA_pVal.txt",
"Protein.txt", "Protein_pVal.txt", "protein_impNA_pVal.txt")}. For analyses
require the fields of significance p-values, the \code{df} will be one of
\code{c("Peptide_pVal.txt", "Peptide_impNA_pVal.txt", "Protein_pVal.txt",
"protein_impNA_pVal.txt")}.}

\item{filepath}{Use system default.}

\item{filename}{The name of a output file.}

\item{...}{\code{filter_}: Variable argument statements for the row
filtration against data in a primary file linked to \code{df}. See also
\code{\link{normPSM}} for the format of \code{filter_} statements. \cr \cr
\code{arrange_}: Variable argument statements for the row ordering against
data in a primary file linked to \code{df}. See also \code{\link{prnHM}}
for the format of \code{arrange_} statements.}
}
\description{
\code{anal_KinSub} adds the data of \href{https://www.phosphosite.org/}{PSP}
kinase-substrate interactions to peptide or protein results.
}
\details{
OUtputs under folder \code{KinSub}.
}
\examples{
\donttest{
# ===================================
# PhosphoSitePlus (PSP)
# ===================================

## !!!require the brief working example in `?load_expts`

library(proteoQ)

# expression data for kinases and their substrates
# (need to first download the ".txt" from PSP)
anal_KinSub("~/proteoQ/dbs/psp/Kinase_Substrate_Dataset.txt")

# `human` only ('unknown' species will be removed)
anal_KinSub(
  filter_by_sp = exprs(species == "human"),
  filter_prots_by = exprs(prot_n_pep >= 2),
  filename = human_2peps.txt,
)

# proteins
anal_KinSub(type = protein)
}

\dontrun{
  # peptides: CDK1 substrates
  anal_KinSub(
    filter_by_gene = exprs(species == "human", gene == "CDK1"),
    filename = hu_CDK1.txt,
  )
  
  # proteins: CDK1 substrates
  anal_KinSub(
    type = protein,
    filter_by_gene = exprs(species == "human", gene == "CDK1"),
    filename = hu_CDK1.txt,
  )
}
}
