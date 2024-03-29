% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gspa.R
\name{pepGSPA}
\alias{pepGSPA}
\alias{prnGSPA}
\title{GSPA of peptide data}
\usage{
pepGSPA(
  gset_nms = c("go_sets", "c2_msig", "kinsub"),
  method = "mean",
  scale_log2r = TRUE,
  complete_cases = FALSE,
  impute_na = FALSE,
  pval_cutoff = 0.05,
  logFC_cutoff = log2(1.2),
  gspval_cutoff = 0.05,
  gslogFC_cutoff = log2(1.2),
  min_size = 10,
  max_size = Inf,
  min_delta = 4,
  min_greedy_size = 1,
  use_adjP = FALSE,
  fml_nms = NULL,
  df = NULL,
  filepath = NULL,
  filename = NULL,
  ...
)

prnGSPA(
  gset_nms = c("go_sets", "c2_msig", "kinsub"),
  method = "mean",
  scale_log2r = TRUE,
  complete_cases = FALSE,
  impute_na = FALSE,
  pval_cutoff = 0.05,
  logFC_cutoff = log2(1.2),
  gspval_cutoff = 0.05,
  gslogFC_cutoff = log2(1.2),
  min_size = 10,
  max_size = Inf,
  min_delta = 4,
  min_greedy_size = 1,
  use_adjP = FALSE,
  fml_nms = NULL,
  df = NULL,
  filepath = NULL,
  filename = NULL,
  ...
)
}
\arguments{
\item{gset_nms}{Character string or vector containing the shorthanded name(s),
full file path(s), or both, to gene sets for enrichment analysis. For
species among \code{"human", "mouse", "rat"}, the default of
\code{c("go_sets", "c2_msig", "kinsub")} will utilize terms from gene
ontology (\code{GO}), molecular signatures (\code{MSig}) and
kinase-substrate network (\code{PSP Kinase-Substrate}). Custom \code{GO},
\code{MSig} and other data bases at given species are also supported. See
also: \code{\link{prepGO}} for the preparation of custom \code{GO};
\code{\link{prepMSig}} for the preparation of custom \code{MSig}. For other
custom data bases, follow the same format of list as \code{GO} or
\code{MSig}.}

\item{method}{Character string or vector; the method to assess the p-values of
GSPA. The default is \code{mean} and the alternative is \code{limma}. See
also section \code{Details} for the calculations.}

\item{scale_log2r}{Logical; if TRUE, adjusts \code{log2FC} to the same scale
of standard deviation across all samples. The default is TRUE. At
\code{scale_log2r = NA}, the raw \code{log2FC} without normalization will
be used.}

\item{complete_cases}{Logical; if TRUE, only cases that are complete with no
missing values will be used. The default is FALSE.}

\item{impute_na}{Logical; if TRUE, data with the imputation of missing values
will be used. The default is FALSE.}

\item{pval_cutoff}{Numeric value or vector; the cut-off in protein
significance \code{pVal}. Entries with \code{pVals} less significant than
the threshold will be excluded from enrichment analysis. The default is 0.05
for all formulas matched to or specified in argument \code{fml_nms}.
Formula-specific threshold is allowed by supplying a vector of cut-off
values.}

\item{logFC_cutoff}{Numeric value or vector; the cut-off in protein
\code{log2FC}. Entries with absolute \code{log2FC} smaller than the
threshold will be excluded from enrichment analysis. The default magnitude
is \code{log2(1.2)} for all formulas matched to or specified in argument
\code{fml_nms}. Formula-specific threshold is allowed by supplying a vector
of absolute values in \code{log2FC}.}

\item{gspval_cutoff}{Numeric value or vector; the cut-off in gene-set
significance \code{pVal}. Only enrichment terms with \code{pVals} more
significant than the threshold will be reported. The default is 0.05 for all
formulas matched to or specified in argument \code{fml_nms}.
Formula-specific threshold is allowed by supplying a vector of cut-off
values.}

\item{gslogFC_cutoff}{Numeric value or vector; the cut-off in gene-set
enrichment fold change. Only enrichment terms with absolute fold change
greater than the threshold will be reported. The default magnitude is
\code{log2(1.2)} for all formulas matched to or specified in argument
\code{fml_nms}. Formula-specific threshold is allowed by supplying a vector
of absolute values in \code{log2FC}.}

\item{min_size}{Numeric value or vector; minimum number of protein entries for
consideration in gene set tests. The number is after data filtration by
\code{pval_cutoff}, \code{logFC_cutoff} or varargs expressions under
\code{filter_}. The default is 10 for all formulas matched to or specified
in argument \code{fml_nms}. Formula-specific threshold is allowed by
supplying a vector of sizes.}

\item{max_size}{Numeric value or vector; maximum number of protein entries for
consideration in gene set tests. The number is after data filtration by
\code{pval_cutoff}, \code{logFC_cutoff} or varargs expressions under
\code{filter_}. The default in infinite for all formulas matched to or
specified in argument \code{fml_nms}. Formula-specific threshold is allowed
by supplying a vector of sizes.}

\item{min_delta}{Numeric value or vector; the minimum count difference between
the up- and the down-expressed group of proteins for consideration in gene
set tests. For example at \code{min_delta = 4}, a gene set will 6
upregulated proteins and 2 down-expressed proteins, or vice versa, will be
assessed. The number is after data filtration by \code{pval_cutoff},
\code{logFC_cutoff} or varargs expressions under \code{filter_}. The default
is 4 for all formulas matched to or specified in argument \code{fml_nms}.
Formula-specific threshold is allowed by supplying a vector of sizes.}

\item{min_greedy_size}{Numeric value or vector; minimum number of unique
protein entries for a gene set to be considered essential. The default in
\code{1} for all formulas matched to or specified in argument
\code{fml_nms}. Formula-specific threshold is allowed by supplying a vector
of sizes.}

\item{use_adjP}{Logical; if TRUE, use Benjamini-Hochberg pVals. The default is
FALSE.}

\item{fml_nms}{Character string or vector; the formula name(s). By default,
the formula(s) will match those used in \code{\link{pepSig}} or
\code{\link{prnSig}}.}

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

\item{filename}{A representative file name to outputs. By default, it will be
determined automatically by the name of the current call.}

\item{...}{\code{filter_}: Logical expression(s) for the row filtration
against data in a primary file of \code{/Model/Protein[_impNA]_pVals.txt}.
See also \code{\link{normPSM}} for the format of \code{filter_} statements.
\cr \cr \code{arrange_}: Variable argument statements for the row ordering
against data in a primary file of \code{/Model/Protein[_impNA]_pVals.txt}.
See also \code{\link{prnHM}} for the format of \code{arrange_} statements.}
}
\description{
\code{pepGSPA} performs the analysis of Gene Set Probability Asymmetricity
(GSPA) against peptide \code{log2FC} data.

\code{prnGSPA} performs the analysis of Gene Set Probability Asymmetricity
(GSPA) against protein \code{log2FC} data.
}
\details{
The significance \code{pVals} of individual proteins are first obtained from
\code{\link{prnSig}}, followed by log10 transformation and separation into up-
or down-expressed groups for each gene set. At the default of \code{method =
mean}, the geometric mean of \code{pVals}, \eqn{P}, are each calculated for
the groups of up or down regulated proteins, with a penalty-like term

\deqn{-log10(P)=(\sum_{i=1}^{n}-log10(p)+m)/(n+m)}{-log10(P)=(-log10(p_1*p_2*...*p_n)+m)/(n+m)}

where \eqn{n} and \eqn{m} are the numbers of entries with \eqn{p} values
\eqn{\le} or \eqn{>} a significance cut-off, respectively. The quotient of the
two \eqn{P} values is then used to represent the significance of gene set
enrichment. The arguments \code{pval_cutoff} and \code{logFC_cutoff} are used
to discriminate low influence genes. Additional subsetting of data via the
\code{vararg} approach of \code{filter_} is feasible. At \code{method =
limma}, moderated t-tests are performed against \code{-log10(pVals)} between
the up and the down groups via \code{\link[limma]{eBayes}}.
}
\examples{
\donttest{
# ===================================
# GSPA tests
# ===================================

## !!!require the brief working example in `?load_expts`

## global option
scale_log2r <- TRUE

## base
# prerequisites in significance tests (impute_na = FALSE)
# as `pVals` are required
pepSig(
  impute_na = FALSE, 
  W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", 
                  "(W2.JHU.TMT2-W2.JHU.TMT1)", 
                  "(W2.PNNL.TMT2-W2.PNNL.TMT1)"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = FALSE)

# GSPA analysis
prnGSPA(
  impute_na = FALSE,
  pval_cutoff = 5E-2, # protein pVal threshold
  logFC_cutoff = log2(1.2), # protein log2FC threshold
  gspval_cutoff = 5E-2, # gene-set pVal threshold
  gslogFC_cutoff = log2(1.2), # gene-set log2FC threshold
  gset_nms = c("go_sets", "c2_msig"), 
)

# GSPA visualization under volcano plots 
# forced lines of `pval_cutoff` and `logFC_cutoff`  
#   from the corresponding `prnGSPA`;
# more examples of aesthetics in `?gspaMap`
gspaMap(
  impute_na = FALSE,
  topn_gsets = 50, 
  show_sig = pVal, 
  xco = 1.2, # optional x-axis cut-off lines
  yco = 0.05, # a optional y-axis cut-off line 
)

## row filtration
prnGSPA(
  impute_na = FALSE,
  gset_nms = c("go_sets", "c2_msig"),
  filter_prots_by_npep = exprs(prot_n_pep >= 2),
  filter_prots_by_species = exprs(species == "human"), 
  filename = row_filters.txt,
)

# by default, volcano-plot visualization for all GSPA result files, 
#   which are "Protein_GSPA_Z.txt" and "row_filters_Protein_GSPA_Z.txt" 
#   up to this point
# will match the primary `fliter_` varargs to 
#   those in the corresponding `prnGSPA(...)`
gspaMap(
  impute_na = FALSE,
  topn_labels = 0, 
  topn_gsets = 50, 
  show_sig = pVal, 
)

# or manually supply specific secondary file(s) with `df2`; 
gspaMap(
  df2 = "row_filters_Protein_GSPA_Z.txt", 
  impute_na = FALSE,
  topn_gsets = 10, 
  show_sig = pVal, 
)

# may apply additionally secondary `fliter2_` varargs against `df2`
gspaMap(
  df2 = "row_filters_Protein_GSPA_Z.txt", 
  filter2_ess_size = exprs(ess_size >= 1),   
  impute_na = FALSE,
  topn_gsets = 10, 
  topn_labels = 0, 
  show_sig = pVal, 
)

# can also visualize results for specific formula(e) and/or 
#   specific `df2`
gspaMap(
  fml_nms = "W16_vs_W2",
  df2 = "row_filters_Protein_GSPA_Z.txt", 
  filter2_ess_size = exprs(ess_size >= 1),   
  impute_na = FALSE,
  topn_gsets = 20, 
  show_sig = pVal, 
)

## customized thresholds
# (formula names defined a priori in `pepSig(...)`)
# (higher thresholds for "W16_vs_W2"; also set `max_size` to 
#  avoid trapping with large terms, such as "cell parts"...)
prnGSPA(
  impute_na = FALSE, 
  fml_nms = c("W2_bat", "W2_loc", "W16_vs_W2"),
  pval_cutoff = c(5E-2, 5E-2, 1E-6), 
  logFC_cutoff = c(log2(1.2)), 
  gslogFC_cutoff = c(log2(1.2)), 
  gspval_cutoff = c(5E-2, 5E-2, 1E-4), 
  max_size = c(Inf, Inf, 120), 
  gset_nms = c("go_sets"), 
  filter_by_npep = exprs(prot_n_pep >= 2), 
  filter_prots_by_species = exprs(species == "human"), 
  filename = thresholds.txt,
)

gspaMap(
  df2 = "thresholds_Protein_GSPA_Z.txt", 
  impute_na = FALSE,
  topn_gsets = 10, 
  topn_labels = 0, 
  show_sig = pVal, 
)

# "limma" for "W16_vs_W2"
prnGSPA(
  impute_na = FALSE, 
  fml_nms = c("W2_bat", "W2_loc", "W16_vs_W2"),
  method = c("mean", "mean", "limma"),
  max_size = c(Inf, Inf, 120), 
  pval_cutoff = c(5E-2), 
  logFC_cutoff = c(log2(1.2)), 
  gspval_cutoff = c(5E-2), 
  gset_nms = c("go_sets", "c2_msig"), 
  filter_by = exprs(prot_n_pep >= 2, species == "human"), 
  filename = methods.txt,
)

gspaMap(
  df2 = "methods_Protein_GSPA_Z.txt", 
  impute_na = FALSE,
  topn_gsets = 25, 
  topn_labels = 0, 
  show_sig = pVal, 
)

## NA imputation
# if not yet, run prerequisitive NA imputation
pepImp(m = 2, maxit = 2)
prnImp(m = 5, maxit = 5)

# if not yet, run prerequisitive significance tests at `impute_na = TRUE`
pepSig(
  impute_na = TRUE, 
  W2_bat = ~ Term["(W2.BI.TMT2-W2.BI.TMT1)", 
                  "(W2.JHU.TMT2-W2.JHU.TMT1)", 
                  "(W2.PNNL.TMT2-W2.PNNL.TMT1)"],
  W2_loc = ~ Term_2["W2.BI-W2.JHU", 
                    "W2.BI-W2.PNNL", 
                    "W2.JHU-W2.PNNL"],
  W16_vs_W2 = ~ Term_3["W16-W2"], 
)

prnSig(impute_na = TRUE)

prnGSPA(
  impute_na = TRUE,
  gset_nms = c("go_sets"),
  filter_prots_by_npep = exprs(prot_n_pep >= 2),
)

# ===================================
# Distance heat maps of gene sets
# (also interactive networks)
# ===================================
# a `term` is a subset of an `ess_term` if the distance is zero
prnGSPAHM(
  filter2_by = exprs(distance <= .7),
  annot_cols = "ess_idx",
  annot_colnames = "Eset index",
  annot_rows = "ess_size",
  filename = show_some_redundancy.png,
)

# human terms only
prnGSPAHM(
  filter2_num = exprs(distance <= .95),
  filter2_sp = exprs(start_with_str("hs", term)),
  annot_cols = "ess_idx",
  annot_colnames = "Eset index",
  filename = show_more_connectivity.png,
)

# custom color palette
prnGSPAHM(
  annot_cols = c("ess_idx", "ess_size"),
  annot_colnames = c("Eset index", "Size"),
  filter2_by = exprs(distance <= .95),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  filename = custom_colors.png,
)

# can also visualize results for specific `df2` and formula(e)
prnGSPAHM(
  df2 = "thresholds_Protein_GSPA_Z_essmap.txt", 
  fml_nms = "W16_vs_W2", 
  annot_cols = c("ess_idx", "ess_size"),
  annot_colnames = c("Eset index", "Size"),
  filter2_by = exprs(distance <= .95),
  filename = w16_vs_w2.png,
)


# ===================================
# custom GO databases
# ===================================
# start over
unlink(file.path(dat_dir, "Protein/GSPA"), recursive = TRUE, force = TRUE)

# prepare GO databases (see also ?prepGO)
prepGO(
  species = human,
  db_path = "~/proteoQ/dbs/go",
  gaf_url = "http://current.geneontology.org/annotations/goa_human.gaf.gz",
  obo_url = "http://purl.obolibrary.org/obo/go/go-basic.obo",
  filename = go_hs.rds,
)

prepGO(
  species = mouse,
  db_path = "~/proteoQ/dbs/go",
  gaf_url = "http://current.geneontology.org/annotations/mgi.gaf.gz",
  obo_url = "http://purl.obolibrary.org/obo/go/go-basic.obo",
  filename = go_mm.rds,
)

head(readRDS(file.path("~/proteoQ/dbs/go", "go_hs.rds")))
head(readRDS(file.path("~/proteoQ/dbs/go", "go_mm.rds")))

# analysis
prnGSPA(
  impute_na = FALSE,
  gset_nms = c("~/proteoQ/dbs/go/go_hs.rds",
               "~/proteoQ/dbs/go/go_mm.rds"),
)

# It is possible that the system "go_sets" may be partially different to the
# custom GO. Be explicit on `gset_nms`; otherwise, the default of `gset_nms =
# c("go_sets", "c2_msig")` will be applied.
gspaMap(
  gset_nms = c("~/proteoQ/dbs/go/go_hs.rds",
               "~/proteoQ/dbs/go/go_mm.rds"),
  impute_na = FALSE,
  topn_gsets = 20, 
  show_sig = pVal, 
)

# ===================================
# custom MSig databases
# ===================================
# start over
unlink(file.path(dat_dir, "Protein/GSPA"), recursive = TRUE, force = TRUE)

# prepare MSig databases (see also ?prepMSig)
prepMSig(species = human)
prepMSig(species = mouse)

head(readRDS(file.path("~/proteoQ/dbs/msig", "msig_hs.rds")))
head(readRDS(file.path("~/proteoQ/dbs/msig", "msig_mm.rds")))

# analysis
prnGSPA(
  impute_na = FALSE,
  gset_nms = c("~/proteoQ/dbs/msig/msig_hs.rds",
               "~/proteoQ/dbs/msig/msig_mm.rds"),
)

# It is possible that the system "c2_msig" may be partially different to the
# custom databases. Be explicit on `gset_nms`; otherwise, the default of
# `gset_nms = c("go_sets", "c2_msig")` will be applied.
gspaMap(
  gset_nms = c("~/proteoQ/dbs/msig/msig_hs.rds",
               "~/proteoQ/dbs/msig/msig_mm.rds"),
  impute_na = FALSE,
  topn_gsets = 20, 
  show_sig = pVal, 
)

## put custom `prepGO`, `prepMSig` preparation into `prnGSPA` analysis
# start over with fresh database preparation
unlink(file.path("~/proteoQ/dbs/go"), recursive = TRUE, force = TRUE)
unlink(file.path("~/proteoQ/dbs/msig"), recursive = TRUE, force = TRUE)

library(magrittr)

prnGSPA(
  impute_na = FALSE,
  pval_cutoff = 5E-2,
  logFC_cutoff = log2(1.2),
  gspval_cutoff = 5E-2,
  gset_nms = c(prepGO(human), prepGO(mouse), prepMSig(human), prepMSig(mouse)),
  filename = cmbn.txt,
) \%>\% 
  gspaMap(
    # gset_nms = ., 
    df2 = "cmbn_Protein_GSPA_Z.txt", 
    impute_na = FALSE,
    topn_gsets = 20, 
    topn_labels = 0, 
    show_sig = pVal, 
  )
}


}
\seealso{
\emph{Metadata} \cr 
 \code{\link{load_expts}} for metadata preparation and a reduced working example in data normalization \cr

 \emph{Data normalization} \cr 
 \code{\link{normPSM}} for extended examples in PSM data normalization \cr
 \code{\link{PSM2Pep}} for extended examples in PSM to peptide summarization \cr 
 \code{\link{mergePep}} for extended examples in peptide data merging \cr 
 \code{\link{standPep}} for extended examples in peptide data normalization \cr
 \code{\link{Pep2Prn}} for extended examples in peptide to protein summarization \cr
 \code{\link{standPrn}} for extended examples in protein data normalization. \cr 
 \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples in data purging \cr
 \code{\link{pepHist}} and \code{\link{prnHist}} for extended examples in histogram visualization. \cr 
 \code{\link{extract_raws}} and \code{\link{extract_psm_raws}} for extracting MS file names \cr 
 
 \emph{Variable arguments of `filter_...`} \cr 
 \code{\link{contain_str}}, \code{\link{contain_chars_in}}, \code{\link{not_contain_str}}, 
 \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}}, 
 \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and 
 \code{\link{ends_with_chars_in}} for data subsetting by character strings \cr 
 
 \emph{Missing values} \cr 
 \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation \cr 
 
 \emph{Informatics} \cr 
 \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr 
 \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization \cr 
 \code{\link{prnGSPA}} for gene set enrichment analysis by protein significance pVals \cr 
 \code{\link{gspaMap}} for mapping GSPA to volcano plot visualization \cr 
 \code{\link{prnGSPAHM}} for heat map and network visualization of GSPA results \cr 
 \code{\link{prnGSVA}} for gene set variance analysis \cr 
 \code{\link{prnGSEA}} for data preparation for online GSEA. \cr 
 \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS visualization \cr 
 \code{\link{pepPCA}} and \code{\link{prnPCA}} for PCA visualization \cr 
 \code{\link{pepLDA}} and \code{\link{prnLDA}} for LDA visualization \cr 
 \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
 \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
 \code{\link{prnCorr_logInt}}  for correlation plots \cr 
 \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend analysis and visualization \cr 
 \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
 \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
 \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr 
 
 \emph{Custom databases} \cr 
 \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
 \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
 \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
 ontology}} \cr 
 \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
 signatures} \cr 
 \code{\link{prepString}} and \code{\link{anal_prnString}} for STRING-DB \cr
 
 \emph{Column keys in PSM, peptide and protein outputs} \cr 
 system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
 system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
 system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
}
