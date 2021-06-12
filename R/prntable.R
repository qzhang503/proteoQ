#'Standardize protein results
#'
#'\code{standPrn} standardizes protein results from \code{\link{Pep2Prn}} with
#'additional choices in data alignment. The utility further supports iterative
#'normalization against data under selected sample columns, data rows or both.
#'
#'In the primary output file, "\code{Protein.txt}", values under columns
#'\code{log2_R...} are logarithmic ratios at base 2 in relative to the
#'\code{reference(s)} within each multiplex TMT set, or to the row means within
#'each plex if no \code{reference(s)} are present. Values under columns
#'\code{N_log2_R...} are aligned \code{log2_R...} according to
#'\code{method_align} without scaling normalization. Values under columns
#'\code{Z_log2_R...} are \code{N_log2_R...} with additional scaling
#'normalization. Values under columns \code{I...} are reporter-ion or LFQ
#'intensity before normalization. Values under columns \code{N_I...} are
#'normalized \code{I...}. Values under columns \code{sd_log2_R...} are the
#'standard deviation of the \code{log2FC} of proteins from ascribing peptides.
#'
#'@param cache Not currently used.
#'@param ... \code{slice_}: Variable argument statements for the identification
#'  of row subsets. The partial data will be taken for parameterizing the
#'  alignment of \code{log2FC} across samples. The full data set will be updated
#'  subsequently with the newly derived parameters. Note that there is no data
#'  entry removals from the complete data set with the \code{slice_} procedure.
#'  \cr \cr The variable argument statements should be in the following format:
#'  each of the statement contains a list of logical expression(s). The
#'  \code{lhs} needs to start with \code{slice_}. The logical condition(s) at
#'  the \code{rhs} needs to be enclosed in \code{exprs} with round parenthesis.
#'  For example, \code{prot_n_pep} is a column key present in
#'  \code{Protein.txt}. The \code{slice_prns_at = exprs(prot_n_pep >= 5)} will
#'  extract protein entries with five or more identifying peptide sequences for
#'  \code{log2FC} alignment. Protein entries with less than five identifying
#'  sequences will remain in \code{Protein.txt} but not used in the
#'  parameterization. See also \code{\link{normPSM}} for the variable arguments
#'  of \code{filter_}. \cr \cr Additional parameters from
#'  \code{\link[mixtools]{normalmixEM}}, i.e., \cr \code{maxit}, the maximum
#'  number of iterations allowed; \cr \code{epsilon}, tolerance limit for
#'  declaring algorithm convergence.
#'@inheritParams normPSM
#'@inheritParams standPep
#'@inheritParams mergePep
#'@seealso 
#'  \emph{Metadata} \cr 
#'  \code{\link{load_expts}} for metadata preparation and a reduced working example in data normalization \cr
#'
#'  \emph{Data normalization} \cr 
#'  \code{\link{normPSM}} for extended examples in PSM data normalization \cr
#'  \code{\link{PSM2Pep}} for extended examples in PSM to peptide summarization \cr 
#'  \code{\link{mergePep}} for extended examples in peptide data merging \cr 
#'  \code{\link{standPep}} for extended examples in peptide data normalization \cr
#'  \code{\link{Pep2Prn}} for extended examples in peptide to protein summarization \cr
#'  \code{\link{standPrn}} for extended examples in protein data normalization. \cr 
#'  \code{\link{purgePSM}} and \code{\link{purgePep}} for extended examples in data purging \cr
#'  \code{\link{pepHist}} and \code{\link{prnHist}} for extended examples in histogram visualization. \cr 
#'  \code{\link{extract_raws}} and \code{\link{extract_psm_raws}} for extracting MS file names \cr 
#'  
#'  \emph{Variable arguments of `filter_...`} \cr 
#'  \code{\link{contain_str}}, \code{\link{contain_chars_in}}, \code{\link{not_contain_str}}, 
#'  \code{\link{not_contain_chars_in}}, \code{\link{start_with_str}}, 
#'  \code{\link{end_with_str}}, \code{\link{start_with_chars_in}} and 
#'  \code{\link{ends_with_chars_in}} for data subsetting by character strings \cr 
#'  
#'  \emph{Missing values} \cr 
#'  \code{\link{pepImp}} and \code{\link{prnImp}} for missing value imputation \cr 
#'  
#'  \emph{Informatics} \cr 
#'  \code{\link{pepSig}} and \code{\link{prnSig}} for significance tests \cr 
#'  \code{\link{pepVol}} and \code{\link{prnVol}} for volcano plot visualization \cr 
#'  \code{\link{prnGSPA}} for gene set enrichment analysis by protein significance pVals \cr 
#'  \code{\link{gspaMap}} for mapping GSPA to volcano plot visualization \cr 
#'  \code{\link{prnGSPAHM}} for heat map and network visualization of GSPA results \cr 
#'  \code{\link{prnGSVA}} for gene set variance analysis \cr 
#'  \code{\link{prnGSEA}} for data preparation for online GSEA. \cr 
#'  \code{\link{pepMDS}} and \code{\link{prnMDS}} for MDS visualization \cr 
#'  \code{\link{pepPCA}} and \code{\link{prnPCA}} for PCA visualization \cr 
#'  \code{\link{pepLDA}} and \code{\link{prnLDA}} for LDA visualization \cr 
#'  \code{\link{pepHM}} and \code{\link{prnHM}} for heat map visualization \cr 
#'  \code{\link{pepCorr_logFC}}, \code{\link{prnCorr_logFC}}, \code{\link{pepCorr_logInt}} and 
#'  \code{\link{prnCorr_logInt}}  for correlation plots \cr 
#'  \code{\link{anal_prnTrend}} and \code{\link{plot_prnTrend}} for trend analysis and visualization \cr 
#'  \code{\link{anal_pepNMF}}, \code{\link{anal_prnNMF}}, \code{\link{plot_pepNMFCon}}, 
#'  \code{\link{plot_prnNMFCon}}, \code{\link{plot_pepNMFCoef}}, \code{\link{plot_prnNMFCoef}} and 
#'  \code{\link{plot_metaNMF}} for NMF analysis and visualization \cr 
#'  
#'  \emph{Custom databases} \cr 
#'  \code{\link{Uni2Entrez}} for lookups between UniProt accessions and Entrez IDs \cr 
#'  \code{\link{Ref2Entrez}} for lookups among RefSeq accessions, gene names and Entrez IDs \cr 
#'  \code{\link{prepGO}} for \code{\href{http://current.geneontology.org/products/pages/downloads.html}{gene 
#'  ontology}} \cr 
#'  \code{\link{prepMSig}} for \href{https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.0/}{molecular 
#'  signatures} \cr 
#'  \code{\link{prepString}} and \code{\link{anal_prnString}} for STRING-DB \cr
#'  
#'  \emph{Column keys in PSM, peptide and protein outputs} \cr 
#'  system.file("extdata", "psm_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "peptide_keys.txt", package = "proteoQ") \cr
#'  system.file("extdata", "protein_keys.txt", package = "proteoQ") \cr
#'
#'@return The primary output is in \code{.../Protein/Protein.txt}.
#'
#'@example inst/extdata/examples/normPrn_.R
#'@import stringr dplyr tidyr purrr
#'@importFrom magrittr %>% %T>% %$% %<>%
#'@importFrom plyr ddply
#'@export
standPrn <- function (method_align = c("MC", "MGKernel"), 
                      range_log2r = c(10, 90), range_int = c(5, 95), 
                      n_comp = NULL, seed = NULL, 
                      col_select = NULL, cut_points = Inf, cache = TRUE, ...) {
  dat_dir <- get_gl_dat_dir()
  dir.create(file.path(dat_dir, "Protein/Histogram"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Protein/cache"), 
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dat_dir, "Protein/log"), 
             recursive = TRUE, showWarnings = FALSE)
  
  old_opts <- options()
  on.exit(options(old_opts), add = TRUE)
  options(warn = 1)
  
  on.exit(
    if (exists(".savecall", envir = rlang::current_env())) {
      if (.savecall) {
        mget(names(formals()), envir = rlang::current_env(), inherits = FALSE) %>% 
          c(dots) %>% 
          save_call("standPrn")
      }
    }, 
    add = TRUE
  )
  
  
  reload_expts()
  
  method_align <- rlang::enexpr(method_align)
  if (method_align == rlang::expr(c("MC", "MGKernel"))) {
    method_align <- "MC"
  } else {
    method_align <- rlang::as_string(method_align)
    stopifnot(method_align %in% c("MC", "MGKernel"), 
              length(method_align) == 1)
  }
  
  range_log2r <- prep_range(range_log2r)
  range_int <- prep_range(range_int)

  id <- match_call_arg(normPSM, group_pep_by)
  pep_id <- match_call_arg(normPSM, group_psm_by)
  
  col_select <- rlang::enexpr(col_select)
  col_select <- ifelse(is.null(col_select), rlang::expr(Sample_ID), rlang::sym(col_select))
  load(file = file.path(dat_dir, "label_scheme.rda"))
  
  ok_existing_params(file.path(dat_dir, "Protein/Histogram/MGKernel_params_N.txt"))
  
  if (is.null(label_scheme[[col_select]])) {
    col_select <- rlang::expr(Sample_ID)
    warning("Column \'", rlang::as_string(col_select), "\' does not exist.
			Use column \'Sample_ID\' instead.", call. = FALSE)
  } else if (sum(!is.na(label_scheme[[col_select]])) == 0) {
    col_select <- rlang::expr(Sample_ID)
    warning("No samples were specified under column \'", rlang::as_string(col_select), "\'.
			Use column \'Sample_ID\' instead.", call. = FALSE)
  }
  
  dots <- rlang::enexprs(...)
  
  filename <- file.path(dat_dir, "Protein/Protein.txt")
  
  if (!file.exists(filename)) {
    stop(filename, " not found; run `Pep2Prn` first.")
  }
  
  df <- read.csv(filename, sep = "\t", check.names = FALSE, 
                 header = TRUE, comment.char = "#") 

  local({
    if (sum(grepl("^log2_R[0-9]+ ", names(df))) <= 1) {
      stop("Need more than one sample for `standPep` or `standPrn`.", 
           call. = FALSE)
    }
  })
  
  message("Primary column keys in `Protein/Protein.txt` for `slice_` varargs.")
  
  df <- normMulGau(
    df = df, 
    method_align = method_align, 
    n_comp = n_comp, 
    seed = seed, 
    range_log2r = range_log2r, 
    range_int = range_int, 
    filepath = file.path(dat_dir, "Protein/Histogram"), 
    col_select = col_select, 
    cut_points = cut_points, 
    !!!dots,
  )
  
  df <- df %>% 
    dplyr::filter(!nchar(as.character(.[["prot_acc"]])) == 0) %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
    dplyr::mutate_at(vars(grep("I[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 0)) %>% 
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), as.numeric) %>% 
    dplyr::mutate_at(vars(grep("log2_R[0-9]{3}[NC]*", names(.))), ~ round(.x, digits = 3)) %T>% 
    write.table(file.path(dat_dir, "Protein", "Protein.txt"), 
                sep = "\t", col.names = TRUE, row.names = FALSE)
  
  .saveCall <- TRUE
  
  invisible(df)
}

