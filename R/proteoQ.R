#' proteoQ: A package for processing mass spectrometry data.
#'
#' The proteoQ package provides five categories of functions in (1) database
#' searches, (2) data preprocessing, (3) annotation, (4) quality assurance
#' assessments and (5) informatic analysis.
#'
#' @section Main utility in database searches: matchMS (additional tool kits:
#'   load_fasta2, parse_unimod, calc_monopeptide, calc_ms2ionseries)
#' @section Utilities in data preprocessing: normPSM -> (purgePSM) -> PSM2Pep ->
#'   mergePep -> standPep -> (purgePep) -> standPrn
#' @section Utilities in custom database preparation: prepGO, prepMSig,
#'   prepString, Uni2Entrez, Ref2Entrez
#' @section Utilities in data QA and QC: pepHist, prnHist, pepCorr_logFC,
#'   pepCorr_logInt, prnCorr_logFC, prnCorr_logInt,
#' @section Utilities in informatic analysis: pepMDS, prnMDS, pepPCA, prnPCA,
#'   pepLDA, prnLDA, pepEucDist, prnEucDist, pepHM, prnHM, anal_pepNMF,
#'   anal_prnNMF, plot_pepNMFCon, plot_prnNMFCon, plot_pepNMFCoef,
#'   plot_prnNMFCoef, plot_metaNMF, anal_prnTrend, plot_prnTrend, pepSig,
#'   prnSig, pepVol, prnVol, prnGSPA, prnGSPAHM, gspaMap, anal_prnString,
#'   pepImp, prnImp, prnGSVA, prnGSEA, cluego
#'
#' @docType package
#' @name proteoQ
NULL
#> NULL
