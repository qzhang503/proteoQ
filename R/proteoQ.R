#' proteoQ: A package for processing mass spectrometry data.
#'
#' The proteoQ provides four categories of functions in data preprocessing, (2)
#' annotation, (3) quality assurance assessments and (4) informatic analysis.
#' The companion package, mzion, performs database searches.
#'
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
