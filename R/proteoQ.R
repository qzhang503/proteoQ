#' proteoQ: A package for processing mass spectrometry data using tandem mass
#' tags (\url{https://en.wikipedia.org/wiki/Tandem_mass_tag}).
#'
#' The proteoQ package provides four categories of functions in data
#' preprocessing, annotation, quality assurance assessments and informatic
#' analysis.
#'
#' @section functions in data preprocessing: normPSM -> (purgePSM) ->
#'   PSM2Pep -> mergePep -> standPep -> (purgePep) -> standPrn
#' @section functions in custom database preparation: prepGO, prepMSig,
#'   prepString, Uni2Entrez, Ref2Entrez
#' @section functions in data QA: pepHist, prnHist, pepCorr_logFC,
#'   pepCorr_logInt, prnCorr_logFC, prnCorr_logInt,
#' @section functions in informatic analysis: pepMDS, prnMDS, pepPCA,
#'   prnPCA, pepEucDist, prnEucDist, pepHM, prnHM, anal_pepNMF, anal_prnNMF,
#'   plot_pepNMFCon, plot_prnNMFCon, plot_pepNMFCoef, plot_prnNMFCoef,
#'   plot_metaNMF, anal_prnTrend, plot_prnTrend, pepSig, prnSig, pepVol, prnVol,
#'   prnGSPA, prnGSPAHM, gspaMap, anal_prnString, pepImp, prnImp, prnGSVA,
#'   prnGSEA, cluego
#'
#' @docType package
#' @name proteoQ
NULL
#> NULL
