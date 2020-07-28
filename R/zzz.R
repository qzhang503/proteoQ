
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!\n\n",
                        "======================================================================\n",
                        "NEW features:\n",
                        "See ?pepLDA and ?prnLDA for reduced-rank linear discriminant analysis.\n",
                        "See https://proteoq.netlify.app/post/metadata-files-for-lfq/ \n",
                        "\tfor data mining against MaxQuant LFQ.\n", 
                        "======================================================================\n")
}
