
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!\n\n",
                        "======================================================================\n",
                        "NEW features:\n",
                        "See ?pepLDA and ?prnLDA for reduced-rank linear discriminant analysis.\n",
                        "======================================================================\n")
}
