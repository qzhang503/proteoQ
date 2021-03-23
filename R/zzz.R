
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!\n\n",
                        "======================================================================\n",
                        "NEW features in MS1 peptide masses:\n",
                        "See ?calc_pepmasses and ?calc_monopeptide.\n",
                        "======================================================================\n")
}
