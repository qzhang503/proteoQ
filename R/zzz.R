
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!\n\n",
                        "======================================================================\n",
                        "NEW features:\n",
                        "See ?calc_pepmasses and ?calc_monopeptide for peptide masses.\n",
                        # "See ?matchMS for ion matches.\n", 
                        "Notes:\n",
                        "remotes::install_version(\"RSQLite\", version = \"2.2.5\")\n", 
                        "======================================================================\n")
}
