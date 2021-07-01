
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!\n\n",
                        "======================================================================\n",
                        "NEW features:\n",
                        "?matchMS for database searches (requires MGFs as inputs).\n", 
                        # "?calc_monopeptide for peptide masses.\n",
                        # "?calc_ms2ionseries for MS2 ions.\n", 
                        # "Notes:\n",
                        # "remotes::install_version(\"RSQLite\", version = \"2.2.5\")\n", 
                        "======================================================================\n")
}
