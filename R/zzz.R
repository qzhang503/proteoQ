
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!\n\n",
                        "======================================================================\n",
                        "NEW features:\n",
                        "?calc_monopeptide for peptide masses.\n",
                        "?calc_ms2ionseries for MS2 ions.\n", 
                        "?matchMS (a light-weight search engine).\n", 
                        # "Notes:\n",
                        # "remotes::install_version(\"RSQLite\", version = \"2.2.5\")\n", 
                        "======================================================================\n")
}
