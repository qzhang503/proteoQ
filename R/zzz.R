
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!\n\n",
                        "=============================================================================\n",
                        "NEW features (v1.6.1.2):\n",
                        "`matchMS` database search takes MGFs from MSConvert or Proteome Discoverer.\n", 
                        "\n", 
                        "(the utility is under active development and may require \n", 
                         "the removals of caches under `.MSearch` and `pepmasses`).\n", 
                        # "?calc_monopeptide for peptide masses.\n",
                        # "Notes:\n",
                        # "remotes::install_version(\"RSQLite\", version = \"2.2.5\")\n", 
                        "=============================================================================\n")
}
