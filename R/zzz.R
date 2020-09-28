
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!\n\n",
                        "======================================================================\n",
                        "New features:\n",
                        "* MSFragger LFQ.\n",
                        "* Multiple and/or custom FASTA files.\n", 
                        "* ?pepLDA and ?prnLDA for reduced-rank linear discriminant analysis.\n",
                        "\n", 
                        "Notes:\n", 
                        "* New preferences in Mascot PSM exorts\n", 
                        "  (https://proteoq.netlify.app/post/exporting-mascot-psms/).\n",
                        "* ?PSM2Pep and ?Pep2Prn for new default in LFQ normalization\n", 
                        "======================================================================\n")
}
