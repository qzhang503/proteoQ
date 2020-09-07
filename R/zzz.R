
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!\n\n",
                        "======================================================================\n",
                        "NEW features:\n",
                        "See ?pepLDA and ?prnLDA for reduced-rank linear discriminant analysis.\n",
                        "See https://proteoq.netlify.app/post/metadata-files-for-lfq/ \n",
                        "\tfor data mining against MaxQuant LFQ.\n", 
                        "See ?mergePep and ?Pep2Prn for data alignment by sections:\n",
                        "\tmergePep(cut_points = seq(4, 7, .5))\n", 
                        "\tPep2Prn(cut_points = seq(4, 7, .5))\n", 
                        "Support multiple and/or custom FASTA files.\n", 
                        "======================================================================\n")
}
