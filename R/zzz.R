
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!\n\n", 
                        "=====================================================================\n", 
                        "NEW features:\n", 
                        "See ?pepLDA and ?prnLDA in reduced-rank linear discriminant analysis.\n", 
                        "=====================================================================\n")
}
