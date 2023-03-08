
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!\n\n",
                        "=============================================================================\n",
                        "NEW features (v1.7.5.0):\n",
                        "[x] Updated significant tests for LFQ\n", 
                        "[x] Added match between runs for LFQ\n", 
                        "[x] Database searches: ", "https://github.com/qzhang503/proteoM\n",
                        "=============================================================================\n")
}
