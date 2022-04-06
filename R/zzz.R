
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!\n\n",
                        "=============================================================================\n",
                        "NEW features (v1.7.3.0):\n",
                        # "Supports SILAC-like analysis (One sample with multiple MS1 groups).\n", 
                        "Database search engine migrated to: ", 
                        "https://github.com/qzhang503/proteoM\n",
                        "=============================================================================\n")
}
