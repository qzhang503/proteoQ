
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!\n\n",
                        "=============================================================================\n",
                        "NEW features (v1.7.4.0):\n",
                        "[x] \"scale_log2r = NA\" to work against data without normalization\n", 
                        "[x] Database search engine migrated to: ", 
                        "https://github.com/qzhang503/proteoM\n",
                        "=============================================================================\n")
}
