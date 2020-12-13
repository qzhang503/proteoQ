
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to proteoQ!\n\n",
                        "==================================================================\n",
                        "Major rewrite in data preprocessing.\n\n", 
                        "Updates:\n",
                        "* Integrated (TMT and LFQ) analysis for various search engines.\n",
                        "* Versatile choices in Mascot PSM exports.\n", 
                        "==================================================================\n")
}
