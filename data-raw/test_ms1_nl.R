foo_test_ms1_nls <- function () {
   # (1) "amods- tmod- vnl- fnl-"
   aa_masses_all <- calc_aamasses(fixedmods = NULL, varmods = NULL)
   out <- calc_monopep("HQGVMNVGMGQKMN", aa_masses_all[[1]])
   stopifnot(out == 1529.69012)
   
   aa_masses_all <- calc_aamasses(fixedmods = c("Deamidated (N)"), varmods = c("Acetyl (N-term)"))
   # (1) "amods- tmod- vnl- fnl-"
   out <- calc_monopep("HQGVMNVGMGQKMNS", aa_masses_all[[1]])
   # (2) "amods- tmod+ vnl- fnl-"
   out <- calc_monopep("HQGVMNVGMGQKMNS", aa_masses_all[[2]])
   
   # (3) "amods- tmod+ vnl+ fnl-"
   # not existed, if amods- -> vnl+ is NULL
   aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (K)"), varmods = c("Fake-Oxidation (N-term)"))
   out <- calc_monopep("MHQGVMNVGMGQKMNS", aa_masses_all[[2]])
   # aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (K)"), varmods = c("Acetyl (N-term)", "Met-loss (Protein N-term = M)"))
   # out <- calc_monopep("MHQGVMNVGMGQKMNS", aa_masses_all[[3]])
   
   # (4) "amods- tmod- vnl+ fnl-"
   # not existed, if amods- -> vnl+ is NULL
   
   # (5) "amods- tmod+ vnl- fnl+"
   aa_masses_all <- calc_aamasses(fixedmods = c("Oxidation (M)", "dHex (S)"), varmods = "Acetyl (Protein N-term)")
   out <- map(aa_masses_all, ~ calc_monopep("MHQGVMNVGMGQKMNS", .x, include_insource_nl = TRUE))
   # bare = 1747.76264
   # [[1]] 16*4+(233-87);  -64*4; -146; -64*4 - 146
   # MHQGVMNVGMGQKMNS MHQGVMNVGMGQKMNS MHQGVMNVGMGQKMNS MHQGVMNVGMGQKMNS 
   # 1957.80021       1701.80707       1811.74230       1555.74916 
   # [[2]]
   # MHQGVMNVGMGQKMNS MHQGVMNVGMGQKMNS MHQGVMNVGMGQKMNS MHQGVMNVGMGQKMNS 
   # 1999.81077       1743.81763       1853.75286       1597.75972 
   
   # (6) "amods- tmod- vnl- fnl+"
   aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "Oxidation (M)", "dHex (S)"), varmods = NULL)
   out <- calc_monopep("HQGVMNVGMGQKMNS", aa_masses_all[[1]], include_insource_nl = TRUE)
   # HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS 
   # 2039.92774      1847.93288 (-64*3) 1893.86983 (-146) 1701.87497 (-64*3-146)
   
   # vmods_combi NULL
   # calc_monopep("MHQGVNVGGQKNS", aa_masses_all[[1]])
   
   # (7) "amods+ tmod- vnl- fnl-"
   aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (N-term)"), varmods = c("Deamidated (N)", "Carbamidomethyl (C)"))
   out <- calc_monopep("HQGVMCNVGMGQKMNSC", aa_masses_all[[4]])
   # HQGVMCNVGMGQKMNSC HQGVMCNVGMGQKMNSC HQGVMCNVGMGQKMNSC HQGVMCNVGMGQKMNSC 
   # 2109.909          2166.930          2110.893          2167.914 
   
   aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (N-term)"), varmods = c("Deamidated (N)", "Carbamidomethyl (S)"))
   out <- calc_monopep("HQGVMNVGMGQKSMNS", aa_masses_all[[4]])
   # HQGVMNVGMGQKSMNS HQGVMNVGMGQKSMNS HQGVMNVGMGQKSMNS HQGVMNVGMGQKSMNS 
   # 1990.92259       1991.90661       2047.94406       2048.92807 
   
   # (8) "amods+ tmod+ vnl- fnl-"
   aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (K)"), 
         varmods = c("Deamidated (N)", "Carbamidomethyl (S)", "Acetyl (Protein N-term)"))
   out <- map(aa_masses_all, ~ calc_monopep("HQGVMNVGMGQKSMNS", .x))
   out <- calc_monopep("HQGVMNVGMGQKSMNS", aa_masses_all[[8]])
   # HQGVMNVGMGQKSMNS HQGVMNVGMGQKSMNS HQGVMNVGMGQKSMNS HQGVMNVGMGQKSMNS 
   # 2032.933         2033.917         2089.955         2090.939 
   
   
   # (9) "amods+ tmod- vnl+ fnl-"
   aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (N-term)"), 
                                  varmods = c("dHex (S)", "Deamidated (N)"))
   out <- calc_monopep("HQGVMNVGMGQKMNS", aa_masses_all[[4]], include_insource_nl = TRUE)
   # HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS 
   # 1992.927        1846.869        1993.911        1847.853 
   
   
   # (10) "amods+ tmod+ vnl+ fnl-"
   aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (K)"), 
         varmods = c("dHex (S)", "Oxidation (M)", "Deamidated (N)", "Acetyl (Protein N-term)"))
   out <- calc_monopep("HQGVMNVGMGQKMNS", aa_masses_all[[16]])
   # HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS 
   # 2050.932        1986.934        1904.875        1840.876        2066.927        1938.931        1920.870 
   # HQGVMNVGMGQKMNS 
   # 1792.873 
   
   
   # (11) "amods+ tmod- vnl- fnl+"
   aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "dHex (S)", "Oxidation (M)"), varmods = c("Deamidated (N)"))
   out <- calc_monopep("HQGVMNVGMGQKMNS", aa_masses_all[[2]])
   # HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS 
   # 2040.912        1894.854        1848.917        1702.859        2041.896        1895.838        1849.901 
   # HQGVMNVGMGQKMNS 
   # 1703.843 
   
   aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (K)", "dHex (S)", "Oxidation (M)"), 
         varmods = c("Deamidated (N)", "Acetyl (Protein N-term)"))
   out <- calc_monopep("HQGVMNVGMGQKMNS", aa_masses_all[[4]])
   # HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS 
   # 2082.922        1936.864        1890.927        1744.870        2083.906        1937.848        1891.911 
   # HQGVMNVGMGQKMNS 
   # 1745.854 
   
   
   # (12) "amods+ tmod+ vnl- fnl+"
   aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (K)", "Oxidation (M)", "dHex (S)"), 
         varmods = c("Deamidated (N)", "Acetyl (Protein N-term)"))
   out <- calc_monopep("HQGVMNVGMGQKMNSC", aa_masses_all[[4]])
   # HQGVMNVGMGQKMNSC HQGVMNVGMGQKMNSC HQGVMNVGMGQKMNSC HQGVMNVGMGQKMNSC HQGVMNVGMGQKMNSC HQGVMNVGMGQKMNSC 
   # 2185.932         1993.937         2039.874         1847.879         2186.916         1994.921 
   # HQGVMNVGMGQKMNSC HQGVMNVGMGQKMNSC 
   # 2040.858         1848.863 
   
   aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (K)", "Oxidation (M)", "dHex (S)"), 
         varmods = c("Deamidated (N)", "Carbamidomethyl (C)", "Acetyl (Protein N-term)"))
   out <- calc_monopep("HQGVMNVGMGQKMNSC", aa_masses_all[[8]])
   # HQGVMNVGMGQKMNSC HQGVMNVGMGQKMNSC HQGVMNVGMGQKMNSC HQGVMNVGMGQKMNSC HQGVMNVGMGQKMNSC HQGVMNVGMGQKMNSC 
   # 2242.953         2050.958         2096.895         1904.900         2243.937         2051.942 
   # HQGVMNVGMGQKMNSC HQGVMNVGMGQKMNSC 
   # 2097.879         1905.884 
   
   
   # (13-16) ignored
   aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (N-term)", "Oxidation (M)", "Deamidated (N)"), 
         varmods = c("dHex (S)"))
   out <- calc_monopep("HQGVMNVGMGQKMNS", aa_masses_all[[2]])
}
