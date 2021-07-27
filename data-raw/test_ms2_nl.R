foo_test_ms2_nls <- function () {
   # (1) "amods- tmod- vnl- fnl-"
   fixedmods = NULL
   varmods = NULL
   mod_indexes <- seq_along(c(fixedmods, varmods)) %>% 
      as.hexmode() %>% 
      `names<-`(c(fixedmods, varmods))
   aa_masses_all <- calc_aamasses(fixedmods, varmods)
   out <- calc_monopep("HQGVMNVGMGQKMN", aa_masses_all[[1]])
   stopifnot(out == 1529.69012)
   out2 <- calc_ms2ions(names(out), ms1_mass = NULL, aa_masses_all[[1]], mod_indexes)
   out3 <- calc_ms2ions(names(out), ms1_mass = out, aa_masses_all[[1]], mod_indexes)


   # (2) "amods- tmod+ vnl- fnl-"
   fixedmods = c("Deamidated (N)")
   varmods = c("Acetyl (N-term)")
   mod_indexes <- seq_along(c(fixedmods, varmods)) %>% as.hexmode() %>% `names<-`(c(fixedmods, varmods))
   aa_masses_all <- calc_aamasses(fixedmods, varmods)
   out <- calc_monopep("HQGVMNVGMGQKMNS", aa_masses_all[[2]])
   # HQGVMNVGMGQKMNS 
   # 1660.70075 
   out2 <- calc_ms2ions(names(out), ms1_mass = NULL, aa_masses_all[[2]], mod_indexes)
   out3 <- calc_ms2ions(names(out), ms1_mass = out, aa_masses_all[[2]], mod_indexes)

   # (3) "amods- tmod+ vnl+ fnl-"
   # not existed, if amods- -> vnl+ is NULL
   
   # (4) "amods- tmod- vnl+ fnl-"
   # not existed, if amods- -> vnl+ is NULL
   
   # (5) "amods- tmod+ vnl- fnl+"
   fixedmods = c("Oxidation (M)", "dHex (S)")
   varmods = c("Acetyl (Protein N-term)")
   mod_indexes <- seq_along(c(fixedmods, varmods)) %>% as.hexmode() %>% `names<-`(c(fixedmods, varmods))
   aa_masses_all <- calc_aamasses(fixedmods, varmods)
   out <- calc_monopep("MHQGVMNVGMGQKMNS", aa_masses_all[[2]])
   # out <- map(aa_masses_all, ~ calc_monopep("MHQGVMNVGMGQKMNS", .x))
   out2 <- calc_ms2ions(names(out[1]), out, aa_masses_all[[2]], mod_indexes)
   duplicated(out2)
   
   # bare = 1747.76264
   # [[1]] 16*4+(233-87);  -64*4; -146; -64*4 - 146
   # MHQGVMNVGMGQKMNS MHQGVMNVGMGQKMNS MHQGVMNVGMGQKMNS MHQGVMNVGMGQKMNS 
   # 1957.80021       1701.80707       1811.74230       1555.74916 
   # [[2]]
   # MHQGVMNVGMGQKMNS MHQGVMNVGMGQKMNS MHQGVMNVGMGQKMNS MHQGVMNVGMGQKMNS 
   # 1999.81077       1743.81763       1853.75286       1597.75972 
   
   calc_ms2ionseries(aa_seq = "MHQGVMNVGMGQKMNS", mass = 1999.81077, 
                     aa_masses = aa_masses_all[[2]], mod_indexes = NULL)
      
      
   # (6) "amods- tmod- vnl- fnl+"
   fixedmods = c("TMT6plex (N-term)", "Oxidation (M)", "dHex (S)")
   varmods = NULL
   aa_masses_all <- calc_aamasses(fixedmods, varmods)
   out <- calc_monopep("HQGVMNVGMGQKMNS", aa_masses_all[[1]])
   # HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS 
   # 2039.92774      1847.93288 (-64*3) 1893.86983 (-146) 1701.87497 (-64*3-146)
   out2 <- calc_ms2ions(names(out[1]), out[1], aa_masses_all[[1]], mod_indexes)
   
   # (7) "amods+ tmod- vnl- fnl-"
   fixedmods = c("TMT6plex (N-term)")
   varmods = c("Deamidated (N)", "Carbamidomethyl (C)")
   mod_indexes <- seq_along(c(fixedmods, varmods)) %>% as.hexmode() %>% `names<-`(c(fixedmods, varmods))
   aa_masses_all <- calc_aamasses(fixedmods, varmods)
   out <- calc_monopep("HQGVMCNVGMGQKMNSC", aa_masses_all[[4]])
   out2 <- calc_ms2ions("HQGVMCNVGMGQKMNSC", out[1], aa_masses_all[[4]], mod_indexes)
   duplicated(out2)
   # HQGVMCNVGMGQKMNSC HQGVMCNVGMGQKMNSC HQGVMCNVGMGQKMNSC HQGVMCNVGMGQKMNSC 
   # 2109.909          2166.930          2110.893          2167.914 
   
   aa_masses_all <- calc_aamasses(fixedmods = c("TMT6plex (N-term)"), varmods = c("Deamidated (N)", "Carbamidomethyl (S)"))
   out <- calc_monopep("HQGVMNVGMGQKSMNS", aa_masses_all[[4]])
   # HQGVMNVGMGQKSMNS HQGVMNVGMGQKSMNS HQGVMNVGMGQKSMNS HQGVMNVGMGQKSMNS 
   # 1990.92259       1991.90661       2047.94406       2048.92807 
   
   # (8) "amods+ tmod+ vnl- fnl-"
   fixedmods = c("TMT6plex (K)")
   varmods = c("Deamidated (N)", "Carbamidomethyl (S)", "Acetyl (Protein N-term)")
   mod_indexes <- seq_along(c(fixedmods, varmods)) %>% as.hexmode() %>% `names<-`(c(fixedmods, varmods))
   aa_masses_all <- calc_aamasses(fixedmods, varmods)
   # out <- map(aa_masses_all, ~ calc_monopep("HQGVMNVGMGQKSMNS", .x))
   out <- calc_monopep("HQGVMNVGMGQKSMNS", aa_masses_all[[8]])
   out2 <- calc_ms2ions(names(out[1]), out[1], aa_masses_all[[8]], mod_indexes)
   duplicated(out2)
   # HQGVMNVGMGQKSMNS HQGVMNVGMGQKSMNS HQGVMNVGMGQKSMNS HQGVMNVGMGQKSMNS 
   # 2032.933         2033.917         2089.955         2090.939 
   
   
   # (9) "amods+ tmod- vnl+ fnl-"
   fixedmods = c("TMT6plex (N-term)")
   # varmods = c("dHex (S)", "Deamidated (N)")
   varmods = c("dHex (S)", "Oxidation (M)")
   aa_masses_all <- calc_aamasses(fixedmods, varmods)
   mod_indexes <- seq_along(c(fixedmods, varmods)) %>% as.hexmode() %>% `names<-`(c(fixedmods, varmods))
   out <- calc_monopep("HQGVMNVGMGQKMNS", aa_masses_all[[4]])
   # HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS 
   # 1992.927        1846.869        1993.911        1847.853 
   out2 <- calc_ms2ions(names(out[1]), out[1], aa_masses_all[[4]], mod_indexes)
   duplicated(out2)
   
   # (10) "amods+ tmod+ vnl+ fnl-"
   fixedmods = c("TMT6plex (K)")
   varmods = c("dHex (S)", "Oxidation (M)", "Deamidated (N)", "Acetyl (Protein N-term)")
   mod_indexes <- seq_along(c(fixedmods, varmods)) %>% as.hexmode() %>% `names<-`(c(fixedmods, varmods))
   aa_masses_all <- calc_aamasses(fixedmods, varmods)
   out <- calc_monopep("HQGVMNVGMGQKMNS", aa_masses_all[[16]])
   out2 <- calc_ms2ions("HQGVMNVGMGQKMNS", 1987.91822, aa_masses_all[[16]], mod_indexes)
   out2 <- calc_ms2ions("HQGVMNVGMGQKMNS", 2050.93249, aa_masses_all[[16]], mod_indexes)
   duplicated(out2)
   
   # HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS 
   # 2050.93249      1986.93420      1904.87458      1840.87630      2066.92740      2002.92912 
   # HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS 
   # 1938.93083      1920.86950      1856.87121      1792.87293      2082.92232      2018.92403 
   # HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS 
   # 1954.92575      1890.92746      1936.86441      1872.86613      1808.86784      1744.86956 
   # HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS 
   # 2051.91651      1987.91822      1905.85860      1841.86031      2067.91142      2003.91314 
   # HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS 
   # 1939.91485      1921.85351      1857.85523      1793.85694      2083.90634      2019.90805 
   # HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS HQGVMNVGMGQKMNS 
   # 1955.90977      1891.91148      1937.84843      1873.85014      1809.85186      1745.85357  
   
   
   
   # (11) "amods+ tmod- vnl- fnl+"
   fixedmods = c("TMT6plex (N-term)", "dHex (S)", "Oxidation (M)")
   varmods = c("Deamidated (N)")
   aa_masses_all <- calc_aamasses(fixedmods, varmods)
   mod_indexes <- seq_along(c(fixedmods, varmods)) %>% as.hexmode() %>% `names<-`(c(fixedmods, varmods))
   out <- calc_monopep("HQGVMNVGMGQKMNS", aa_masses_all[[2]])
   out2 <- calc_ms2ions(names(out[1]), out[1], aa_masses_all[[2]], mod_indexes)
   duplicated(out2)
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
   fixedmods = c("TMT6plex (K)", "Oxidation (M)", "dHex (S)")
   varmods = c("Deamidated (N)", "Acetyl (Protein N-term)")
   aa_masses_all <- calc_aamasses(fixedmods, varmods)
   mod_indexes <- seq_along(c(fixedmods, varmods)) %>% as.hexmode() %>% `names<-`(c(fixedmods, varmods))
   out <- calc_monopep("HQGVMNVGMGQKMNSC", aa_masses_all[[4]])
   out2 <- calc_ms2ions(names(out[1]), out[1], aa_masses_all[[4]], mod_indexes)
   duplicated(out2)
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
   fixedmods = c("TMT6plex (N-term)", "Oxidation (M)", "Deamidated (N)")
   varmods = c("dHex (S)")
   aa_masses_all <- calc_aamasses(fixedmods, varmods)
   mod_indexes <- seq_along(c(fixedmods, varmods)) %>% as.hexmode() %>% `names<-`(c(fixedmods, varmods))
   out <- calc_monopep("HQGVMNVGMGQKMNSC", aa_masses_all[[2]])
   out2 <- calc_ms2ions(names(out[1]), aa_masses_all[[2]], mod_indexes)
}
