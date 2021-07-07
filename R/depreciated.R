#' Groups proteins by shared peptides.
#'
#' Adds columns \code{prot_hit_num} and \code{prot_family_member} to
#' \code{psmQ.txt}.
#'
#' @param df Interim results from \link{matchMS}.
#' @param out_path The output path.
groupProts_orig <- function (df, out_path = "~/proteoQ/outs") {
  
  # pep_seq in df are all from target and significant;
  # yet target pep_seq can be assigned to both target and decoy proteins
  #
  #    prot_acc     pep_seq
  #  1 -GOG8C_HUMAN EEQERLR
  #  2 -GOG8D_HUMAN EEQERLR
  # 11 MNT_HUMAN    EEQERLR
  
  # --- protein ~ peptide map ---
  mat <- map_pepprot(df[c("prot_accc", "pep_seq")], out_path) 
  
  # --- set aside df0 ---
  message("Grouping proteins by families.")
  
  sets <- mat %>% 
    greedysetcover3() %T>%
    saveRDS(file.path(out_path, "prot_sets.rds")) %>%
    `[[`("prot_acc") %>%
    unique()
  
  df <- df %>% dplyr::mutate(prot_isess = prot_acc %in% sets)
  df0 <- df %>% filter(!prot_isess)
  df1 <- df %>% filter(prot_isess)
  rm(list = c("df"))
  
  mat_ess <- mat[colnames(mat) %in% df1$prot_acc]
  
  peps_uniq <- local({
    rsums <- rowSums(mat)
    rsums2 <- rowSums(mat_ess)
    
    peps <- data.frame(pep_seq = rownames(mat)) %>%
      dplyr::mutate(pep_literal_unique = (rsums == 1L)) %>%
      dplyr::mutate(pep_razor_unique = (rsums2 == 1L))
  })
  
  # --- protein ~ protein distance map
  grps <- cut_protgrps(mat_ess, out_path)
  
  # --- put together
  df0 <- df0 %>%
    dplyr::mutate(prot_hit_num = NA, prot_family_member = NA)
  
  df1 <- df1 %>%
    dplyr::left_join(grps, by = "prot_acc") %>%
    dplyr::bind_rows(df0) %>%
    dplyr::left_join(peps_uniq, by = "pep_seq")
  
  invisible(df1)
}

