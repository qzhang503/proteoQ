#' Groups proteins by shared peptides.
#'
#' Adds columns \code{prot_hit_num} and \code{prot_family_member} to
#' \code{psm.txt}.
#'
#' @param df Interim results from \link{matchMS}.
#' @param out_path The output path.
groupProts2 <- function (df, out_path = NULL) {
  
  # `pep_seq` in `df` are all from target and significant;
  # yet target `pep_seq` can be assigned to both target and decoy proteins
  #
  #    prot_acc     pep_seq
  #  1 -GOG8C_HUMAN EEQERLR
  #  2 -GOG8D_HUMAN EEQERLR
  # 11 MNT_HUMAN    EEQERLR
  
  # --- (1) protein ~ peptide map ---
  mat <- map_pepprot2(df[, c("prot_acc", "pep_seq")], out_path)
  gc()
  
  # --- (2) protein ~ protein groups by distance map ---
  grps <- cut_protgrps2(mat, out_path)
  gc()
  
  # --- (3) set covers by groups ---
  sets <- greedysetcover3(mat)
  gc()
  
  if (!is.null(out_path)) {
    saveRDS(sets, file.path(out_path, "prot_sets.rds")) 
  }
  
  sets <- sets %>% 
    `[[`("prot_acc") %>%
    unique()
  gc()
  
  # --- set aside df0 ---
  df <- df %>% dplyr::mutate(prot_isess = prot_acc %in% sets)
  df0 <- df %>% filter(!prot_isess)
  df1 <- df %>% filter(prot_isess)

  mat_ess <- mat[, colnames(mat) %in% unique(df1$prot_acc)]
  
  peps_uniq <- local({
    rsums <- Matrix::rowSums(mat)
    rsums2 <- Matrix::rowSums(mat_ess)
    
    peps <- data.frame(pep_seq = rownames(mat)) %>%
      dplyr::mutate(pep_literal_unique = (rsums == 1L)) %>%
      dplyr::mutate(pep_razor_unique = (rsums2 == 1L))
  })
  
  # --- put together ---
  df0 <- df0 %>%
    dplyr::mutate(prot_hit_num = NA, prot_family_member = NA)
  
  df1 <- df1 %>%
    dplyr::left_join(grps, by = "prot_acc") %>%
    dplyr::bind_rows(df0) %>%
    dplyr::left_join(peps_uniq, by = "pep_seq")
  
  gc()
  
  invisible(df1)
}


#' Helper of \link{groupProts}.
#'
#' Builds the logical map between peptide (in rows) and proteins (in columns).
#'
#' @param df The data frame from upstream steps. It must contains the two
#'   columns of \code{prot_acc} and \code{pep_seq}.
#' @param out_path An output path.
#' @examples 
#' \donttest{
#' df <- data.frame(prot_acc = character(2000), pep_seq = character(2000))
#' set.seed(100)
#' df$prot_acc <- sample(LETTERS[1:20], 2000, replace = TRUE)
#' df$pep_seq <- sample(letters[1:26], 20, replace = TRUE)
#' df <- df[!duplicated(df), ] 
#' 
#' out <- map_pepprot2(df)
#' }
map_pepprot2 <- function (df, out_path = NULL) {
  message("Building protein-peptide maps.")
  
  df <- df[, c("prot_acc", "pep_seq")] %>%
    tidyr::unite(pep_prot., pep_seq, prot_acc, sep = "@", remove = FALSE) %>%
    dplyr::filter(!duplicated(pep_prot.)) %>%
    dplyr::select(-pep_prot.)
  gc()
  
  peps <- df$pep_seq
  
  mat = Matrix::sparse.model.matrix(~ -1 + prot_acc, df)
  colnames(mat) <- gsub("prot_acc", "", colnames(mat))
  mat <- mat == 1L
  rownames(mat) <- peps
  gc()
  
  # ---
  dpeps <- peps[duplicated(peps)]
  drows <- rownames(mat) %in% dpeps
  mat0 <- mat[!drows, ]
  mat <- mat[drows, ]
  
  # ---
  mpeps <- unique(rownames(mat))
  len <- length(mpeps)
  
  ncol <- ncol(mat)
  out <- rep(0, len * ncol)
  
  start <- 1
  end <- ncol

  for (i in seq_len(len)) {
    pep <- rownames(mat)[[1]]
    rows <- rownames(mat) == pep
    
    mati <- mat[rows, ]
    out[start:end] <- Matrix::colSums(mati)

    mat <- mat[!rows, ]
    
    start <- start + ncol
    end <- end + ncol
  }
  
  rm(list = c("mat", "mati"))
  gc()
  
  out <- Matrix::Matrix(out, ncol = ncol, byrow = TRUE, sparse = TRUE)
  rownames(out) <- mpeps
  gc()
  
  out <- rbind2(out, mat0)
  gc()
  
  # ---

  # collapse rows of the same pep_seq; may use `sum`
  #
  #      pep_seq prot_acc
  # 1       A        X
  # 2       A        Y
  # 3       B        X
  # 4       C        Y
  #
  #   X Y
  # 1 1 0
  # 2 0 1
  # 3 1 0
  # 4 0 1
  #
  # A tibble: 3 x 3
  # pep_seq X     Y
  # <chr>   <lgl> <lgl>
  # 1 A       TRUE  TRUE
  # 2 B       TRUE  FALSE
  # 3 C       FALSE TRUE
  
  if (!is.null(out_path)) {
    Matrix::writeMM(out, file.path(out_path, "prot_pep_map.mtx"))
  }
  
  invisible(out)
}


#' Cuts proteins into groups.
#'
#' By the number of shared peptides.
#'
#' @param mat A logical matrix; peptides in rows and proteins in columns.
#' @param out_path A file pth to outputs.
cut_protgrps2 <- function (mat, out_path = NULL) {

  out = proxyC::simil(mat, margin = 2)
  gc()
  
  out <- as.matrix(out)
  gc()
  
  out <- round(out, digits = 2)
  gc()

  # f <- function(x, y) sum(x & y)
  # out <- proxy::dist(mat, by_rows = FALSE, method = f)

  cns <- colnames(mat)
  
  # --- finds protein groups
  out[out == 0] <- 2
  out[out < 2] <- 0
  gc()
  
  out <- out %>% as.dist(diag = TRUE, upper = TRUE)
  gc()
  
  hc <- hclust(out, method = "single")
  gc()
  
  grps <- data.frame(prot_hit_num = cutree(hc, h = 1)) %>%
    tibble::rownames_to_column("prot_acc") %>%
    dplyr::group_by(prot_hit_num) %>%
    dplyr::mutate(prot_family_member = row_number()) %>%
    dplyr::ungroup() 
  
  if (!is.null(out_path)) {
    saveRDS(grps, file.path(out_path, "prot_grps.rds"))
  }
  
  stopifnot(identical(grps$prot_acc, cns))
  
  invisible(grps)
}

