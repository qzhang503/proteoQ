#' Helper of \link{groupProts}.
#'
#' @param out The data frame from upstream steps.
#' @param out_path The output path.
#' @export
grp_prots <- function (out, out_path = NULL) {
  
  dir.create(file.path(out_path), recursive = TRUE, showWarnings = FALSE)
  
  out <- out %>% dplyr::arrange(pep_seq)
  
  # essential entries
  rows <- (out$pep_issig & (!out$pep_isdecoy) & (out$pep_rank <= 3L))
  
  df <- out[rows, ]
  
  if (nrow(df) > 1L) {
    df <- df %>% proteoQ::groupProts2(out_path)
  } else {
    df <- df %>%
      dplyr::mutate(prot_isess = TRUE,
                    prot_hit_num = 1L,
                    prot_family_member = 1L)
  }
  
  # non-essential entries
  prot_accs <- df %>%
    dplyr::filter(!duplicated(prot_acc), prot_isess) %>%
    `[[`("prot_acc")
  
  df2 <- out[!rows, ] %>%
    dplyr::filter(!duplicated(prot_acc)) %>%
    dplyr::select(prot_acc) %>%
    dplyr::mutate(prot_isess = ifelse(prot_acc %in% prot_accs, TRUE, FALSE)) %>%
    dplyr::left_join(df[, c("prot_acc", "prot_hit_num", "prot_family_member")] %>%
                       dplyr::filter(!duplicated(prot_acc)),
                     by = "prot_acc") %>%
    dplyr::right_join(out[!rows, ], by = "prot_acc") %>%
    dplyr::mutate(pep_literal_unique = NA, pep_razor_unique = NA) %>%
    dplyr::select(names(df))
  
  rm(list = c("prot_accs"))
  
  out <- dplyr::bind_rows(df, df2) %>%
    dplyr::select(-which(names(.) %in% c("prot_n_psm", "prot_n_pep")))
}


#' Groups proteins by shared peptides.
#'
#' Adds columns \code{prot_hit_num} and \code{prot_family_member} to
#' \code{psm.txt}.
#'
#' @param df Interim results from \link{matchMS}.
#' @param out_path The output path.
#' @export
groupProts2 <- function (df, out_path = NULL) {
  
  # `pep_seq` in `df` are all from target and significant;
  # yet target `pep_seq` can be assigned to both target and decoy proteins
  #
  #    prot_acc     pep_seq
  #  1 -GOG8C_HUMAN EEQERLR
  #  2 -GOG8D_HUMAN EEQERLR
  # 11 MNT_HUMAN    EEQERLR
  
  # --- (1) protein ~ peptide map ---
  mat <- proteoQ::map_pepprot2(df[, c("prot_acc", "pep_seq")], out_path)
  gc()
  
  # --- (2) protein ~ protein groups by distance map ---
  grps <- proteoQ::cut_protgrps2(mat, out_path)
  gc()
  
  # --- (3) set covers by groups ---
  sets <- proteoQ::greedysetcover3(mat)
  gc()
  
  if (!is.null(out_path)) {
    saveRDS(sets, file.path(out_path, "prot_pep_setcover.rds")) 
  }
  
  sets <- sets %>% 
    `[[`("prot_acc") %>%
    unique()
  gc()
  
  # --- set aside df0 ---
  df <- df %>% dplyr::mutate(prot_isess = prot_acc %in% sets)
  df0 <- df %>% dplyr::filter(!prot_isess)
  df <- df %>% dplyr::filter(prot_isess)

  mat_ess <- mat[, colnames(mat) %in% unique(df$prot_acc)]
  
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
  
  df <- df %>%
    dplyr::left_join(grps, by = "prot_acc") %>%
    dplyr::bind_rows(df0) %>%
    dplyr::left_join(peps_uniq, by = "pep_seq")
  
  gc()
  
  invisible(df)
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
#' @export
map_pepprot2 <- function (df, out_path = NULL) {

  df <- df[, c("prot_acc", "pep_seq")] 
  df <- unique(df)
  
  gc()
  
  peps <- df$pep_seq
  
  mat <- Matrix::sparse.model.matrix(~ -1 + prot_acc, df)
  colnames(mat) <- gsub("prot_acc", "", colnames(mat))
  mat <- mat == 1L
  rownames(mat) <- peps
  gc()
  
  # ---
  dpeps <- peps[duplicated(peps)]
  drows <- (rownames(mat) %in% dpeps)
  mat0 <- mat[!drows, ]
  mat <- mat[drows, ]
  
  # ---
  mpeps <- unique(rownames(mat))
  
  len <- as.numeric(length(mpeps))
  ncol <- as.numeric(ncol(mat))
  len2 <- len * ncol
  out <- rep(0L, len2)
  
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
    
    if (i %% 100 == 0) gc()
  }
  
  rm(list = c("mat", "mati"))
  gc()
  
  # ---
  if (object.size(out)/1024^3 > 5) {
    size <- 10000
    n_chunks <- ceiling(len/size)
    
    x0 <- NULL
    
    for (i in 1:n_chunks) {
      x <- out[(ncol*(size*(i-1))+1):min(len2, (ncol*(size*i)))]
      x <- Matrix::Matrix(x, ncol = ncol, byrow = TRUE, sparse = TRUE)
      gc()
      x0 <- rbind2(x0, x)
    }
    
    rownames(x0) <- mpeps
    out <- x0
    
    rm(list = c("x", "x0"))
    gc()
  } else {
    out <- Matrix::Matrix(out, ncol = ncol, byrow = TRUE, sparse = TRUE)
    rownames(out) <- mpeps
    gc()
  }

  out <- out == 1L
  gc()
  
  out <- rbind2(out, mat0)
  gc()
  
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
    Matrix::writeMM(out, file = file.path(out_path, "prot_pep_map.mtx"))
    saveRDS(colnames(out), file.path(out_path, "prot_pep_map_col.rds"))
    saveRDS(rownames(out), file.path(out_path, "prot_pep_map_row.rds"))
  }
  
  invisible(out)
}


#' Cuts proteins into groups.
#'
#' By the number of shared peptides.
#'
#' @param mat A logical matrix; peptides in rows and proteins in columns.
#' @param out_path A file pth to outputs.
#' @export
cut_protgrps2 <- function (mat = NULL, out_path = NULL) {

  dista = proxyC::simil(mat, margin = 2) # sparse distance matrix
  cns <- colnames(mat)
  rm(list = c("mat"))
  gc()
  
  # ---
  ncol <- ncol(dista)
  cols <- 1:ncol

  if (ncol > 10000) {
    mat2 <- matrix(nrow = ncol, ncol = ncol)
    colnames(mat2) <- cns
    rownames(mat2) <- cns
    
    max_rc <- 2500 * 40000 # 5000 out of memory
    gc()
    
    size <- ceiling(max_rc/ncol)
    n_chunks <- ceiling(ncol/size)
    
    for (i in 1:n_chunks) {
      rows <- (1+(i-1)*size):min(i*size, ncol)

      x <- dista[rows, cols] 
      gc()
      
      x <- (x == 0) # sparse logical matrix
      gc()
      
      mat2[rows, cols] <- as.matrix(x) # regular matrix
      gc()
    }
    
    rm(list = c("x"))
    gc()
  } else {
    dista <- (dista == 0) # sparse logical matrix
    gc()
    
    dista <- as.matrix(dista) # regular matrix
    gc()
    
    mat2 <- dista
    
    rm(list = c("dista"))
    gc()
  }
  
  # Diagonal values are `FALSE`
  # TRUE - orthogonal (without shared peptides)
  # FALSE - with shared peptides
  # 
  #             KKA1_ECOLX NP_000005 NP_000007
  # KKA1_ECOLX      FALSE      TRUE      TRUE
  # NP_000005        TRUE     FALSE      TRUE
  # NP_000007        TRUE      TRUE     FALSE
  
  # --- finds protein groups
  mat2 <- proteoQ::as_lgldist(mat2, diag = FALSE, upper = FALSE)
  gc()
  
  hc <- hclust(mat2, method = "single")
  gc()
  
  grps <- data.frame(prot_hit_num = cutree(hc, h = .9)) %>% 
    tibble::rownames_to_column("prot_acc") %>%
    dplyr::group_by(prot_hit_num) %>%
    dplyr::mutate(prot_family_member = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  if (!is.null(out_path)) {
    saveRDS(grps, file.path(out_path, "prot_grps.rds"))
  }
  
  stopifnot(identical(grps$prot_acc, cns))
  
  invisible(grps)
}


#' Simplified \link[stats]{as.dist} for memory efficiency.
#' 
#' Assumed the input is already a symmetric matrix.
#' 
#' @inheritParams stats::as.dist
as_dist <- function (m, diag = FALSE, upper = FALSE) {
  
  p <- nrow(m)
  
  ans <- m[row(m) > col(m)]
  gc()
  attributes(ans) <- NULL
  
  if (!is.null(rownames(m))) {
    attr(ans, "Labels") <- rownames(m)
  } else if (!is.null(colnames(m))) {
    attr(ans, "Labels") <- colnames(m)
  }
  
  attr(ans, "Size") <- p
  attr(ans, "call") <- match.call()
  class(ans) <- "dist"
  
  if (is.null(attr(ans, "Diag")) || !missing(diag)) {
    attr(ans, "Diag") <- diag
  }
    
  if (is.null(attr(ans, "Upper")) || !missing(upper)) {
    attr(ans, "Upper") <- upper
  }

  ans
}


#' Simplified \link[stats]{as.dist} for memory efficiency.
#' 
#' Assumed the input is already a symmetric matrix.
#' 
#' @inheritParams stats::as.dist
#' @export
as_lgldist <- function(m, diag = FALSE, upper = FALSE) {
  
  d = proteoCpp::to_lgldistC(m)
  
  if (!is.null(rownames(m))) {
    attr(d, "Labels") <- rownames(m)
  } else if (!is.null(colnames(m))) {
    attr(d, "Labels") <- colnames(m)
  }

  attr(d, "class") = "dist"
  attr(d, "Size") = nrow(m)
  attr(d, "call") = match.call()
  attr(d, "Diag") = diag
  attr(d, "Upper") = upper

  d
}


#' Greedy set cover.
#' 
#' A bool matrix input. Output both essential sets and elements.
#' 
#' @param mat A bool matrix of protein (cols)-peptide (rows) map. 
#' 
#' @return A two-column data frame of prot_acc and pep_seq. 
#' @export
greedysetcover3 <- function (mat) {
  
  if (is.matrix(mat)) {
    mat <- Matrix::Matrix(as.matrix(mat), sparse = TRUE)
    gc()
  }
  
  prot_acc <- NULL
  pep_seq <- NULL
  
  while(nrow(mat)) {
    max <- which.max(Matrix::colSums(mat, na.rm = TRUE))
    
    if (max == 0L) break
    
    prot <- names(max)
    rows <- which(mat[, max])
    # peps <- names(rows) # name dropped if only one row
    peps <- rownames(mat)[rows]
    
    prot_acc <- c(prot_acc, rep(prot, length(peps)))
    pep_seq <- c(pep_seq, peps)
    
    mat <- mat[-rows, -max, drop = FALSE]
  }
  
  rm(list = c("mat"))
  gc()
  
  dplyr::bind_cols(prot_acc = prot_acc, pep_seq = pep_seq)
}


