#' Matches the current id to the id in normPep or normPrn
#'
#' @param id.
#'
#' @import plyr dplyr purrr rlang
#' @importFrom magrittr %>%
match_identifier <- function (id = c("pep_seq", "pep_seq_mod", "prot_acc", "gene")) {
  if (id %in% c("prot_acc", "gene")) {
    fn_pars <- "normPrn.txt"
  } else if (id %in% c("pep_seq", "pep_seq_mod")) {
    fn_pars <- "normPep.txt"
  }
  
  call_pars <- tryCatch(read.csv(file.path(dat_dir, "Calls", fn_pars), check.names = FALSE,
                                 header = TRUE, sep = "\t", comment.char = "#"),
                        error = function(e) NA)
  
  if (!is.null(dim(call_pars))) {
    id <- call_pars %>%
      dplyr::filter(var == "id") %>%
      dplyr::select("value.1") %>%
      unlist() %>%
      as.character()
  }
  
  return(id)
}


#' Add Z_log2_R, sd_N_log2_R and sd_Z_log2_R
calc_more_psm_sd <- function (df, group_psm_by, range_log2r, range_int, set_idx, injn_idx) {
  # SD columns for "N_log2_R"
  df <- df %>% 
    calcSD_Splex(group_psm_by, "N_log2_R") %>% 
    `names<-`(gsub("^N_log2_R", "sd_N_log2_R", names(.))) %>% 
    dplyr::right_join(df, by = group_psm_by)
  
  # "Z_log2_R" columns and their SD
  cf_SD <- calc_sd_fcts_psm(df, range_log2r, range_int, set_idx, injn_idx)
  
  cf_x <- df %>%
    dplyr::select(matches("^N_log2_R[0-9]{3}")) %>%
    `colnames<-`(gsub(".*\\s*\\((.*)\\)$", "\\1", names(.))) %>%
    dplyr::summarise_all(~ median(.x, na.rm = TRUE)) 
  
  df <- mapply(normSD, df[,grepl("^N_log2_R[0-9]{3}", names(df))],
               center = cf_x, SD = cf_SD$fct, SIMPLIFY = FALSE) %>%
    data.frame(check.names = FALSE) %>%
    `colnames<-`(gsub("N_log2", "Z_log2", names(.))) %>%
    cbind(df, .)
  
  df_z <- df %>% dplyr::select(grep("^Z_log2_R", names(.)))
  nan_cols <- purrr::map_lgl(df_z, is_all_nan, na.rm = TRUE)
  df_z[, nan_cols] <- 0
  df[, grep("^Z_log2_R", names(df))] <- df_z
  rm(df_z, nan_cols)
  
  df <- df %>% 
    calcSD_Splex(group_psm_by, "Z_log2_R") %>% 
    `names<-`(gsub("^Z_log2_R", "sd_Z_log2_R", names(.))) %>% 
    dplyr::right_join(df, by = group_psm_by)
  
  df <- dplyr::bind_cols(
    df %>% dplyr::select(-grep("[RI]{1}[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("I[0-9]{3}[NC]*", names(.))), 
    df %>% dplyr::select(grep("R[0-9]{3}[NC]*", names(.))), 
  )
}


#' Median-centering normalization
#' 
#' # not currently used?
#' 
#' @import dplyr purrr rlang  magrittr
logfcPep <- function(df, label_scheme, set_idx) {
  label_scheme_sub <- label_scheme %>% 
    .[.$TMT_Set == set_idx & .$LCMS_Injection == 1, ]
  
  channelInfo <- channelInfo(label_scheme_sub, set_idx)
  
  col_sample <- grep("^I[0-9]{3}", names(df))
  
  if (length(channelInfo$refChannels) > 0) {
    ref_index <- channelInfo$refChannels
  } else {
    ref_index <- channelInfo$labeledChannels
  }
  
  df <- sweep(df[, col_sample, drop = FALSE], 1,
              rowMeans(df[, col_sample[ref_index], drop = FALSE], na.rm = TRUE), "/") %>%
    log2(.) %>%
    `colnames<-`(gsub("I", "log2_R", names(.)))	%>%
    cbind(df, .) %>%
    dplyr::mutate_at(.vars = grep("[I|R][0-9]{3}", names(.)), ~ replace(.x, is.infinite(.x), NA))
  
  col_log2Ratio <- grepl("^log2_R[0-9]{3}", names(df))
  cf <- apply(df[, col_log2Ratio, drop = FALSE], 2, median, na.rm = TRUE)
  
  df <- sweep(df[, col_log2Ratio, drop = FALSE], 2, cf, "-") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    cbind(df, .)
  
  df <- sweep(df[, grepl("^I[0-9]{3}", names(df)), drop = FALSE], 2, 2^cf, "/") %>%
    `colnames<-`(paste("N", names(.), sep="_"))	%>%
    cbind(df, .)
  
  df <- df %>%
    na_zeroIntensity()
}


