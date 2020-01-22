prnVol_ed <- function (id = "gene", anal_type = "Volcano", df = NULL, 
                       scale_log2r = TRUE, complete_cases = FALSE, impute_na = FALSE, 
                       filepath = NULL, filename = NULL, fml_nms = NULL, 
                       adjP = FALSE, show_labels = TRUE, 
                       logFC_cutoff = log2(1.2), 
                       gset_nms = c("go_sets", "kegg_sets"), ...) {
  
  check_dots(c("id", "anal_type"), ...)
  dir.create(file.path(dat_dir, "Protein\\Volcano\\log"), recursive = TRUE, showWarnings = FALSE)

  dots <- rlang::enexprs(...)

  id <- match_call_arg(normPSM, group_pep_by)
  stopifnot(rlang::as_string(id) %in% c("prot_acc", "gene"))
  
  scale_log2r <- match_logi_gv("scale_log2r", scale_log2r)
  
  df <- rlang::enexpr(df)
  filepath <- rlang::enexpr(filepath)
  filename <- rlang::enexpr(filename)	
  show_sig <- rlang::as_string(rlang::enexpr(show_sig))
  
  stopifnot(is_logical(adjP))
  stopifnot(is_logical(show_labels))
  stopifnot(is_double(pval_cutoff))
  stopifnot(is_double(logFC_cutoff))
  
  gset_nms <- local({
    file <- file.path(dat_dir, "Calls\\anal_prnGSPA.rda")
    
    if (file.exists(file)) {
      gset_nms <- gset_nms %>% 
        .[. %in% match_call_arg(anal_prnGSPA, gset_nms)]
      
      if (is.null(gset_nms)) stop ("Unknown gene sets.")
    }
    
    return(gset_nms)
  })
  
  dots <- rlang::enexprs(...)
  fmls <- dots %>% .[grepl("^\\s*~", .)]
  dots <- dots[!names(dots) %in% names(fmls)]
  dots <- concat_fml_dots(fmls = fmls, fml_nms = fml_nms, dots = dots)
  
  reload_expts()
  
  info_anal(df = !!df, id = !!id, filepath = !!filepath, filename = !!filename, 
            scale_log2r = scale_log2r, complete_cases = complete_cases, impute_na = impute_na, 
            anal_type = "Volcano")(gset_nms = gset_nms, 
                                var_cutoff = 1000, 
                                pval_cutoff = pval_cutoff, 
                                logFC_cutoff = logFC_cutoff, 
                                !!!dots)
}

