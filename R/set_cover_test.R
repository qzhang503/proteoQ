set_cover_test <- function() {
  df <- tryCatch(read.csv(file.path(dat_dir, "Protein", "Protein.txt"), check.names = FALSE, header = TRUE, sep = "\t",
                          comment.char = "#"), error = function(e) NA)
  
  gset_nms = "go_sets"
  species = "human"
  gsets <- load_dbs(gset_nms, species)
  
  gsets <- gsets %>% 
    .[!grepl("molecular_function$", names(.))] %>% 
    .[!grepl("cellular_component$", names(.))] %>% 
    .[!grepl("biological_process$", names(.))]
    
  
  gs <- purrr::map2(gsets, names(gsets), ~ {
    .x %>% 
      tibble() %>% 
      dplyr::mutate(set = .y)
  }) %>% 
    bind_rows() %>% 
    `names<-` (c("id", "set")) %>% 
    dplyr::filter(id %in% unique(df$entrez)) %>% 
    dplyr::select(c("set", "id"))
  
  # filter by sets from GSPA...
  pval_cutoff = .05
  logFC_cutoff = 1.2
  subdir = "W2_bat"
  gsea_res <- do.call(rbind,
                      purrr::map(
                        list.files(path = file.path(dat_dir, "Protein\\GSPA", subdir), 
                                   pattern = "Protein_GSPA_.*\\.csv$", , full.names = TRUE),
                        readr::read_csv
                      )) %>% 
    dplyr::filter(.$term %in% names(gsets))
  
  gsea_res <- gsea_res %>% 
    dplyr::filter(p_val <= pval_cutoff, abs(log2fc) >= log2(logFC_cutoff))
  
  gs <- gs %>% 
    dplyr::filter(set %in% gsea_res$term)
  
  res <- RcppGreedySetCover::greedySetCover(gs, FALSE)
}