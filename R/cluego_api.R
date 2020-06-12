#' ClueGO API
#' 
#' \code{cluego} sends \code{\link{anal_prnTrend}} findings in, e.g.
#' \code{Protein_Trend_[...].txt} etc. to
#' \code{\href{http://apps.cytoscape.org/apps/cluego}{ClueGO}} via API.
#' @param df2 Character string; the name of a secondary data file from
#'   \code{\link{anal_prnTrend}}, such as Protein_Trend_Z_nclust6.txt etc.
#' @param species Named character string . The default is \code{c(human = "Homo
#'   Sapiens")}. In this example, \code{human} is the value that can be found
#'   under the \code{species} column in a Trend result table and \code{"Homo
#'   Sapiens"} is the corresponding species name in ClueGO.
#' @param n_clust Numeric vector; the cluster ID(s) corresponding to
#'  \code{\link{anal_prnTrend}} for visualization. 
#' @import httr RJSONIO xml2 rlang dplyr purrr 
#' @importFrom magrittr %>%
#' @examples \donttest{
#' # Make sure CytoSpace is opened, and 
#' # yFiles layouts and ClueGO plug-in are installed.
#' 
#' # Human
#' cluego(
#'   df2 = Protein_Trend_Z_nclust5.txt, 
#'   species = c(human = "Homo Sapiens"), 
#'   n_clust = c(3, 5)
#' )
#' 
#' # Mouse
#' cluego(
#'   df2 = Protein_Trend_Z_nclust5.txt, 
#'   species = c(mouse = "Mus Musculus"), 
#'   n_clust = c(3:4)
#' )
#' }
cluego <- function (df2 = "Protein_Trend_Z_nclust5.txt", species = c("human" = "Homo Sapiens"), 
                    n_clust = 3) {
  
  df2 <- rlang::as_string(rlang::enexpr(df2))
  filepath <- file.path(dat_dir, "Protein", "Trend")

  if (length(species) > 1) {
    stop("Only one species at a time.", call. = FALSE)
  }
  
  if (!rlang::is_integerish(n_clust)) {
    stop("`n_clust` needs to be integers.", call. = FALSE)
  }
  
  df <- tryCatch(
    readr::read_tsv(file.path(filepath, df2), col_types = cols(group = col_factor())), 
    error = function(e) NA
  )
  
  if (is.null(dim(df))) stop("File `", df2, "` not found.", call. = FALSE)
  
  df <- df %>% 
    dplyr::filter(.data$species == names(.env$species), cluster %in% n_clust)
  
  if (nrow(df) == 0) {
    stop("Zero data entries at species `", species, 
         "` and clusters ", purrr::reduce(n_clust, paste, sep = ", "), 
         call. = FALSE)
  }
  
  n_clust <- unique(df$cluster)

  # 0 --- start up
  cytoscape_base_url <- paste("http://localhost:1234", "/v1", sep = "")
  cluego_base_url <- paste(cytoscape_base_url, "apps", "cluego", "cluego-manager", sep = "/")
  
  response <- tryCatch(
    httr::POST(url = paste(cytoscape_base_url, "apps", "cluego", "start-up-cluego", sep = "/"), 
               encode = "json"), 
    error = function(e) NA
  )
  
  if (all(is.na(response))) stop("Start Cytoscape first!", call. = FALSE)
  
  if (grepl("HTTP ERROR 500",response)) {
    message("Initiating ClueGO...")
    Sys.sleep(2)
  }
  
  # 1 --- "organisms"
  response <- httr::PUT(
    url = paste(cluego_base_url, "organisms", "set-organism", 
                URLencode(species), sep = "/"), encode = "json"
  )

  # 2 --- "max-input-panel"
  max_input_panel_number <- length(n_clust)
  response <- httr::PUT(
    url = paste(cluego_base_url, "cluster", "max-input-panel", max_input_panel_number, sep = "/"), 
    encode = "json"
  )

  # 3 --- "set-analysis-properties"
  input_panel_indexes <- as.character(seq_along(n_clust)) 
  
  node_shape <- "Ellipse"
  cluster_color <- "#ff0000"
  min_number_of_genes_per_term <- 3
  min_percentage_of_genes_mapped <- 4
  no_restrictions <- FALSE  
  
  purrr::map(input_panel_indexes, ~ {
    response <- httr::PUT(
      url = paste(cluego_base_url, "cluster", "set-analysis-properties", 
                  .x, node_shape, URLencode(cluster_color, reserved = TRUE), 
                  min_number_of_genes_per_term, min_percentage_of_genes_mapped, 
                  no_restrictions, sep = "/"), 
      encode = "json"
    )
  })
  
  #  4 --- "upload-ids-list"
  purrr::map2(input_panel_indexes, n_clust, ~ {
    genes <- df %>% 
      dplyr::filter(cluster == .y) %>% 
      .[[match_call_arg(normPSM, group_pep_by)]] %>% 
      unique() %>% 
      RJSONIO::toJSON()

    response <- httr::PUT(
      url = paste(cluego_base_url, "cluster", "upload-ids-list", URLencode(.x), 
                  sep = "/"), 
      body = genes, encode = "json", httr::content_type_json()
    )
  })

  # 5 --- select visual style
  if (max_input_panel_number > 1) {
    visual_style <- "ShowClusterDifference"
  } else {
    visual_style = "ShowGroupDifference"
  }
  
  response <- httr::PUT(
    url = paste(cluego_base_url, "cluster", "select-visual-style", visual_style, sep = "/"), 
    encode = "json"
  )
  
  # 6 --- "set-ontologies"
  selected_ontologies <- RJSONIO::toJSON(c("3;Ellipse", "8;Triangle", "9;Rectangle"))
  
  response <- httr::PUT(
    url = paste(cluego_base_url, "ontologies", "set-ontologies", sep = "/"), 
    body = selected_ontologies, encode = "json", httr::content_type_json()
  )
  
  # 7 --- "set-min-max-levels"
  min_level <- 5
  max_level <- 6
  all_levels <- FALSE  
  
  response <- httr::PUT(
    url = paste(cluego_base_url, "ontologies", "set-min-max-levels", 
                min_level, max_level, all_levels, sep = "/"), 
    encode = "json"
  )
  
  # # 3.1 Get all available Ontologies
  # response <- httr::GET(paste(cluego_base_url, "ontologies", "get-ontology-info", sep = "/"))
  # content(response)

  # 3.2 Select Evidence Codes
  evidence_codes <- RJSONIO::toJSON(c("All"))
  
  response <- httr::PUT(
    url = paste(cluego_base_url, "ontologies", "set-evidence-codes", sep = "/"), 
    body = evidence_codes, encode = "json", httr::content_type_json()
  )
  
  # # 3.3 Get all available Evidence Codes
  # response <- httr::GET(paste(cluego_base_url, "ontologies", "get-evidence-code-info",sep = "/"))
  
  # 3.5 Use GO term significance cutoff
  p_value_cutoff <- 0.05
  use_significance_cutoff <- TRUE
  
  response <- httr::PUT(
    url = paste(cluego_base_url, "ontologies", use_significance_cutoff, p_value_cutoff, sep = "/"), 
    encode = "json"
  )
  
  # 3.6 Use GO term fusion
  use_go_term_fusion <- TRUE
  
  response <- httr::PUT(
    url = paste(cluego_base_url, "ontologies", use_go_term_fusion, sep = "/"), 
    encode = "json"
  )
  
  # 3.7 Set statistical parameters
  enrichment_type <- "Enrichment/Depletion (Two-sided hypergeometric test)" 
  multiple_testing_correction <- TRUE
  use_mid_pvalues <- FALSE
  use_doubling <- FALSE
  
  response <- httr::PUT(
    url = paste(cluego_base_url, "stats", enrichment_type, multiple_testing_correction, 
                use_mid_pvalues, use_doubling, sep = "/"), encode = "json")
  
  # 3.8 Set the Kappa Score level
  kappa_score <- 0.4
  
  response <- httr::PUT(
    url = paste(cluego_base_url, "ontologies", "set-kappa-score-level", kappa_score, sep = "/"), 
    encode = "json")
  
  # # 3.9 Set grouping parameters
  do_grouping <- TRUE
  coloring_type <- "Random"
  group_leading_term <- "Highest Significance" 
  grouping_type <- "Kappa Score"
  init_group_size <- 1
  perc_groups_for_merge <- 50
  perc_terms_for_merge <- 50
  
  response <- httr::PUT(
    url = paste(cluego_base_url, "grouping", do_grouping, coloring_type, group_leading_term, 
                grouping_type, init_group_size, perc_groups_for_merge, perc_terms_for_merge, 
                sep = "/"), encode = "json")

  # 8 --- "analysis"
  analysis_name <- paste(df2, reduce(sort(n_clust), paste, sep = "-"))
  
  response <- httr::GET(paste(cluego_base_url, URLencode(analysis_name), 
                              URLencode("Cancel and refine selection"), 
                              sep = "/"))
}

