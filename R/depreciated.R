#'Downloads STRING databases
#'
#'@param species Character string; the species. The currently supported species
#'  include \code{human, mouse, rat, fly, cow, dog}. The default is
#'  \code{human}.
#'@param overwrite Logical; if TRUE, overwrite the downloaded database(s). The
#'  default is FALSE.
#'@inheritParams anal_prnString
#'@inheritParams prepString
#'@import rlang dplyr magrittr purrr fs downloader
#'@seealso \code{\link{anal_prnString}} for protein-protein interaction
#'  networks.
#'@export
dl_stringdbs <- function(species = "human", db_path = "~\\proteoQ\\dbs\\string", overwrite = FALSE) {
  species <- rlang::as_string(rlang::enexpr(species))
  db_path <- create_db_path(db_path)
  
  stop("`dl_stringdbs` depreciated; instead use `prepString`.", call. = FALSE)
  
  urls <- switch(species, 
                 human = c(
                   "9606.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz",
                   "9606.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/9606.protein.aliases.v11.0.txt.gz", 
                   "9606.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/9606.protein.info.v11.0.txt.gz"
                 ), 
                 mouse = c(
                   "10090.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/10090.protein.links.full.v11.0.txt.gz",
                   "10090.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/10090.protein.aliases.v11.0.txt.gz",
                   "10090.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/10090.protein.info.v11.0.txt.gz"
                 ), 
                 rat = c(
                   "10116.protein.links.full.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.links.full.v11.0/10116.protein.links.full.v11.0.txt.gz",
                   "10116.protein.aliases.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.aliases.v11.0/10116.protein.aliases.v11.0.txt.gz", 
                   "10116.protein.info.v11.0.txt.gz" = "https://stringdb-static.org/download/protein.info.v11.0/10116.protein.info.v11.0.txt.gz"
                 ), 
                 "unknown"
  )
  
  for(i in seq_along(urls)) {
    url <- urls[[i]]
    fn_zip <- names(urls)[[i]]
    fn_tsv <- gsub("\\.gz$", "", fn_zip)
    filepath <- file.path(db_path, fn_zip)
    
    if ((!file.exists(filepath)) | overwrite) {
      downloader::download(url, filepath, mode = "wb")
      con <- gzfile(path.expand(filepath))
      
      if (grepl("protein\\.links", filepath)) {
        read.csv(con, sep = " ", check.names = FALSE, header = TRUE, comment.char = "#") %>% 
          write.table(file.path(db_path, fn_tsv), sep = "\t", col.names = TRUE, row.names = FALSE)
        # close(con)
      } else if (grepl("protein\\.info", filepath)) {
        local({
          temp <- readLines(con)
          for (idx in seq_along(temp)) {
            temp[idx] <- gsub("^(.*)\t[^\t].*$", "\\1", temp[idx])
          }
          temp[1] <- "protein_external_id\tpreferred_name\tprotein_size"
          writeLines(temp, file.path(db_path, fn_tsv))          
        })
        # close(con)
      } else if (grepl("protein\\.alias", filepath)) {
        local({
          temp <- read.csv(con, sep = "\t", check.names = FALSE, header = TRUE, comment.char = "#")
          
          col_nms <- c("string_protein_id", "alias", "source")
          first_row <- names(temp) %>% 
            data.frame() %>% 
            t() %>% 
            `colnames<-`(col_nms)
          
          temp %>% 
            `colnames<-`(col_nms) %>% 
            dplyr::mutate_all(as.character) %>% 
            rbind(first_row, .) %>% 
            write.table(file.path(db_path, fn_tsv), sep = "\t", col.names = TRUE, row.names = FALSE)          
        })
        # close(con)
      }
      
    }
  }
  
}

