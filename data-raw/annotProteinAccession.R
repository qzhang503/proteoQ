# library(proteoQ)
# library(UniProt.ws)
library(plyr)
library(reshape2)
library(dplyr)
library(purrr)
library(stringr)
library(magrittr)

#' Prefix form of colnames(x)[c(2, 5, ...)] for use in pipes
#'
#' \code{names_pos<-} rename the columns at the indeces of \code{pos}.
#'
#' @param x A data frame.
#' @param pos Numeric.  The index of coloumns for name change.
#' @param value Characters.  The new column names.
#' @return The data frame with new names.
#'
#' @import dplyr
#' @importFrom magrittr %>%
`names_pos<-` <- function(x, pos, value) {
  names(x)[pos] <- value
  x
}


#'Species lookup
sp_lookup <- function(species) {
  switch(species, 
         human = "hs",
         mouse = "mm",
         rat = "rn",
         fly = "dm", 
         bovine = "bt",
         dog = "cf", 
         crap = "crap", 
         stop("Unknown `species`.", Call. = FALSE)
  )    
}

sp_lookup_annotdbi <- function(species) {
  switch(species, 
         human = "Hs",
         mouse = "Mm",
         rat = "Rn",
         fly = "Dm", 
         bovine = "Bt",
         dog = "Cf", 
         stop("Unknown `species`.", Call. = FALSE)
  )    
}

#'Toxonomy lookup
taxid_lookup <- function(species) {
  switch (species,
          human = 9606, 
          mouse = 10090,
          rat = 10116, 
          fly = 7227, 
          bovine = 9913,
          dog = 9612, 
          crap = 000000, 
          stop("Unknown `species`.", Call. = FALSE)
  )
}



# a case of same GO term but different species
# add species name to GO data sets
devtools::document("c:/results/r/proteoq")
species <- "rat"
abbr_sp <- purrr::map_chr(species, sp_lookup)
filename <- paste0("go_sets_", abbr_sp)
data(package = "proteoQ", list = filename)

temp <- get(filename)
if (!grepl(paste0("^", abbr_sp), temp[1])) {
  temp <- temp %>% 
    `names<-`(paste(abbr_sp, names(.), sep = "_"))
  assign(filename, temp)
  save(list = filename, file = file.path("C:/Results/R/proteoQ/data", paste0(filename, ".RData")))
}

# new_name <- paste0("go_sets_", abbr_species_lwr)
# assign(new_name, gsets)
# save(list = new_name, file = file.path(db_path, paste0(new_name, ".rda")))



























  


# ----------------------------------------------
devtools::document(pkg  = "C:/Results/R/proteoQ")
kinase_lookup <- kinase_lookup %>% 
					dplyr::rename(
						refseq_acc = REFSEQ_PROTEIN, 
						gene = Gene, 
						kinase = Kinase, 
						classification = Classification, 
						order = Order, 
						uniprot_acc = Uniprot_accession, 
						uniprot_id = Uniprot_id, 
						prot_desc = Description
					)

kinase_lookup <- kinase_lookup %>% 
					dplyr::rename(
						kin_attr = kinase, 
						kin_class = classification, 
						kin_order = order
					)

filename <- "kinase_lookup"
save(list = filename, file = file.path("C:/Results/R/proteoQ/data", paste0(filename, ".rda")), compress = "xz")

	


	
	
	


# ----------------------------------------------
library(seqinr)
library(biomaRt)
library(UniProt.ws)
library(plyr)
library(dplyr)
library(reshape2)
library(data.table)
library(proteoQ)

#' map UNIPROT to REFSEQ
#' the "from" is always UNIPROT
map_refseq <- function(species) {
	if (species == "Hs") {
		up <- UniProt.ws(taxId = 9606)
	} else if (species == "Mm") {
		up <- UniProt.ws(taxId = 10090)
	} else if (species == "Rn") {
		up <- UniProt.ws(taxId = 10116)
	} else if (species == "Dm") {
		up <- UniProt.ws(taxId = 7227)
	}
	
	res <- UniProt.ws::select(up,
		keys = up@taxIdUniprots,
		keytype = "UNIPROTKB", 
		columns = c("REFSEQ_PROTEIN")) %>% 
		mutate(REFSEQ_PROTEIN = gsub("\\..*", "", REFSEQ_PROTEIN)) %>% 
		dplyr::rename(Uniprot_accession = UNIPROTKB)
}





# combine all .R files into one
# foo_combine_codes(filepath = file.path("C:/Results/R/proteoQ/inst/extdata/examples"))
foo_combine_codes <- function (filepath = file.path("C:/Results/R/proteoQ/R")) {
  filenames <- dir(filepath, pattern = ".R$")
  
  dir.create(file.path(filepath, "temp"))
  
  purrr::map(file.path(filepath, filenames), readLines) %>% 
    purrr::reduce(`c`, init = NULL) %>% 
    writeLines(file.path(filepath, "temp/all.R"))  
}


# zipped base protein.txt, peptide.txt and call pars
foo_zip <- function () {
  zip("my.zip", "C:/proteoQ/examples")
  unzip("./my.zip", exdir = "~/proteoQ_demo", overwrite  = FALSE)
}

