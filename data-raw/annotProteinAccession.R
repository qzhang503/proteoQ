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




#' Map uniprot or refseq to entrez
#' 
#' @examples
#' \dontrun{
#' species <- "rat"
#' res <- map_entrez(species, from = "egUNIPROT")
#' dir.create(file.path("~\\proteoQ\\dbs\\entrez\\to_unirpot"), recursive = TRUE, showWarnings = FALSE)
#' write.table(res, file.path("~\\proteoQ\\dbs\\entrez\\to_unirpot", paste0("uniprot_entrez_", species, ".txt")), sep = "\t", col.names = TRUE, row.names = FALSE)
#' 
#' res <- map_entrez(species, from = "egREFSEQ")
#' dir.create(file.path("~\\proteoQ\\dbs\\entrez\\to_refseq"), recursive = TRUE, showWarnings = FALSE)
#' write.table(res, file.path("~\\proteoQ\\dbs\\entrez\\to_refseq", paste0("refseq_entrez_", species, ".txt")), sep = "\t", col.names = TRUE, row.names = FALSE)
#' }
#' @export
map_entrez <- function(species, from) {
  sp_lookup_org <- function(species) {
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
  

  taxid <- taxid_lookup(species)
  abbr_species <- sp_lookup_org(species) 
  lwr_species <- tolower(abbr_species)
  
  if (!requireNamespace(paste("org", abbr_species, "eg.db", sep = "."), quietly = TRUE)) 
    BiocManager::install(paste("org", abbr_species, "eg.db", sep = "."))
  
  library(paste("org", abbr_species, "eg.db", sep = "."), character.only = TRUE)
  
  x <- get(paste("org", abbr_species, from, sep = "."))
  mapped_genes <- mappedkeys(x) # "from" to entrez
  
  accessions <- as.list(x[mapped_genes]) %>% 
    plyr::ldply(., rbind) %>% 
    `names_pos<-`(., 1, c("entrez")) %>% 
    `names_pos<-`(., 2:ncol(.), paste(from, 1:(length(.)-1), sep = ".")) %>% 
    mutate_at(.vars = grep("^eg", names(.)), funs(as.character)) %>% 
    melt(id = "entrez") %>% 
    filter(!is.na(entrez), !is.na(value)) %>% 
    dplyr::select(-c("variable")) 
  
  if (from == "egUNIPROT") {
    accessions <- accessions %>% dplyr::rename(uniprot_acc = value)
    path <- "~\\proteoQ\\dbs\\entrez\\to_unirpot"
    out_nm <- paste0("uniprot_entrez_", lwr_species, ".txt")
  } else if (from == "egREFSEQ") {
    accessions <- accessions %>% dplyr::rename(refseq_acc = value)
    path <- "~\\proteoQ\\dbs\\entrez\\to_refseq"
    out_nm <- paste0("refseq_entrez_", lwr_species, ".txt")
  }
  
  dir.create(file.path(path), recursive = TRUE, showWarnings = FALSE)
  write.table(accessions, file.path(path, out_nm), sep = "\t", col.names = TRUE, row.names = FALSE)

  # invisible(accessions)
}


species <- c("human", "mouse", "rat", "fly", "bovine", "dog")
purrr::walk(species, map_entrez, "egUNIPROT")
purrr::walk(species, map_entrez, "egREFSEQ")

#' Make rda files
foo <- function () {
  path <- "~\\proteoQ\\dbs\\entrez\\to_unirpot"
  filelist <- list.files(path = file.path(path), pattern = "^uniprot_entrez_.*\\.txt$")
  # path <- "~\\proteoQ\\dbs\\entrez\\to_refseq"
  # filelist <- list.files(path = file.path(path), pattern = "^refseq_entrez_.*\\.txt$")
  
  fn_prx <- gsub("(.*)\\.txt$", "\\1", filelist)
  
  purrr::walk(fn_prx, ~ {
    assign(.x, read.csv(file.path(path, paste0(.x, ".txt")), 
                        check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#"))
    save(list = .x, file = file.path(path, paste0(.x, ".rda")))
  })
}

foo()





#' Map from acc_x to acc_y, acc_z... 
#'
#' \code{map_from_to} annotates the \code{from} accession to additional fields.
#'
#' @examples
#' species <- c("human", "mouse")
#' fasta <- c("~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta", 
#'            "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_mm_2013_07.fasta") 
#' purrr::walk2(as.list(species), as.list(fasta), map_from_to)
#' @import dplyr purrr rlang seqinr
#' @importFrom magrittr %>% %$%
map_from_to <- function (species = "human", 
                         fasta = "~\\proteoQ\\dbs\\fasta\\refseq\\refseq_hs_2013_07.fasta", 
                         from = "REFSEQ", to = c("REFSEQ", "ENTREZID", "SYMBOL")) {

  stopifnot(file.exists(fasta))
  
  to <- to %>% .[! . %in% from]
  
  taxid <- taxid_lookup(species)
  abbr_species <- sp_lookup_annotdbi(species) 
  lwr_species <- tolower(abbr_species)
  
  if (!requireNamespace(paste("org", abbr_species, "eg.db", sep = "."), quietly = TRUE)) 
    BiocManager::install(paste("org", abbr_species, "eg.db", sep = "."))
  
  library(paste("org", abbr_species, "eg.db", sep = "."), character.only = TRUE)
  
  keys <- fasta %>% 
    seqinr::read.fasta(seqtype = "AA", as.string = TRUE, set.attributes = TRUE) %>% 
    names() %>% 
    gsub("\\.[^\\.]*$", "", .) %>% 
    unique()

  accessions <- AnnotationDbi::select(get(paste("org", abbr_species, "eg.db", sep = ".")), 
                                      keys = keys,
                                      columns = to,
                                      keytype = from)
  
  key_lookup <- c(
    "REFSEQ" = "refseq_acc",
    "ENTREZID" = "entrez",
    "SYMBOL" = "gene"
  )
  
  org_lookup <- c(
    "human" = "Homo sapiens",
    "mouse" = "Mus musculus", 
    "rat" = "Rattus norvegicus", 
    "fly" = "Drosophila melanogaster", 
    "bovine" = "Bos taurus", 
    "dog" = "Canis lupus familiaris", 
    "crap" = "crap"
  )
  
  # org_lookup[species]

  accessions <- accessions %>% 
    `colnames<-`(key_lookup[names(accessions)]) %>% 
    dplyr::mutate(organism = org_lookup[species])
  
  if (from == "REFSEQ") {
    out_path <- "~\\proteoQ\\dbs\\crossref\\refseq"
    out_nm <- paste0("refseq_", lwr_species, "_crossref.txt")
  }
  
  dir.create(file.path(out_path), recursive = TRUE, showWarnings = FALSE)
  write.table(accessions, file.path(out_path, out_nm), sep = "\t", col.names = TRUE, row.names = FALSE)
  
}


#' Make rda files
foo <- function () {
  path <- "~\\proteoQ\\dbs\\crossref\\refseq"
  filelist <- list.files(path = file.path(path), pattern = "^refseq_.*\\.txt$")

  fn_suffix <- gsub("^.*\\.([^.]*)$", "\\1", filelist)
  fn_prefix <- gsub("\\.[^.]*$", "", filelist)
  
  purrr::walk(fn_prefix, ~ {
    assign(.x, read.csv(file.path(path, paste0(.x, ".txt")), 
                        check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#"))
    save(list = .x, file = file.path(path, paste0(.x, ".rda")))
  })
}

foo()


# a case of same GO term but different species
# add species name to GO data sets
devtools::document("c:\\results\\r\\proteoq")
species <- "rat"
abbr_sp <- purrr::map_chr(species, sp_lookup)
filename <- paste0("go_sets_", abbr_sp)
data(package = "proteoQ", list = filename)

temp <- get(filename)
if (!grepl(paste0("^", abbr_sp), temp[1])) {
  temp <- temp %>% 
    `names<-`(paste(abbr_sp, names(.), sep = "_"))
  assign(filename, temp)
  save(list = filename, file = file.path("C:\\Results\\R\\proteoQ\\data", paste0(filename, ".RData")))
}





















# ------------------------------
#' map "from" to entrez
#' res_entrez <- map_entrez(species, from = "egUNIPROT")
map_entrez_v1 <- function(species, from) {
  sp_lookup_org <- function(species) {
    switch(species, 
           human = "Hs",
           mouse = "Mm",
           rat = "Rn",
           fly = "Dm", 
           bovine = "Bt",
           dog = "Cf", 
           crap = "crap", 
           stop("Unknown `species`.", Call. = FALSE)
    )    
  }
  
  abbr_species <- sp_lookup_org(species) 
  taxid <- taxid_lookup(species)
  
  if (!requireNamespace(paste("org", abbr_species, "eg.db", sep = "."), quietly = TRUE)) 
    BiocManager::install(paste("org", abbr_species, "eg.db", sep = "."))
  
  library(paste("org", abbr_species, "eg.db", sep = "."), character.only = TRUE)
  
  x <- get(paste("org", abbr_species, from, sep = "."))
  mapped_genes <- mappedkeys(x) # "from" to entrez
  
  accessions <- as.list(x[mapped_genes]) %>% 
    plyr::ldply(., rbind) %>% 
    `names_pos<-`(., 1, c("entrez")) %>% 
    `names_pos<-`(., 2:ncol(.), paste(from, 1:(length(.)-1), sep = ".")) %>% 
    mutate_at(.vars = grep("^eg", names(.)), funs(as.character)) %>% 
    melt(id = "entrez") %>% 
    filter(!is.na(entrez), !is.na(value)) %>% 
    dplyr::select(-c("variable")) # %>% 
  #	`colnames<-`(c("entrez", "Uniprot_accession")) 
  
  if (from == "egUNIPROT") {
    accessions <- accessions %>% dplyr::rename(Uniprot_accession = value)
  } else if (from == "egREFSEQ") {
    accessions <- accessions %>% dplyr::rename(Uniprot_id = value)
  }
  
  return(accessions)
}



#'Cross-ref file name lookup
#'
#'fn <- crossref_fn("human")
crossref_fn <- function(species) {
  switch(species, 
         hs = "~\\proteoQ\\dbs\\crossref\\hs\\uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_.tab",
         mm = "~\\proteoQ\\dbs\\crossref\\mm\\uniprot-filtered-organism__Mus+musculus+(Mouse)+[10090]_.tab",
         rn = "~\\proteoQ\\dbs\\crossref\\rn\\uniprot-filtered-organism__Rattus+norvegicus+(Rat)+[10116]_.tab",
         dm = "~\\proteoQ\\dbs\\crossref\\dm\\uniprot-dm-filtered-organism__Drosophila+melanogaster+(Fruit+fly)+[722--.tab", 
         bt = "~\\proteoQ\\dbs\\crossref\\bt\\uniprot-reviewed_yes+AND+organism__Bos+taurus+[9913]_.tab", 
         cf = "~\\proteoQ\\dbs\\crossref\\cf\\uniprot-canis+lupus+familiaris-filtered-organism__Canis+lupus+familiaris+(--.tab", 
         crap = "~\\proteoQ\\dbs\\crossref\\crap\\uniprot_craps.tab", 
         stop("Unknown `species`.", Call. = FALSE)
  )    
}


#'Re-arrange cross-reference table
#'
#'crossrefs <- arrange_crossrefs("human")
arrange_crossrefs <- function (species) {
  abbr_species <- sp_lookup(species) 
  taxid <- taxid_lookup(species)
  
  crossrefs <- crossref_fn(abbr_species) %>% 
    read.csv(check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>% 
    dplyr::select(-grep("^yourlist\\:", names(.)))
  
  refseq <- crossrefs %>% 
    dplyr::select("Entry name", "Cross-reference (RefSeq)") %$% 
    str_split(.$"Cross-reference (RefSeq)", ";", simplify = TRUE) %>% 
    data.frame() %>% 
    dplyr::bind_cols(crossrefs[, "Entry name", drop = FALSE], .) %>% 
    tidyr::gather(-"Entry name", key = Col, value = RefSeq) %>% 
    dplyr::select(-Col)	%>% 
    dplyr::filter(nchar(RefSeq) > 0) %>% 
    dplyr::filter(grepl("^NP_", RefSeq)) %>% 
    dplyr::mutate(RefSeq = gsub("\\s+\\[.*\\]\\s*", "", .$RefSeq)) %>% 
    dplyr::mutate(RefSeq = gsub("\\.\\d*", "", .$RefSeq))
  
  crossrefs <- crossrefs %>% 
    dplyr::select(-"Cross-reference (RefSeq)") %>% 
    dplyr::mutate("Gene names" = gsub("\\s+.*", "", .$"Gene names")) %>% 
    dplyr::full_join(refseq, by = "Entry name")
  
  rm(refseq)
  
  if (species == "crap") {
    res_entrez <- data.frame(Entrez = NA, Entry = unique(crossrefs$Entry))
  } else {
    res_entrez <- map_entrez(species, from = "egUNIPROT") %>% 
      dplyr::rename("Entry" = "Uniprot_accession", "Entrez" = "entrez")
  }
  
  crossrefs <- left_join(crossrefs, res_entrez, by = "Entry")
  # "uniprot_acc" "uniprot_id"  "status"      "prot_desc"   "gene"        "organism"    "length"      "refseq_acc"  "entrez" 
  
  
  # entrez <- read.csv(file.path("~\\proteoQ\\dbs\\crossref", abbr_species, paste0("entrez_", abbr_species, ".txt")), 
  #                    check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>% 
  #   dplyr::rename("Entry name" = "From") %>% 
  #   dplyr::rename("Entrez" = "To")
  
  # crossrefs <- crossrefs %>% 
  #   dplyr::left_join(entrez, by = "Entry name")
  
  filename <- paste0("prn_annot_", abbr_species) %>% tolower()
  assign(filename, 
         crossrefs %>% 
           dplyr::rename(uniprot_acc = Entry, 
                         uniprot_id = "Entry name", 
                         status = "Status", 
                         prot_desc = "Protein names", 
                         gene = "Gene names", 
                         organism = Organism, 
                         length = Length, 
                         refseq_acc = RefSeq, 
                         entrez = Entrez))
  
  write.table(get(filename), 
              file.path("~\\proteoQ\\dbs\\crossref", abbr_species, paste0("uniprot_crossref_", abbr_species, ".txt")), 
              sep = "\t", col.names = TRUE, row.names = FALSE)
  save(list = filename, file = file.path("C:\\Results\\R\\proteoQ\\data", paste0(filename, ".RData")))
}


# crossrefs <- arrange_crossrefs("dog")








# -----------------------------------------
#' Make rda files
foo <- function () {
  path <- "~\\proteoQ\\dbs\\entrez\\to_unirpot"
  filelist <- list.files(path = file.path(path), pattern = "^uniprot_entrez_.*\\.txt$")
  # path <- "~\\proteoQ\\dbs\\entrez\\to_refseq"
  # filelist <- list.files(path = file.path(path), pattern = "^refseq_entrez_.*\\.txt$")
  
  fn_prx <- gsub("(.*)\\.txt$", "\\1", filelist)

  purrr::walk(fn_prx, ~ {
    assign(.x, read.csv(file.path(path, paste0(.x, ".txt")), 
           check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#"))
    save(list = .x, file = file.path(path, paste0(.x, ".rda")))
  })
}

foo()








  


# ----------------------------------------------
devtools::document(pkg  = "C:\\Results\\R\\proteoQ")
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
save(list = filename, file = file.path("C:\\Results\\R\\proteoQ\\data", paste0(filename, ".RData")))

	



library(GSVAdata)
data(c2BroadSets, package = "GSVAdata")
c2_msig_hs <- map(c2BroadSets, ~ .@geneIds)
setNames <- map(c2BroadSets, ~ .@setName)
names(c2_msig_hs) <- setNames
save(list = "c2_msig_hs", file = file.path("C:\\Results\\R\\proteoQ\\data", paste0("c2_msig_hs", ".RData")))

	
	
	
	

	
	
	
	
	
	


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

#' map REFSEQ to gene
map_refseq_gene <- function(species) {
	if (species == "Hs") {
		mart <- useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
		symbol <- "hgnc_symbol"
	} else if (species == "Mm") {
		mart <- useMart(biomart = 'ensembl', dataset = 'mmusculus_gene_ensembl')
		symbol <- "mgi_symbol"
	} else if (species == "Rn") {
		mart <- useMart(biomart = 'ensembl', dataset = 'rnorvegicus_gene_ensembl')
		symbol <- "rgd_symbol"
	} else if (species == "Dm") {
		mart <- useMart(biomart = 'ensembl', dataset = 'dmelanogaster_gene_ensembl')
		symbol <- "external_gene_name"
	}
	
	# listAttributes(mart) %>% filter(grepl("gene", name, ignore.case = TRUE))

	getBM(attributes = c("refseq_peptide", symbol), filters = "refseq_peptide", values = accessions$prot_acc, mart = mart)
}
# res_gene <- map_refseq_gene(species) %>% `names_pos<-`(1:2, c("Uniprot_id", "Gene"))

#' map "from" to entrez
#' from = c("egUNIPROT", "egREFSEQ")
map_entrez <- function(species, from) {
  library(paste("org", species, "eg.db", sep = "."), character.only = TRUE)
  
  x <- get(paste("org", species, from, sep = "."))
  mapped_genes <- mappedkeys(x) # "from" to entrez
  
  accessions <- as.list(x[mapped_genes]) %>% 
    plyr::ldply(., rbind) %>% 
    `names_pos<-`(., 1, c("entrez")) %>% 
    `names_pos<-`(., 2:ncol(.), paste(from, 1:(length(.)-1), sep = ".")) %>% 
    mutate_at(.vars = grep("^eg", names(.)), funs(as.character)) %>% 
    melt(id = "entrez") %>% 
    filter(!is.na(entrez), !is.na(value)) %>% 
    dplyr::select(-c("variable")) # %>% 
  #	`colnames<-`(c("entrez", "Uniprot_accession")) 
  
  if (from == "egUNIPROT") {
    accessions <- accessions %>% dplyr::rename(Uniprot_accession = value)
  } else if (from == "egREFSEQ") {
    accessions <- accessions %>% dplyr::rename(Uniprot_id = value)
  }
  
  return(accessions)
}



#' subset fasta by species
db <- read.fasta(paste0("C:\\Results\\R\\proteoQ\\data\\SwissProt_2014_07.fasta"), seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
for (species in c("Hs", "Mm", "Rn", "Dm")) {
	if (species == "Hs") {
		db_sub <- db[grepl("_HUMAN$", names(db))]
	} else if (species == "Mm") {
		db_sub <- db[grepl("_MOUSE$", names(db))]
	} else if (species == "Rn") {
		db_sub <- db[grepl("_RAT$", names(db))]
	} else if (species == "Dm") {
		db_sub <- db[grepl("_DROME$", names(db))]
	}
	
	write.fasta(sequences = db_sub, names = seqinr::getAnnot(db_sub), nbchar = 80, 
		file.out = file.path("C:\\Results\\R\\proteoQ\\data", paste("SwissProt", species, "2014_07.fasta", sep = "_")))
		
	# temp <- read.fasta(file.path("C:\\Results\\R\\proteoQ\\data", paste("SwissProt", species, "2014_07.fasta", sep = "_")), seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
	# temp_x <- read.fasta(file.path("C:\\Results\\R\\proteoQ\\data\\original\\", paste("SwissProt", species, "2014_07.fasta", sep = "_")), seqtype = "AA", as.string = TRUE, set.attributes = TRUE)
}





SwissProt_annot <- NULL
for (species in c("", "Mm", "Rn", "Dm")) {
	# obtain entrez from UNIPROT
	res_entrez <- map_entrez(species, from = "egUNIPROT")

	# obtain REFSEQ from UNIPROT
	res_refseq <- map_refseq(species)
		
	# ---------------------------------------------------------------
	# obtain Description and Gene from fasta
	db <- read.fasta(paste0("C:\\Results\\R\\proteoQ\\data\\SwissProt_", species, "_2014_07.fasta"), 
		seqtype = "AA", as.string = TRUE, set.attributes = TRUE)

	db_annot <- seqinr::getAnnot(db) 
	
	# protein description from annotation
	db_desc <- db_annot %>% 
		unlist() %>% 
		plyr::ldply(data.frame, check.names = FALSE) %>% 
		`colnames<-`("Description")	
	
	# gene names from annotation
	db_gene <- data.frame(do.call('rbind', strsplit(as.character(db_annot), 'GN=', fixed = TRUE))) %>% 
		dplyr::select(c(2)) %>% 
		`colnames<-`("Gene") %>% 
		mutate(Gene = gsub("\\s+.*|.*\\|", "", Gene))		
	
	# accessions, description and gene names from fasta names
	accessions <- data.frame(do.call('rbind', strsplit(as.character(names(db)), '|', fixed = TRUE))) %>% 
		dplyr::select(c(2:3)) %>% 
		dplyr::rename(Uniprot_accession = X2, Uniprot_id = X3) %>% 
		bind_cols(db_desc, db_gene) %>% 
		mutate(Description = sub(".*? (.+)", "\\1", Description))

	rm(db, db_annot, db_desc, db_gene)
	# ---------------------------------------------------------------


	# combined fasta, entrez and refseq results
	assign(paste("SwissProt", species, "annot", sep = "_"), 
		list(accessions, res_entrez, res_refseq) %>% 
		purrr::reduce(left_join, by = "Uniprot_accession")
	)
	
	rm(accessions, res_entrez, res_refseq)
	

	# save results for individual species
	saveRDS(get(paste("SwissProt", species, "annot", sep = "_")), 
		file.path("C:\\Results\\R\\proteoQ\\data", paste("SwissProt", species, "annot.rds", sep = "_")))

	# combine results over all species 
	SwissProt_annot <- rbind(SwissProt_annot, get(paste("SwissProt", species, "annot", sep = "_")))
}

# save results for all species
saveRDS(SwissProt_annot, file.path("C:\\Results\\R\\proteoQ\\data", paste("SwissProt", "annot.rds", sep = "_")))





