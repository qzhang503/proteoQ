library(dplyr)
library(stringr)
library(magrittr)

for (species in c("hs", "mm", "rn", "dm")) {
	if(species == "hs") {
		fn <- "C:\\Results\\DB\\Swissprot\\20190306\\hs\\uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_.tab"
	} else if(species == "mm") {
		fn <- "C:\\Results\\DB\\Swissprot\\20190306\\mm\\uniprot-filtered-organism__Mus+musculus+(Mouse)+[10090]_.tab"
	} else if(species == "rn") {
		fn <- "C:\\Results\\DB\\Swissprot\\20190306\\rn\\uniprot-filtered-organism__Rattus+norvegicus+(Rat)+[10116]_.tab"
	} else if(species == "dm") {
		fn <- "C:\\Results\\DB\\Swissprot\\20190306\\dm\\uniprot-dm-filtered-organism__Drosophila+melanogaster+(Fruit+fly)+[722--.tab"
	}
	
	crossrefs <- read.csv(fn, check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#")
	
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
	# write.table(crossrefs, file.path("C:\\Results\\DB\\Swissprot\\20190306\\Dm\\temp.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
	
	entrez <- read.csv(file.path("C:\\Results\\DB\\Swissprot\\20190306", species, paste0("Entrez_", species, ".txt")) , check.names = FALSE, header = TRUE, sep = "\t", comment.char = "#") %>% 
		dplyr::rename("Entry name" = "From") %>% 
		dplyr::rename("Entrez" = "To")
	
	crossrefs <- crossrefs %>% 
		dplyr::left_join(entrez, by = "Entry name")

	filename <- paste0("prn_annot_", species) %>% tolower()
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
	
	write.table(get(filename), file.path("C:\\Results\\DB\\Swissprot\\20190306", species, paste0("UniProt_CrossRef_", species, ".txt")), sep = "\t", col.names = TRUE, row.names = FALSE)
	save(list = filename, file = file.path("C:\\Results\\R\\proteoQ\\data", paste0(filename, ".RData")))
}


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





