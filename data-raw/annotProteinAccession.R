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


# --- craps lookup --- 
foo_craps <- function () {
  data(package = "proteoQ", prn_annot_crap)
  
  prn_annot_crap[] <- prn_annot_crap %>% 
    purrr::map(~ {
      if (is.factor(.x)) as.character(.x) else (.x)
    })

  save(prn_annot_crap, file = "~/proteoQ/prn_annot_crap.rda", compress = "xz")
}




# ----------------------------------------------
make_kin_lookup <- function (include_mouse = FALSE) {
  to_title_case <- function(x) {
    paste0(substr(x, 1, 1), substr(x, 2, nchar(x)) %>% tolower())
  }
  
  
  load_expts("~/proteoQ/examples_mascot")
  
  fasta = c("~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2014_07.fasta", 
            "~/proteoQ/dbs/fasta/uniprot/uniprot_mm_2014_07.fasta")
  
  kin_order <- c("TK" = 1, "TKL" = 2, "STE" = 3, "CK1" = 4, "AGC" = 5, "CAMK" = 6, "CMGC" = 7, 
                 "Atypical" = 8, "Other" = 9, "RGC" = 10, "Unclassified" = 11, "FAM20" = 12, 
                 "Lipid" = 13, "Metabolic" = 14)
  
  ## --- three columns in `df`
  # prot_acc	kin_class	gene
  # P31749	AGC 	AKT1
  # P31751	AGC 	AKT2
  # Q9Y243	AGC 	AKT3
  
  df <- readr::read_tsv(file.path("~/proteoQ/kin_lookup_08282020.txt")) %>% 
    dplyr::mutate(kin_class = gsub("Lipid Kinase", "Lipid", kin_class), 
                  kin_class = gsub("Metabolic Enzyme", "Metabolic", kin_class))
  
  # mouse
  if (include_mouse) {
    df <- local({
      df_mm <- df %>% 
        dplyr::select(-prot_acc) %>% 
        dplyr::mutate(gene = to_title_case(gene))
      
      annot_from_to(abbr_species = "Mm", 
                    keys = unique(df_mm$gene), 
                    from = "SYMBOL", 
                    to = "UNIPROT") %>% 
        dplyr::filter(!is.na(UNIPROT)) %>% 
        dplyr::rename(gene = SYMBOL, prot_acc = UNIPROT) %>% 
        dplyr::mutate(gene2 = toupper(gene)) %>% 
        dplyr::left_join(df[, c("kin_class", "gene")], by = c("gene2" = "gene")) %>% 
        dplyr::select(prot_acc, kin_class, gene) 
    }) %>% 
      rbind(df, .)
  }

  # [1] "refseq_acc"  "gene"        "kin_attr"    "kin_class"   "kin_order"   "uniprot_acc" "uniprot_id" 
  # [8] "prot_desc"   "entrez"     

  df <- df %>% 
    dplyr::mutate(kin_order = .env$kin_order[.data$kin_class]) %>% 
    annotPrn(fasta, entrez = NULL) %>% 
    dplyr::rename(uniprot_acc = prot_acc) %>% 
    dplyr::mutate(kin_attr = TRUE)
  
  # refseq
  refseq_hs <- local({
    up <- UniProt.ws(taxId = 9606)
    
    refseq_hs <- UniProt.ws::select(up,
                       keys = up@taxIdUniprots %>% .[. %in% df$uniprot_acc],
                       keytype = "UNIPROTKB", 
                       columns = c("REFSEQ_PROTEIN")) %>% 
      mutate(REFSEQ_PROTEIN = gsub("\\..*$", "", REFSEQ_PROTEIN)) %>% 
      dplyr::rename(uniprot_acc = UNIPROTKB, refseq_acc = REFSEQ_PROTEIN)
  })
  
  if (include_mouse) {
    refseq_mm <- local({
      up <- UniProt.ws(taxId = 10090)
      
      UniProt.ws::select(up,
                         keys = up@taxIdUniprots %>% .[. %in% df$uniprot_acc],
                         keytype = "UNIPROTKB", 
                         columns = c("REFSEQ_PROTEIN")) %>% 
        mutate(REFSEQ_PROTEIN = gsub("\\..*$", "", REFSEQ_PROTEIN)) %>% 
        dplyr::rename(uniprot_acc = UNIPROTKB, refseq_acc = REFSEQ_PROTEIN)
    })
    
    refseq <- rbind(refseq_hs, refseq_mm)
  } else {
    refseq <- refseq_hs
  }

  kinase_lookup <- dplyr::left_join(df, refseq, by = "uniprot_acc") %>% 
    dplyr::select(refseq_acc, gene, kin_attr, kin_class, kin_order, 
                  uniprot_acc, uniprot_id, prot_desc, entrez)

  save(kinase_lookup, file = "~/proteoQ/kinase_lookup.rda", compress = "xz")
}

	
#' PSP kinase-substrate
#'
#' @param db_nms  Character string(s) to the name(s) of PSP database(s) with
#'   prepended directory path(s). Users need to download the kinase-substrate
#'   table, e.g. \code{Kinase_Substrate_Dataset.txt} directly from the PSP
#'   website and supply the corresponding file path(s) and name(s). Currently
#'   assume single database file.
#' @param match_orgs Logical; if TRUE, matches the organism between kinases and
#'   their acting substrates. The default is TRUE.
#' @param type A character string for pathway IDs. The default is \code{gene}.
KinSub <- function (db_nms = "~/proteoQ/dbs/psp/Kinase_Substrate_Dataset.txt", 
                    type = c("gene", "entrez"), match_orgs = TRUE) {

  dbs <- readr::read_tsv(file.path(db_nms), 
                         col_types = cols(GENE = col_character(), 
                                          KIN_ORGANISM = col_character(),
                                          SUB_GENE_ID = col_character(), 
                                          SUB_GENE = col_character(),
                                          SUB_ORGANISM = col_character())) %>% 
    { if (match_orgs) dplyr::filter(., KIN_ORGANISM == SUB_ORGANISM) else . } 
  
  type <- rlang::enexpr(type)
  if (length(type) > 1) {
    type <- "gene"
  } else {
    type <- rlang::as_string(type)
    stopifnot(type %in% c("gene", "entrez"), 
              length(type) == 1)
  }
  
  if (type == "entrez") {
    filelist <- c("uniprot_entrez_hs", "uniprot_entrez_mm", "uniprot_entrez_rn")
    suppressWarnings(data(package = "proteoQ", list = filelist))
    
    up_en <- purrr::map(filelist, ~ try(get(.x), silent = TRUE)) %>% 
      dplyr::bind_rows() %>% 
      dplyr::select(c("uniprot_acc", "entrez")) # %>% 
      # dplyr::filter(!duplicated(uniprot_acc)) # ok: one uniprot_acc, multiple entrez ids
    
    suppressWarnings(rm(list = filelist, envir = .GlobalEnv))
    
    dbs <- dbs %>% 
      dplyr::left_join(up_en, by = c("KIN_ACC_ID" = "uniprot_acc")) %>% 
      dplyr::filter(!is.na(entrez)) %>% # Non c("human", "mouse", "rat") has no entrez
      dplyr::rename(GENE_ID = entrez) %>% 
      reloc_col_before("GENE_ID", "KIN_ACC_ID") 
    
    # some redundancy of the same GENE different GENE_ID
    # A tibble: 3 x 6
    # GENE  KINASE   GENE_ID KIN_ACC_ID KIN_ORGANISM SUBSTRATE
    #   1 Dyrk2 DYRK2      69181 Q5U4C9     mouse        NDEL1    
    #   2 Pak2  PAK2   100910732 Q64303     rat          MEK1     
    #   3 Pak2  PAK2       29432 Q64303     rat          MEK1  
    
    # list by GENE_ID (entrez)
    out <- dbs %>% 
      split(.$KIN_ORGANISM, drop = TRUE) %>% 
      purrr::map(~ split(.x, .$GENE_ID, drop = TRUE)) %>% 
      purrr::map(~ .x %>% purrr::map(`[[`, "SUB_GENE_ID")) %>% 
      `names<-`(gsub("human", "hs", names(.))) %>% 
      `names<-`(gsub("mouse", "mm", names(.))) %>% 
      `names<-`(gsub("rat", "rn", names(.))) %>% 
      purrr::iwalk(~ {
        fn <- paste0("kinsub_", .y)
        assign(fn, .x)
        save(list = fn, file = file.path(dat_dir, paste0(fn, ".rda")), compress = "xz")
      })
  } else {
    # list by GENE
    out <- dbs %>% 
      split(.$KIN_ORGANISM, drop = TRUE) %>% 
      purrr::map(~ split(.x, .$GENE, drop = TRUE)) %>% 
      purrr::map(~ .x %>% purrr::map(`[[`, "SUB_GENE_ID")) %>% 
      `names<-`(gsub("human", "hs", names(.))) %>% 
      `names<-`(gsub("mouse", "mm", names(.))) %>% 
      `names<-`(gsub("rat", "rn", names(.))) %>% 
      purrr::iwalk(~ {
        fn <- paste0("kinsub_", .y)
        assign(fn, .x)
        save(list = fn, file = file.path(dat_dir, paste0(fn, ".rda")), compress = "xz")
      })
  }
}


	


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
	  dplyr::mutate(REFSEQ_PROTEIN = gsub("\\..*", "", REFSEQ_PROTEIN)) %>% 
	  dplyr::rename(Uniprot_accession = UNIPROTKB)
	
	## ---- easier new approach ----
	# human
	dat_dir <- "~/proteoQ"
	fasta <- read_fasta(acc_pattern = ">..\\|([^\\|]+)\\|.*", 
	                    file = "~/proteoQ/dbs/fasta/uniprot/uniprot_hs_2020_05.fasta")
	
	# (1) uniprot_acc ~ gene
	up_gn <- data.frame(prot_desc = purrr::map_chr(fasta, attr, "header")) %>% 
	  dplyr::mutate(uniprot_acc = gsub(">..\\|([^\\|]+)\\|.*", "\\1", prot_desc)) %>% 
	  dplyr::mutate(gene = gsub("^.*GN=(\\S+)\\s*.*", "\\1", prot_desc)) %>% 
	  dplyr::select(-prot_desc) 
	
	# (2) refseq_acc ~ uniprot_acc
	up_gn %>% 
	  dplyr::select(uniprot_acc) %>% 
	  readr::write_tsv(file.path(dat_dir, "uniprot_acc.txt"))

	# load "uniprot_acc.txt" to uniprot web site
	# ...
	# save results: "uniprot_acc_to_refseq_hs.txt"
	
	up_rs <- readr::read_tsv("~/proteoQ/uniprot_acc_to_refseq_hs.txt") %>% 
	  dplyr::rename(uniprot_acc = From, refseq_acc = To) %>% 
	  dplyr::mutate(refseq_acc = gsub("\\.[0-9]*$", "", refseq_acc))

	# (3) refseq_acc ~ gene
	rs_gn <- up_rs %>% 
	  dplyr::left_join(up_gn, by = "uniprot_acc") %>% 
	  dplyr::select(-uniprot_acc) %>% 
	  dplyr::filter(!duplicated(refseq_acc))
}


foo_spec_lib <- function () {
  library(magrittr)
  library(dplyr)
  
  dat_dir <- "Z:/Skeath/2017/SKEATH-002/Mass Spectrometry/QE2/QZ"
  
  df <- readr::read_tsv(file.path(dat_dir, "lab_contams.txt"))
  
  fasta <- proteoQ::read_fasta(file = "Y:/QZHANG/dbs/fasta/uniprot/uniprot_20201007.fasta")
  names(fasta) <- gsub("^..\\|.*\\|(.*)$", "\\1", names(fasta))
  
  lab_contams <- fasta %>% 
    .[names(.) %in% df$Accession]
  
  proteoQ::write_fasta(lab_contams, file.path(dat_dir, "Lab_Contaminants.fasta"))
}



# combine all .R files into one
# foo_combine_codes(filepath = file.path("C:/Results/R/proteoQ/inst/extdata/examples"))
foo_combine_codes <- function (filepath = file.path(file.path("~/Github/proteoQ/R"))) 
{
  filepath <- proteoM:::find_dir(filepath)
  filenames <- dir(filepath, pattern = ".R$")
  
  dir.create(file.path(filepath, "temp"))
  
  out <- purrr::map(file.path(filepath, filenames), readLines) 
  out <- purrr::reduce(out, `c`, init = NULL)
  writeLines(out, file.path(filepath, "temp/all.R"))
}


foo_list_func <- function (filepath = file.path(file.path("~/Github/proteoQ/R"))) 
{
  filepath <- proteoM:::find_dir(filepath)
  filenames <- dir(filepath, pattern = ".R$")
  
  dir.create(file.path(filepath, "temp"))
  
  lines <- purrr::map(file.path(filepath, filenames), readLines)
  
  fns_all <- lapply(lines, function (x) {
    fn_lines <- x[grepl("<-\\s*function\\s*\\(", x)]
    fns <- gsub("^(.*)\\s*<- function\\s*\\(.*", "\\1", fn_lines)
    gsub("\\s*$", "", fns)
  })
  
  names(fns_all) <- filenames
  
  sink(file.path(filepath, "temp/funs.txt"))
  fns_all
  sink()
  
  ans <- readLines(file.path(filepath, "temp/funs.txt"))
  ans <- paste0("# ", ans)
  writeLines(ans, file.path(filepath, "temp/funs.R"))
}


# zipped base protein.txt, peptide.txt and call pars
foo_zip <- function () {
  zip("my.zip", "C:/proteoQ/examples")
  unzip("./my.zip", exdir = "~/proteoQ_demo", overwrite  = FALSE)
}


# Annotation of subcellular locations
foo_subcellular <- function () {
  dat_dir <- "~/proteoQ/examples"
  
  scc <- readr::read_tsv(
    file.path(dat_dir, "uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_+AND+review--.tab")) %>% 
    dplyr::rename(uniprot_acc = Entry, scc = "Subcellular location [CC]") %>% 
    dplyr::select(uniprot_acc, scc) %>% 
    dplyr::filter(!is.na(scc))
  
  scc <- scc %>% 
    dplyr::mutate(scc = gsub(";.*\\.+?", ".", scc)) %>% 
    dplyr::mutate(scc = gsub("\\s*\\{+?.*\\}+?", "", scc)) %>% 
    dplyr::mutate(scc = gsub("\\s*Note\\=.*$", "", scc)) %>% 
    dplyr::mutate(scc = gsub("\\s*\\[+?.*\\]+?\\:", "", scc)) %>% 
    dplyr::mutate(scc = gsub("SUBCELLULAR LOCATION: ", "", scc)) %>% 
    dplyr::mutate(scc = gsub("\\.\\s*$", "", scc))
  
  scc <- as.character(scc$scc) %>% 
    strsplit("\\.") %>% 
    plyr::ldply(rbind) %>% 
    purrr::map(~ gsub("^.*\\,\\s+(.*)", "\\1", .x)) %>% 
    dplyr::bind_cols() %>% 
    `colnames<-`(paste0("scc_", names(.))) %>% 
    dplyr::bind_cols(scc, .) %>% 
    dplyr::select(-scc)
  
  run_scripts <- FALSE
  if (run_scripts) {
    scc <- scc %>% 
      tidyr::gather("id", "scc", -uniprot_acc) %>% 
      dplyr::select(-id) %>% 
      dplyr::filter(!is.na(scc))
  }

  save(scc, file = file.path(dat_dir, "scc_hs.rda"), compress = "xz")
  # load(file.path(dat_dir, "scc_hs.rda"))
}




# Prepares MaxQuant modifications 
# compile modification table
# fixed and variable mods and counts
parse_mq_mods <- function() {
  dat_dir <- "~/proteoQ/examples"
  contents <- xml2::read_xml(file.path(dat_dir, "mq_modifications.xml")) %>% 
    xml2::xml_contents() 
  
  titles <- contents %>% 
    xml2::xml_attrs() %>% 
    purrr::map_chr(`[`, "title")
  
  children <- purrr::map(contents, xml2::xml_children)
  
  pos_i <- children %>% 
    purrr::map(xml2::xml_name) %>% 
    purrr::map_dbl(~ which(.x == "position"))
  
  positions <- children %>% 
    purrr::map(xml2::xml_text) %>% 
    purrr::map2_chr(pos_i, `[`)
  
  ## --- compositions and masses ---
  compositions <- contents %>% 
    xml2::xml_attrs() %>% 
    purrr::map_chr(`[`, "composition")
  
  element_universe <- compositions %>% 
    gsub("\\([^\\(]*\\)+?", "", .) %>% 
    purrr::map(~ stringr::str_split(.x, " ", simplify = TRUE)) %>% 
    purrr::reduce(`c`) %>% 
    unique() %>% 
    .[. != ""]
  
  # https://www.sisweb.com/referenc/source/exactmas.htm
  masses_lookup <- c(
    "H" = 1.007825,
    "Hx" = 2.014102,
    "C" = 12.000000,
    "Cx" = 13.003355, 
    "N" = 14.003074,
    "Nx" = 15.000109,
    "O" = 15.994915,
    "Ox" = 17.999159,
    "P" = 30.973763,
    "S" = 31.972072, 
    "Na" = 22.989770
  )
  
  counts <- suppressWarnings(
    compositions %>% 
      stringr::str_split(" ") %>% 
      purrr::map(~ gsub("^[A-z]+\\(*(-{0,1}\\d*)\\)*$", "\\1", .x) %>% 
                   as.numeric(.) %>% 
                   replace(is.na(.), 1)) 
  )
  
  elements <- compositions %>% 
    stringr::str_split(" ") %>% 
    purrr::map(~ gsub("\\(*-{0,1}\\d*\\)*", "", .x))
  
  masses <- rep(NA, length(elements))
  for (i in seq_along(elements)) {
    masses[i] <- purrr::map2(masses_lookup[elements[[i]]], counts[[i]], `*`) %>% 
      purrr::reduce(sum)
  }
  rm(i)
  
  ## --- outputs --- 
  mq_mods <- data.frame(title = titles, position = positions, 
                        composition = compositions, mass = masses) %>% 
    dplyr::mutate(position = gsub("anyCterm", "C-term", position), 
                  position = gsub("anyNterm", "N-term", position),
                  position = gsub("proteinNterm", "Protein N-term", position), 
                  position = gsub("proteinCterm", "Protein C-term", position), 
                  position = gsub("notCterm", "Not C-term", position), 
                  position = gsub("notNterm", "Not N-term", position), 
                  position = gsub("anywhere", "Anywhere", position), )
  
  save(mq_mods, file = file.path(dat_dir, "mq_mods.rda"), compress = "xz")
}


foo_fc_tables <- function () {
  library(magrittr)
  library(dplyr)
  library(purrr)
  library(readr)
  library(openxlsx)
  
  scale_log2r <- TRUE
  NorZ_ratios <- paste0(ifelse(scale_log2r, "Z", "N"), "_log2_R")
  
  df <- readr::read_tsv(file.path(dat_dir, "Protein/Protein.txt")) %>% 
    dplyr::select("gene", grep(NorZ_ratios, names(.))) %>% 
    `names<-`(gsub("^.*[NZ]{1}_log2_R[0-9]{3}[NC]{0,1}\\s\\((.*)\\)$", "\\1", names(.))) %T>% 
    write_tsv(file.path(dat_dir, "Protein/protein_logFC.txt")) 
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "logFC")
  openxlsx::writeData(wb, sheet = "logFC", df)
  openxlsx::saveWorkbook(wb, file.path(dat_dir, "Protein/protein_logFC.xlsx"), overwrite = TRUE)
  
  df <- df %>% 
    dplyr::mutate_if(is.numeric, ~ ifelse(.x > 0, 2^.x, -1/(2^.x)) %>% round(digits = 2)) %T>% 
    write_tsv(file.path(dat_dir, "Protein/protein_FC.txt")) 
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheetName = "FC")
  openxlsx::writeData(wb, sheet = "FC", df)
  openxlsx::saveWorkbook(wb, file.path(dat_dir, "Protein/protein_FC.xlsx"), overwrite = TRUE)
  
  # --- to nmb 
  # (not run)
  df <- readr::read_tsv(file.path(dat_dir, "Protein/Protein.txt")) %>% 
    dplyr::select("gene", grep(NorZ_ratios, names(.))) %>% 
    `names<-`(gsub("^.*[NZ]{1}_log2_R[0-9]{3}[NC]{0,1}\\s\\((.*)\\)$", "\\1", names(.))) %>% 
    dplyr::mutate_if(is.numeric, ~ .x - .data$NBM) %>% 
    dplyr::select(-grep("^Ref\\.[0-9]+", names(.))) %T>% 
    write_tsv(file.path(dat_dir, "Protein/protein_logFC_to_nbm.txt")) %>% 
    dplyr::mutate_if(is.numeric, ~ ifelse(.x > 0, 2^.x, -1/(2^.x)) %>% round(digits = 2)) %T>% 
    write_tsv(file.path(dat_dir, "Protein/protein_FC_to_nbm.txt")) 
}


foo_make_fasta_peps <- function () {
  run_scripts <- FALSE
  if (run_scripts) {
    peptide <- paste0("-WASQVSENR@PVCK@AIIQGK@QFEGLVDTGADVSIIALNQWPK@NWPK@QK@", 
                      "AVTGLVGIGTASEVYQSTEILHCLGPDNQESTVQPMITSIPLNLWGR@", 
                      "DLLQQWGAEITMPAPLYSPTSQK@IMTK@MGYIPGK@GLGK@NEDGIK@", 
                      "IPFEAK@INQK@R@EGIGYPF-")
    
    gsub(paste0("^(", 
                rep("[^@]*?@{1}", max_miss+1) %>% purrr::reduce(paste0), ")", 
                "([^@]*?@{1}).*$"), 
         "\\1", peptide)
  }
  
}
