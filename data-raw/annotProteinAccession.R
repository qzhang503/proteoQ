annot_subcellular <- function (db_path = "~/proteoQ/dbs/cc", 
                               species = "human")
{
  #  | Broad compartment     | GO ID      | GO term               |
  #  | --------------------- | ---------- | --------------------- |
  #  | Nucleus               | GO:0005634 | nucleus               |
  #  | Chromatin             | GO:0000785 | chromatin             |
  #  | Cytoplasm             | GO:0005737 | cytoplasm             |
  #  | Cytosol               | GO:0005829 | cytosol               |
  #  | Plasma membrane       | GO:0005886 | plasma membrane       |
  #  | Mitochondrion         | GO:0005739 | mitochondrion         |
  #  | Endoplasmic reticulum | GO:0005783 | endoplasmic reticulum |
  #  | Golgi apparatus       | GO:0005794 | Golgi apparatus       |
  #  | Lysosome              | GO:0005764 | lysosome              |
  #  | Endosome              | GO:0005768 | endosome              |
  #  | Ribosome              | GO:0005840 | ribosome              |
  #  | Proteasome            | GO:0000502 | proteasome complex    |
  #  | Cytoskeleton          | GO:0005856 | cytoskeleton          |
  #  | Centrosome            | GO:0005813 | centrosome            |
  #  | Extracellular region  | GO:0005576 | extracellular region  |
  #  | Vesicle               | GO:0031982 | vesicle               |
  #  | Peroxisome            | GO:0005777 | peroxisome            |
  
  ### Only run once: Compilation details
  library(dplyr)
  library(stringr)
  library(readr)
  library(tidyr)
  
  db_path <- "D:/QZ/dbs/cc"
  # species <- "human"
  species <- "mouse"
  
  abbr_species <- switch(
    species, 
    human = "hs",
    mouse = "mm",
    rat = "rn",
    stop("`species` not in one of `human`, `mouse` or `rat`.")
  )
  
  levels_subcellular <- c("CP", "NP", "ChA")
  
  file_hpa <- "hpa_subcellular_location.tsv"
  file_up  <- paste0("up_", abbr_species, "_subcellular_location.tsv")
  file_go  <- paste0("tbl_go_cc_", species, ".tsv")
  
  df_up <- readr::read_tsv(file.path(db_path, file_up)) |>
    dplyr::rename(CC = "Subcellular location [CC]", 
                  gene = "Gene Names (primary)", 
                  prot_acc = "Entry", 
                  prot_desc = "Protein names")
  
  df_go <- readr::read_tsv(file.path(db_path, file_go)) |>
    tidyr::separate(
      go_name,
      into = c("go_id", "go_term"),
      sep = " ",
      extra = "merge"
    ) |>
    dplyr::rename(location = go_term) |>
    dplyr::mutate(location = tolower(location)) |>
    dplyr::filter(gene %in% df_up$gene)
  
  if (species == "human") {
    df_hpa <- readr::read_tsv(file.path(db_path, file_hpa)) |>
      dplyr::rename(gene = "Gene name", 
                    reliability = "Reliability", 
                    location = "Main location") |>
      dplyr::filter(gene %in% df_up$gene)
  } else {
    df_hpa <- NULL
  }
  
  ### (1) Human protein atlas
  if (species == "human") {
    df_hpa <- df_hpa |>
      dplyr::filter(reliability %in% c("Enhanced", "Approved", "Supported")) |>
      dplyr::select(c("gene", "reliability", "location"))
    
    df_hpa <- df_hpa |>
      tidyr::separate_rows(location, sep = ";") |>
      dplyr::mutate(location = trimws(location)) |> 
      dplyr::filter(!is.na(location))
    
    if (TRUE) {
      df_hpa <- df_hpa |>
        dplyr::mutate(
          broad_compartment = dplyr::case_when(
            location %in% c(
              "Nucleoplasm",
              "Nuclear speckles",
              "Nucleoli fibrillar center",
              "Nuclear membrane",
              "Nuclear bodies",
              "Nucleoli",
              "Nucleoli rim"
            ) ~ "Nucleus",
            location %in% c(
              "Mitotic chromosome",
              "Kinetochore"
            ) ~ "Chromatin",
            location %in% c(
              "Cytoplasmic bodies",
              "Aggresome",
              "Rods & Rings"
            ) ~ "Cytoplasm",
            location == "Cytosol" ~ "Cytosol",
            location %in% c(
              "Plasma membrane",
              "Cell Junctions",
              "Focal adhesion sites",
              "Cleavage furrow"
            ) ~ "Plasma membrane",
            location == "Mitochondria" ~ "Mitochondrion",
            location == "Endoplasmic reticulum" ~ "Endoplasmic reticulum",
            location == "Golgi apparatus" ~ "Golgi apparatus",
            location == "Lysosomes" ~ "Lysosome",
            location == "Endosomes" ~ "Endosome",
            FALSE & location == "Ribosome" ~ "Ribosome",
            FALSE & location == "Proteasome" ~ "Proteasome",
            location %in% c(
              "Microtubules",
              "Intermediate filaments",
              "Actin filaments",
              "Mitotic spindle",
              "Microtubule ends",
              "Cytokinetic bridge",
              "Midbody",
              "Midbody ring"
            ) ~ "Cytoskeleton",
            location %in% c(
              "Centrosome",
              "Centriolar satellite",
              "Basal body",
              "Flagellar centriole"
            ) ~ "Centrosome",
            location %in% c(
              "Acrosome",
              "Perinuclear theca",
              "Calyx",
              "Annulus",
              "Equatorial segment",
              "Connecting piece",
              "End piece",
              "Mid piece",
              "Principal piece"
            ) ~ "Extracellular region",
            location %in% c(
              "Vesicles",
              "Lipid droplets"
            ) ~ "Vesicle",
            location == "Peroxisomes" ~ "Peroxisome",
            location %in% c(
              "Primary cilium",
              "Primary cilium tip",
              "Primary cilium transition zone"
            ) ~ "Cytoskeleton",
            TRUE ~ NA_character_
          )
        )
    }
    
    df_hpa <- df_hpa |>
      dplyr::rename(location_fine = location, 
                    location = broad_compartment)
    
    df_hpa <- df_hpa |>
      dplyr::mutate(
        compartment_2 = dplyr::case_when(
          location %in% c(
            "Nucleus",
            "Chromatin",
            "Centrosome"
          ) ~ "NP",
          location %in% c(
            "Cytoplasm",
            "Cytosol",
            "Cytoskeleton",
            "Mitochondrion",
            "Endoplasmic reticulum",
            "Golgi apparatus",
            "Endosome",
            "Lysosome",
            "Peroxisome",
            "Vesicle",
            "Plasma membrane",
            "Extracellular region"
          ) ~ "CP",
          TRUE ~ NA_character_
        )
      ) |>
      dplyr::mutate(location = tolower(location))
    
    readr::write_tsv(
      df_hpa, file.path(db_path, "hpa_subcellular_location_edited.tsv"))
  }
  
  map_loc <- function (df) {
    df_go_b %>%
      dplyr::mutate(
        location = case_when(
          
          # =========================
          # Chromatin
          # =========================
          stringr::str_detect(location,
                              "chromatin|chromosome|centromere|kinetochore|telomere") ~
            "chromatin",
          
          # =========================
          # Nucleus
          # =========================
          stringr::str_detect(location,
                              "nucleus|nucleolus|nucleoplasm|nuclear|cajal body|pml body|gem") ~
            "nucleus",
          
          # =========================
          # Centrosome
          # =========================
          stringr::str_detect(location,
                              "centrosome|centriole|centriolar|microtubule organizing center|spindle pole body") ~
            "centrosome",
          
          # =========================
          # Cytoskeleton
          # =========================
          stringr::str_detect(location,
                              "cytoskeleton|microtubule|filament|lamellipodium|filopodium|podosome|invadopodium|stress fiber|ruffle|cell cortex|cleavage furrow|spindle|midbody|myofibril|sarcomere|z line|m line|i band|a band") ~
            "cytoskeleton",
          
          # =========================
          # Mitochondrion
          # =========================
          stringr::str_detect(location,
                              "mitochond") ~
            "mitochondrion",
          
          # =========================
          # ER
          # =========================
          stringr::str_detect(location,
                              "endoplasmic reticulum|sarcoplasmic reticulum|microsome") ~
            "endoplasmic reticulum",
          
          # =========================
          # Golgi
          # =========================
          stringr::str_detect(location,
                              "golgi") ~
            "golgi apparatus",
          
          # =========================
          # Lysosome
          # =========================
          stringr::str_detect(location,
                              "lysosome|autolysosome") ~
            "lysosome",
          
          # =========================
          # Endosome
          # =========================
          stringr::str_detect(location,
                              "endosome|multivesicular body") ~
            "endosome",
          
          # =========================
          # Peroxisome
          # =========================
          stringr::str_detect(location,
                              "peroxisome") ~
            "peroxisome",
          
          # =========================
          # Vesicle
          # =========================
          stringr::str_detect(location,
                              "vesicle|granule|phagosome|autophagosome|vacuole|secretory|exosome") ~
            "vesicle",
          
          # =========================
          # Extracellular
          # =========================
          stringr::str_detect(location,
                              "extracellular|secreted|basement membrane|matrix|synaptic cleft|zona pellucida") ~
            "extracellular region",
          
          # =========================
          # Plasma membrane
          # =========================
          stringr::str_detect(location,
                              "cell membrane|plasma membrane|membrane raft|cell surface|sarcolemma|tight junction|gap junction|adherens junction|desmosome|hemidesmosome|cell junction|caveola") ~
            "plasma membrane",
          
          # =========================
          # Cytosol
          # =========================
          stringr::str_detect(location,
                              "^cytosol$") ~
            "cytosol",
          
          # =========================
          # Cytoplasm
          # =========================
          stringr::str_detect(location,
                              "cytoplasm|cytoplasmic|perikaryon|sarcoplasm|p-body|stress granule") ~
            "cytoplasm",
          
          TRUE ~ NA_character_
        )
      ) |>
      dplyr::filter(!is.na(location))
  }
  
  
  ### (3) GO
  # 1896 Go terms
  go_keep <- c(
    "GO:0005634", # nucleus
    "GO:0000785", # chromatin
    "GO:0005737", # cytoplasm
    "GO:0005829", # cytosol
    "GO:0005886", # plasma membrane
    "GO:0005739", # mitochondrion
    "GO:0005783", # ER
    "GO:0005794", # Golgi
    "GO:0005764", # lysosome
    "GO:0005768", # endosome
    "GO:0005840", # ribosome
    "GO:0000502", # proteasome
    "GO:0005856", # cytoskeleton
    "GO:0005813", # centrosome
    "GO:0005576", # extracellular region
    "GO:0031982", # vesicle
    "GO:0005681", # spliceosomal complex, nucleus and CP
    "GO:0005777"  # peroxisome
  )
  
  rows <- df_go$go_id %in% go_keep
  df_go_a <- df_go[rows, ] # gene: 18805 -> 17189
  df_go_b <- df_go[!rows, ]
  df_go_a[["go_id"]] <- NULL
  df_go_b[["go_id"]] <- NULL
  # rm(list = c("rows"))
  
  # 1616 genes
  df_go_b <-local({
    gs_b <- unique(df_go_b$gene)
    gs_a <- unique(df_go_a$gene)
    gs_  <- unique(df_go$gene)
    gnx  <- gs_b[!gs_b %in% gs_a]
    
    df_go_b <- df_go_b |> dplyr::filter(gene %in% gnx)
  }) |> 
    map_loc() |>
    unique()
  
  df_go_a <- df_go_a |>
    dplyr::mutate(location = dplyr::case_when(
      location == "ribosome" ~ "cytosol;endoplasmic reticulum",
      location == "proteasome complex" ~ "cytosol;nucleus",
      TRUE ~ location
    )) |>
    tidyr::separate_rows(location, sep = ";") |>
    dplyr::select(c("gene", "location")) |>
    dplyr::distinct(gene, location, .keep_all = TRUE) |>
    unique()
  
  df_go <- dplyr::bind_rows(df_go_a, df_go_b) |>
    unique()
  
  df_go <- df_go |>
    dplyr::mutate(
      compartment_2 = dplyr::case_when(
        location %in% c("nucleus", "spliceosomal complex") ~ "NP",
        location %in% c("chromatin") ~ "ChA",
        location %in% c(
          "cytoplasm",
          "cytosol",
          "plasma membrane",
          "extracellular region",
          "mitochondrion",
          "endoplasmic reticulum",
          "golgi apparatus",
          "lysosome",
          "endosome",
          "vesicle",
          "cytoskeleton",
          "centrosome",
          "peroxisome",
          "ribosome",
          "proteasome complex"
        ) ~ "CP",
        TRUE ~ NA_character_
      )
    )
  
  # for duplicated assign assignments
  if (FALSE) {
    df_go <- df_go |>
      dplyr::filter(location == "spliceosomal complex") |>
      dplyr::mutate(compartment_2 = "CP")
  }
  
  # 18093 genes
  readr::write_tsv(
    df_go, file.path(db_path, paste0("tbl_go_cc_", species, "_edited.tsv")))
  
  
  ### (3) Uniprot
  df_up <- df_up |>
    dplyr::filter(!is.na(CC)) |>
    dplyr::mutate(CC = gsub("Note=.*", "", CC))
  
  # Replace [Isoform A.1], [Isoform 1] before replacing "."
  df_up <- df_up |>
    dplyr::mutate(
      CC = gsub("\\[[^]]+\\]:\\s*", "", CC), 
      CC = stringr::str_remove(CC, "^SUBCELLULAR LOCATION:\\s*"),
      CC = stringr::str_remove_all(CC, "\\s*\\{[^\\}]+\\}"),
      CC = stringr::str_replace_all(CC, "\\.\\s*$", "")
    ) |>
    tidyr::separate_rows(CC, sep = "\\.\\s*") |>
    dplyr::mutate(CC = tolower(CC)) |>
    dplyr::mutate(
      location = stringr::str_trim(CC)
    ) |>
    dplyr::filter(location != "", ) |>
    dplyr::select(-CC)
  
  df_up <- df_up |>
    tidyr::separate_rows(location, sep = ";\\s*") |>
    dplyr::filter(location != "") |>
    dplyr::mutate(
      location = gsub("subcellular location:\\s*", "", location), 
    ) |>
    tidyr::separate_rows(location, sep = ",\\s*") |>
    unique()
  
  # Collapse location to major compartments
  # 12706 genes
  df_up <- map_loc(df_up)
  
  df_up <- df_up |>
    dplyr::distinct(gene, location, .keep_all = TRUE)
  
  df_up <- df_up |>
    dplyr::mutate(
      compartment_2 = dplyr::case_when(
        location %in% c("nucleus") ~ "NP",
        location %in% c("chromatin") ~ "ChA",
        location %in% c(
          "cytoplasm",
          "cytosol",
          "plasma membrane",
          "extracellular region",
          "mitochondrion",
          "endoplasmic reticulum",
          "golgi apparatus",
          "lysosome",
          "endosome",
          "vesicle",
          "cytoskeleton",
          "centrosome",
          "peroxisome",
          "ribosome",
          "proteasome complex"
        ) ~ "CP",
        TRUE ~ NA_character_
      )
    )
  
  cols_keep <- colnames(df_go)
  df_up <- df_up[, cols_keep]
  
  readr::write_tsv(
    df_up, file.path(db_path, paste0("up_", species, "_subcellular_location_edited.tsv")))
  
  ## Put together
  # df_go$source <- "GO"
  # df_up$source <- "UP"
  if (species == "human") {
    df_hpa <- df_hpa[, cols_keep] |>
      unique()
    # all(df_hpa$location %in% df_go$location)
    
    df <- dplyr::bind_rows(
      df_go  |> dplyr::mutate(source = "GO"), 
      df_up  |> dplyr::mutate(source = "UP"), 
      df_hpa |> dplyr::mutate(source = "HPA"), 
    ) |> 
      unique()
  } else {
    df <- dplyr::bind_rows(
      df_go  |> dplyr::mutate(source = "GO"), 
      df_up  |> dplyr::mutate(source = "UP"), 
    ) |> 
      unique()
  }
  
  readr::write_tsv(df, file.path(db_path, paste0(abbr_species, "_up_go_3_components_union.tsv")))
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
  filepath <- mzion:::find_dir(filepath)
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
