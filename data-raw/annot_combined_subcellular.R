library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(babelgene)
library(ggplot2)
library(proteoQ)

# data("subcell_hs")
# table(subcell_hs$location)

db_path <- "D:/QZ/dbs/cc"
species <- "human"

## 1.1. Load data
df_mcp <- readr::read_tsv(file.path(db_path, "subcell_markers.txt")) |> # 1698
  dplyr::rename(prot_acc = "Protein Group", 
                uniprot_id = "Protein Name", 
                prot_des = "Protein Description", 
                gene = "Genes", 
                location = "Markers") |>
  dplyr::filter(location != "unknown") |>
  dplyr::select(dplyr::one_of(c("gene", "location")))

df_cell <- readr::read_tsv(file.path(db_path, "mmc2.txt")) |> # 2412
  dplyr::rename(gene = gene_name_canonical, location = organelle) |>
  dplyr::select(dplyr::one_of(c("gene", "location")))

## 1.2. Define translation dictionary
map_mcp <- c(
  "Chromatin"     = "Chromatin",
  "Cytosol"       = "Cytosol",
  "Endosome"      = "Endosome",
  "ER"            = "ER",
  "Golgi"         = "Golgi",
  "Lysosome"      = "Lysosome",
  "Mitochondrion" = "Mitochondrion",
  "Nucleus"       = "Nucleus",
  "Peroxisome"    = "Peroxisome",
  "PM"            = "Plasma_membrane",
  "Proteasome"    = "Proteasome",
  "Ribosome 40S"  = "Ribosome",
  "Ribosome 60S"  = "Ribosome"
)

map_cell <- c(
  "actin-binding_protein" = "Cytosol",
  "centrosome"            = "Centrosome",
  "cytosol"               = "Cytosol",
  "early_endosome"        = "Endosome",
  "ER"                    = "ER",
  "ERGIC"                 = "ER",
  "Golgi"                 = "Golgi",
  "lysosome"              = "Lysosome",
  "mitochondria"          = "Mitochondrion",
  "nuclear_pore"          = "Nucleus",
  "nucleus"               = "Nucleus",
  "p-body"                = "P_body",
  "peroxisome"            = "Peroxisome",
  "plasma membrane"       = "Plasma_membrane",
  "ribosome"              = "Ribosome",
  "stress granule"        = "Stress_granule"
)

## 1.3. Apply recoding to dataframes
df_mcp$location <- dplyr::recode(df_mcp$location, !!!map_mcp)
df_cell$location <- dplyr::recode(df_cell$location, !!!map_cell)

## 1.4. Process chromatin proteins
gn_chrom <- df_mcp[with(df_mcp, location == "Chromatin"), ] |> dplyr::pull(gene)
df_cell <- df_cell |> dplyr::filter(!gene %in% gn_chrom) # 2382
rm(list = "gn_chrom")

## 1.5. Harmonize location labels, preserving Chromatin
df_mcp$dataset <- "MCP"
df_cell$dataset <- "Cell"

df_combined <- bind_rows(df_mcp, df_cell)
table(df_combined$location, df_combined$dataset)

## 2.1. Identify all multi-dipped genes and their location assignments
multi_dipped <- df_combined |>
  dplyr::group_by(gene) |>
  dplyr::filter(dplyr::n() > 1) |>
  dplyr::summarize(
    n_occurrences = dplyr::n(),
    datasets      = paste(unique(dataset), collapse = ", "),
    locations     = paste(unique(location), collapse = "; "),
    n_unique_locs = dplyr::n_distinct(location),
    .groups = "drop"
  ) |>
  dplyr::arrange(desc(n_unique_locs), desc(n_occurrences))

## 2.2. Remove confounding entreis
bads <- multi_dipped |> dplyr::filter(n_unique_locs > 1L) |> dplyr::pull(gene)

subcell_hs <- df_combined |> # 3418
  dplyr::filter(!gene %in% bads) |>
  dplyr::select(dplyr::one_of(c("gene", "location"))) |>
  unique()

## 2.3. Combine into CP, NP and ChA
fraction_map <- c(
  "Chromatin"       = "ChA",
  "Centrosome"      = "NP",
  "Nucleus"         = "NP",
  "Cytosol"         = "CP",
  "Endosome"        = "CP",
  "ER"              = "CP",
  "Golgi"           = "CP",
  "Lysosome"        = "CP",
  "Mitochondrion"   = "CP",
  "P_body"          = "CP",
  "Peroxisome"      = "CP",
  "Plasma_membrane" = "CP",
  "Proteasome"      = "CP",
  "Ribosome"        = "CP",
  "Stress_granule"  = "CP"
)

subcell_hs <- subcell_hs |>
  dplyr::mutate(
    location_slim3 = dplyr::recode(location, !!!fraction_map),
    location_slim3 = factor(location_slim3, levels = c("CP", "NP", "ChA"))
  )

## 2.4 Merge Centrosome with Nucleus (not enough resolution)
if (FALSE) {
  subcell_hs <- subcell_hs |>
    dplyr::mutate(
      location = ifelse(location == "Centrosome", "Nucleus", location),
    )
}

readr::write_tsv(subcell_hs, file.path(db_path, "subcell_location_references_hs.tsv"))
# stopifnot(nrow(subcell_hs) == length(unique(subcell_hs$gene)))
save(subcell_hs, file = file.path(db_path, "subcell_hs.rda"))

## 3.1. Convert to mouse genes
hsmm <- proteoQ:::convert_genes_human_to_mouse(df_combined$gene) |>
  dplyr::select(dplyr::one_of(c("human_symbol", "symbol")))

subcell_mm <- subcell_hs |>
  dplyr::left_join(hsmm, by = c("gene" = "human_symbol")) |>
  dplyr::mutate(gene = symbol) |>
  dplyr::select(dplyr::one_of("gene", "location", "location_slim3")) |>
  dplyr::filter(!is.na(gene)) |>
  unique() %T>% 
  readr::write_tsv(file.path(db_path, "subcell_location_references_mm.tsv"))
# stopifnot(nrow(subcell_mm) == length(unique(subcell_mm$gene)))
save(subcell_mm, file = file.path(db_path, "subcell_mm.rda"))
