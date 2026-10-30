\donttest{
### 1. Prepare metadata to workspace
library(proteoQ)

id <- "gene"
species <- "human"
species_lookup <- c(human = "hs", mouse = "mm", rat = "rn")
abbr_species <- species_lookup[[species]]

fasta <- c("fasta/file_1.fasta","fasta/crap.fasta")
dat_dir <- "A/working/directory"

# Metadata file: expt_smry.txt; 
# auto-generated frac_smry.txt if no peptide pre-fractionation
load_expts(
  expt_smry = "expt_smry.txt",
  frac_smry = "frac_smry.txt"
)


### 2. Prepare peptide and protein tables
## 2.1 Process DIANN inputs
makePepDIANN(fasta = fasta,use_lowercase_aa = TRUE,)
makeProtDIANN(fasta = fasta)

## 2.2 Standardize peptide and protein data by human subset
standPep(
  method_align = "MC", 
  range_log2r = c(5, 95), 
  range_int = c(5, 95), 
  n_comp = 3, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
  slice_ = exprs(species == !!species),
)

standPrn(
  method_align = "MC", 
  range_log2r = c(5, 95), 
  range_int = c(5, 95), 
  n_comp = 2, 
  seed = 749662, 
  maxit = 200, 
  epsilon = 1e-05, 
  slice_ = exprs(species == !!species),
)


### 3. Classification of protein subcellular locations
## 3.1 Set up column keys in metadata
load(file.path(dat_dir, "label_scheme.rda"))

col_select <- "Select"
col_order  <- "Order"

col_fraction <- "Fraction"
levels_fraction <- c("CP", "NP", "ChA")

col_subtype <- "Type"
levels_subtype <- c("Veh", "TG")

# Construct group levels: "[levels_subtype]_[levels_fraction]"
col_group  <- "Group"
levels_group <- c(outer(levels_subtype, levels_fraction, paste, sep = "_"))
# stopifnot(all(levels_group %in% label_scheme[[col_group]]))

## 3.2 Load reference markers
# Cellular compartments coincide with experimental fractions (levels_fraction)
col_ref <- "location_slim3"
levels_pred <- levels_fraction

df_ref <- load_ref_markers(
  species = species, col_ref = col_ref, 
  levels_ref_compartments = levels_pred)

## 3.3 Unsupervised classification
n_clust <- 16:36

classify_subcellular_unsupervised(
  species = species, id = id, n_clust = n_clust, 
  dat_dir = dat_dir, label_scheme = label_scheme, 
  col_group = col_group, levels_group = levels_group, 
  col_fraction = col_fraction, levels_fraction = levels_fraction,
  col_subtype = col_subtype, levels_subtype = levels_subtype,
  col_select = col_select, col_order = col_order)

## 3.4 UMAP visualization of the whole proteome
colors_frac <- c(CP = "#E41A1C", NP = "#377EB8", ChA = "#4DAF4A")

colors_ref <- c(
  Cytosol = "#E41A1C", Nucleus = "#984EA3", Chromatin = "#FF7F00", 
  Centrosome = "#00BFFF", ER = "#377EB8", Mitochondrion = "#4DAF4A", 
  Lysosome = "#FFFF33", Golgi = "#008080", PM = "#A65628", 
  Proteasome = "#F781BF", Peroxisome = "#3A1188", Endosome = "#808000", 
  Ribosome = "#882255",P_body = "#99d8c9", Stress_granule = "#fee6ce")

levels_ref_compartments <- names(colors_ref)

stopifnot(all(levels_ref_compartments %in% unique(df_ref[["location"]])))

# Map markers to the whole proteome
method_supervision <- "unsupervised"

plot_prnSubcell_UMAP(
  species = species, 
  
  # In this case, sample fractions match cellular compartments
  method_supervision = method_supervision, 
  
  # Use "log2FC" for consistency with unsupervised learning in 
  # 'classify_subcellular_unsupervised'
  key_col = "log2FC",
  
  # Map markers to the full proteome 
  highlight_markers = TRUE, 
  df_ref = df_ref, col_ref = "location", colors_ref = colors_ref, 
  
  # Define the plot order of reference markers
  levels_ref_compartments = levels_ref_compartments,
  
  # 'col_fraction' encodes experimental fractions at type_fraction = "matched"
  # (it encodes pseudo fractions at type_fraction = "pseudo")
  col_fraction = col_fraction, levels_fraction = levels_fraction, 
  # colors_frac = colors_frac, # Not used at highlight_markers = TRUE
  
  col_group = col_group, levels_group = levels_group, 
  
  col_subtype = col_subtype, levels_subtype = levels_subtype, 
  
  # Define the order of subcellular fractions; 
  # each facet correspond to a pseudo fraction
  wrap_facets = TRUE, wrap_facets_by = col_fraction, 
  levels_pred = levels_fraction, 
  
  # Control the aesthetics
  xmin = -12, xmax = 12, ymin = -12, ymax = 12, 
  nrow = 2, ncol = 3, alpha = 1, 
  break_width = 1, hide_grids = FALSE, drop = TRUE, 
  # width_combined = 4, height_combined = 4, 
  width_sep = 12, height_sep = 4, 
  ext = "png", 
  
  # Use cached UMAP coordinates when available
  use_cache = TRUE,
  seed = 123,
)

# Color-code by sample fractions
if (FALSE) {
  plot_prnSubcell_UMAP(
    species = species, 
    type_fraction = "matched", 
    key_col = "log2FC",
    
    # Map markers to the full proteome 
    highlight_markers = FALSE, 
    # levels_ref_compartments = levels_ref_compartments,
    
    col_fraction = col_fraction, levels_fraction = levels_fraction, 
    colors_frac = colors_frac, 
    
    col_group = col_group, levels_group = levels_group, 
    col_subtype = col_subtype, levels_subtype = levels_subtype, 
    
    wrap_facets = TRUE, wrap_facets_by = col_fraction, 
    levels_pred = levels_fraction, 
    
    xmin = -12, xmax = 12, ymin = -12, ymax = 12, 
    nrow = 2, ncol = 3, alpha = 1, 
    break_width = 5, hide_grids = TRUE, drop = TRUE, 
    width_combined = 8, height_combined = 4, 
    width_sep = 10, height_sep = 6, 
    ext = "png", 
    
    use_cache = TRUE,
    seed = 123,
  )
}

## Plot TG / Veh only
if (FALSE) {
  for (i in seq_along(levels_subtype)) {
    sub_type <- levels_subtype[[i]]
    
    plot_prnSubcell_UMAP(
      species = species, 
      method_supervision = method_supervision, 
      key_col = "log2FC",
      
      # Map markers to the full proteome 
      highlight_markers = FALSE, 
      # levels_ref_compartments = levels_ref_compartments,
      
      col_fraction = col_fraction, levels_fraction = levels_fraction, 
      colors_frac = colors_frac, 
      
      col_group = col_group, levels_group = levels_group, 
      col_subtype = col_subtype, levels_subtype = sub_type, 
      
      wrap_facets = FALSE, wrap_facets_by = col_fraction, 
      levels_pred = levels_fraction, 
      
      xmin = -12, xmax = 12, ymin = -12, ymax = 12, 
      nrow = 2, ncol = 3, alpha = .5, 
      break_width = 5, hide_grids = TRUE, drop = TRUE, 
      width_combined = 3, height_combined = 2, 
      width_sep = 10, height_sep = 3, 
      filename = paste0("umap_", sub_type, ".png"), 
      ext = "png", 
      
      use_cache = TRUE,
      seed = 123,
    )
  }
}


## 4.2 User-defined Polygons of tetherome according to the above marker UMAPs
polygons <- list(
  Veh = list(
    # CP = list(x = c(-5, 0, 0, -5), y = c(-5, -5, -2, -2)), 
    CP = NULL, 
    # NP = list(x = c(-10, 8, 8, -10), y = c(0, 0, 12, 12)), 
    NP = list(x = c(-8, -4.5, 4, -.5), y = c(4, 2, 11.5, 11.5)), 
    ChA = list(x = c(-6.5, -4.5, -1, -1, -4), y = c(-.5, -4, -2, 2, 2))
    ),
  TG = list(
    CP = NULL, 
    # CP = list(x = c(-7.5, 2.5, 2.5, -7.5), y = c(-7.5, -7.5, -5, 0)), 
    # NP = list(x = c(-10, 7.5, 7.5, -10), y = c(.28, 1.8, 10, 10)), 
    NP = list(x = c(-5, 7, -1, -5), y = c(-6.5, -5, -.5, -3)), 
    ChA = list(x = c(-10, -5, -5, -6.5, -10), y = c(-1, -1, 4, 5, 1))
    )
)

## 4.3 Confirmatory visualization of polygons
# Tetherome identifications (marker_tetherome.tsv)
plot_prnSubcell_UMAP(
  species = species, 
  method_supervision = method_supervision, 
  key_col = "log2FC",
  highlight_markers = TRUE, 
  df_ref = df_ref, col_ref = "location", colors_ref = colors_ref, 
  
  # Draw polygons
  polygons = polygons, 
  
  levels_ref_compartments = levels_ref_compartments,
  col_fraction = col_fraction, levels_fraction = levels_fraction, 
  col_group = col_group, levels_group = levels_group, 
  col_subtype = col_subtype, levels_subtype = levels_subtype, 
  wrap_facets = TRUE, wrap_facets_by = col_fraction, 
  levels_pred = levels_fraction, 
  
  xmin = -12, xmax = 12, ymin = -12, ymax = 12, 
  nrow = 2, ncol = 3, alpha = 1, 
  break_width = 1, hide_grids = FALSE, drop = TRUE, 
  width_combined = 8, height_combined = 4, 
  width_sep = 12, height_sep = 4, 
  # ext = "png", 
  
  use_cache = TRUE,
  seed = 123,
)

## 4.4 Update purity, entropy and localization scores
#  (incorporate the new knowledge of subcellular specificity; 
#   tetherome information not used for the updates)
update_subcell_metrics(levels_subtype = levels_subtype)

## 4.5 Identify proteins with subcellular translocation
#  (a pair of conditions, e.g., Veh versus TG at a time)
#  (tetherome information not used)
find_from_to_transloc(
  from_to = list(c("CP", "NP"), c("NP", "ChA")), 
  levels_subtype = levels_subtype[1:2], 
  width = 9, height = 3.4)

## 5. Assess statistical p-values at peptide and protein levels
pepSig(
  frac = ~ Term["TG_CP-Veh_CP", "TG_NP-Veh_NP", "TG_ChA-Veh_ChA"], 
  Xfrac = ~ Term["(TG_NP-Veh_NP)-(TG_CP-Veh_CP)", 
                 "(TG_ChA-Veh_ChA)-(TG_NP-Veh_NP)", 
                 "(TG_ChA-Veh_ChA)-(TG_CP-Veh_CP)"],
)

prnSig()


## 6. Update protein and peptide statistics
#  (by nullifying p-values and log2Ratios for location unspecific entries)

# Choose one of the two
update_subcell_stat(rm_tetherome = TRUE)
# update_subcell_stat(rm_tetherome = FALSE)


### 7. Optional: visualize the updated results via volcano plot
prnVol(
  filter_ = exprs(prot_n_pep >= 2L, species == !!species, !grepl("^Krt", gene)
  ), 
  topn_labels = 6, 
  add_id_table = FALSE, 
  adjP = TRUE, 
  xco = 1.5, yco = 0.05,
  width = 11, height = 4.25,
  xmin = -6, xmax = 6, 
  ymin = 0, ymax = 16, 
)

pepVol(
  filter_ = exprs(species == !!species, !grepl("^Krt", gene)), 
  adjP = TRUE, 
  xco = 1.5, yco = 0.05,
  width = 16, height = 6,
)

## 8. Optional: GSPA pathway analysis
prnGSPA(
  gset_nms = "go_sets", 
  pval_cutoff = 1E-2, # protein pVal threshold
  logFC_cutoff = log2(1.5), # protein log2FC threshold
  use_adjP = FALSE,
  gslogFC_cutoff = log2(1.2),
  is_subcellular = TRUE, # immaterial given 'update_subcell' being executed
)

# volcano plot visualization of GSPA pathways
prnGSPAMap(
  show_passed_only = FALSE, 
  gspval_cutoff = 1E-2, 
  show_sig = "none", 
  topn_labels = 6, 
  add_id_table = FALSE, 
  adjP = FALSE, 
  xco = 1.5, yco = 0.01,
  width = 11, height = 4.5,
  xmin = -6, xmax = 6, 
  ymin = 0, ymax = 16, 
  filter_ = exprs(prot_n_pep >= 2L, species == !!species, !grepl("^KRT", gene)), 
  is_subcellular = TRUE, 
)

}

