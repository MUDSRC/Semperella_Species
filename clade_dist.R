# ~ HEADER --------------------------------------------
#
# ~ Author:         Alfredo Marchiò
# ~ Email:          alfredo.marchio@research.uwa.edu.au
# ~ Organization:   Minderoo-UWA Deep-Sea Research Centre
# 
# ~ Date:           2026-07-03
# ~ Version:        1.0
#
# ~ Script Name:    clade_dist.R
#
# ~ Script Description:
# Calculate distances between and within clades based on
# triangular distance matrices between single samples
# from MEGA12. 
#
# Copyright 2026 - Alfredo Marchiò
#
# ----------------------------------------------------

# ---- Libraries ----------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(writexl)

# ---- Parameters 
dist_file <- "distances_samples.xlsx"
meta_file <- "clade_metadata.csv"

# ---- Helpers
clean_taxon <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("_", " ") %>%
    str_squish()
}

# ---- 1) Input files -----
metadata <- read.csv(meta_file, stringsAsFactors = FALSE) %>%
  select(Taxon, Clade, Subclade) %>%
  filter(!is.na(Taxon), Taxon != "") %>%
  mutate(
    Taxon = clean_taxon(Taxon),
    Clade = str_squish(Clade),
    Subclade = str_squish(Subclade)
  )

sheets <- excel_sheets(dist_file)

pairwise_all <- map_dfr(sheets, function(sh) {
  
  marker <- str_extract(sh, "^[^_]+")
  metric <- str_remove(sh, "^[^_]+_")
  
  mat <- read_excel(dist_file, sheet = sh, col_names = FALSE)
  
  col_taxa <- clean_taxon(unlist(mat[1, -1]))
  row_taxa <- clean_taxon(unlist(mat[-1, 1]))
  
  values <- mat[-1, -1]
  colnames(values) <- col_taxa
  
  values %>%
    mutate(taxon1 = row_taxa) %>%
    pivot_longer(
      cols = -taxon1,
      names_to = "taxon2",
      values_to = "value"
    ) %>%
    filter(!is.na(value)) %>%
    mutate(
      marker = marker,
      metric = metric,
      taxon1 = clean_taxon(taxon1),
      taxon2 = clean_taxon(taxon2),
      value = as.numeric(value)
    )
})

# ---- 2) Data cleaning and joins metadata with distances -----
pairwise_annotated <- pairwise_all %>%
  left_join(metadata, by = c("taxon1" = "Taxon")) %>%
  rename(
    clade1 = Clade,
    subclade1 = Subclade
  ) %>%
  left_join(metadata, by = c("taxon2" = "Taxon")) %>%
  rename(
    clade2 = Clade,
    subclade2 = Subclade
  )

# Check unmatched taxa
unmatched_taxa <- pairwise_annotated %>%
  filter(is.na(clade1) | is.na(clade2)) %>%
  select(marker, metric, taxon1, taxon2, clade1, clade2) %>%
  distinct()

# Keep only ingroup comparisons with metadata
pairwise_ingroup <- pairwise_annotated %>%
  filter(!is.na(clade1), !is.na(clade2))

# ---- 3) Distance between clades -----
between_clades <- pairwise_ingroup %>%
  filter(clade1 != clade2) %>%
  mutate(
    clade_a = pmin(clade1, clade2),
    clade_b = pmax(clade1, clade2),
    comparison = paste(clade_a, clade_b, sep = "_vs_")
  ) %>%
  group_by(marker, metric, comparison, clade_a, clade_b) %>%
  summarise(
    n_pairs = n(),
    mean = mean(value),
    sd = ifelse(n() > 1, sd(value), NA),
    min = min(value),
    median = median(value),
    max = max(value),
    .groups = "drop"
  )

# ---- 4) Distance within clades -----
within_clades <- pairwise_ingroup %>%
  filter(clade1 == clade2) %>%
  group_by(marker, metric, clade = clade1) %>%
  summarise(
    comparison = paste0("within_", first(clade1)),
    n_pairs = n(),
    mean = mean(value),
    sd = ifelse(n() > 1, sd(value), NA),
    min = min(value),
    median = median(value),
    max = max(value),
    .groups = "drop"
  )

# ---- 5) Distance between sub-clades -----
# Defined only within the same main clade.
between_subclades <- pairwise_ingroup %>%
  filter(clade1 == clade2, subclade1 != subclade2) %>%
  mutate(
    subclade_a = pmin(subclade1, subclade2),
    subclade_b = pmax(subclade1, subclade2),
    comparison = paste(subclade_a, subclade_b, sep = "_vs_")
  ) %>%
  group_by(marker, metric, clade = clade1, comparison, subclade_a, subclade_b) %>%
  summarise(
    n_pairs = n(),
    mean = mean(value),
    sd = ifelse(n() > 1, sd(value), NA),
    min = min(value),
    median = median(value),
    max = max(value),
    .groups = "drop"
  )

# ---- 6) Distance within sub-clades -----
within_subclades <- pairwise_ingroup %>%
  filter(clade1 == clade2, subclade1 == subclade2) %>%
  group_by(marker, metric, clade = clade1, subclade = subclade1) %>%
  summarise(
    comparison = paste0("within_", first(subclade1)),
    n_pairs = n(),
    mean = mean(value),
    sd = ifelse(n() > 1, sd(value), NA),
    min = min(value),
    median = median(value),
    max = max(value),
    .groups = "drop"
  )

# ---- 7) Distances among individual Semperella and Semperella-like  -----
semperella_distances <- pairwise_ingroup %>%
  filter(
    subclade1 %in% c("Semperella", "Semperella-like"),
    subclade2 %in% c("Semperella", "Semperella-like"),
    taxon1 != taxon2
  ) %>%
  mutate(
    # Standardise the order of the two specimens
    taxon_a = pmin(taxon1, taxon2),
    taxon_b = pmax(taxon1, taxon2),
    
    subclade_a = if_else(
      taxon1 == taxon_a,
      subclade1,
      subclade2
    ),
    
    subclade_b = if_else(
      taxon2 == taxon_b,
      subclade2,
      subclade1
    ),
    
    comparison = paste(taxon_a, taxon_b, sep = "_vs_"),
    
    comparison_type = case_when(
      subclade_a == "Semperella" &
        subclade_b == "Semperella" ~
        "within Semperella",
      
      subclade_a == "Semperella-like" &
        subclade_b == "Semperella-like" ~
        "within Semperella-like",
      
      TRUE ~ "Semperella vs Semperella-like"
    )
  ) %>%
  select(
    marker,
    metric,
    comparison_type,
    comparison,
    taxon_a,
    taxon_b,
    subclade_a,
    subclade_b,
    value
  ) %>%
  arrange(
    marker,
    metric,
    comparison_type,
    taxon_a,
    taxon_b
  )

# p-distance table
semperella_pdistance_table <- semperella_distances %>%
  filter(metric == "p-distance") %>%
  mutate(
    Distance = round(value, 4)
  ) %>%
  select(
    Gene = marker,
    `Comparison type` = comparison_type,
    `Specimen 1` = taxon_a,
    `Specimen 2` = taxon_b,
    Distance
  )

# bp difference table
semperella_bp_table <- semperella_distances %>%
  filter(metric == "bp") %>%
  mutate(
    Differences = as.integer(value)
  ) %>%
  select(
    Gene = marker,
    `Comparison type` = comparison_type,
    `Specimen 1` = taxon_a,
    `Specimen 2` = taxon_b,
    Differences
  )

# ---- 8) Complete Semperella p-distance summary --------------------------
semperella_summary <- semperella_distances %>%
  filter(metric == "p-distance") %>%
  mutate(
    Gene = case_when(
      str_to_lower(marker) == "16s" ~ "16S",
      str_to_lower(marker) == "28s" ~ "28S",
      str_to_lower(marker) == "coi" ~ "COI",
      TRUE ~ marker
    )
  ) %>%
  select(
    `Comparison type` = comparison_type,
    `Specimen 1` = taxon_a,
    `Specimen 2` = taxon_b,
    Gene,
    Distance = value
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = Gene,
    values_from = Distance
  ) %>%
  arrange(
    factor(
      `Comparison type`,
      levels = c(
        "Within Semperella",
        "Within Semperella-like",
        "Between lineages"
      )
    ),
    `Specimen 1`,
    `Specimen 2`
  ) %>%
  mutate(
    across(
      any_of(c("16S", "28S", "COI")),
      ~ round(.x, 4)
    )
  )

# ---- 9) Genetic-distance gap check --------------------------------------

semperella_gap_check <- semperella_distances %>%
  filter(metric == "p-distance") %>%
  mutate(
    distance_class = if_else(
      subclade_a == subclade_b,
      "within",
      "between"
    )
  ) %>%
  group_by(marker, distance_class) %>%
  summarise(
    n_comparisons = n(),
    minimum = min(value, na.rm = TRUE),
    maximum = max(value, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = distance_class,
    values_from = c(
      n_comparisons,
      minimum,
      maximum,
      mean
    )
  ) %>%
  mutate(
    genetic_gap = minimum_between - maximum_within,
    all_within_lower = maximum_within < minimum_between,
    interpretation = if_else(
      all_within_lower,
      "All within-lineage distances are lower than all between-lineage distances",
      "Distance ranges overlap"
    )
  ) %>%
  mutate(
    across(
      c(
        minimum_within,
        maximum_within,
        mean_within,
        minimum_between,
        maximum_between,
        mean_between,
        genetic_gap
      ),
      ~ round(.x, 4)
    )
  )

# ---- 10) Save outputs -----
write_xlsx(
  list(
    "Between clades" = between_clades,
    "Within clades" = within_clades,
    "Between subclades" = between_subclades,
    "Within subclades" = within_subclades,
    "semperella_distances" = semperella_distances,
    "semperella_summary" = semperella_summary,
    "semperella_gap_check" = semperella_gap_check
  ),
  "genetic_distance_summaries.xlsx"
)
