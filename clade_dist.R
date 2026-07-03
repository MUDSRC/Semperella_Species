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

# ---- 1) Input files
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

pairwise_annotated <- pairwise_all %>%
  left_join(metadata, by = c("taxon1" = "Taxon")) %>%
  rename(clade1 = Clade, subclade1 = Subclade) %>%
  left_join(metadata, by = c("taxon2" = "Taxon")) %>%
  rename(clade2 = Clade, subclade2 = Subclade)

# Check unmatched taxa
unmatched_taxa <- pairwise_annotated %>%
  filter(is.na(clade1) | is.na(clade2)) %>%
  select(marker, metric, taxon1, taxon2, clade1, clade2) %>%
  distinct()

# Keep only ingroup comparisons with metadata
pairwise_ingroup <- pairwise_annotated %>%
  filter(!is.na(clade1), !is.na(clade2))

# 1. BETWEEN CLADES
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

# 2. WITHIN CLADES
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

# 3. BETWEEN SUBCLADES
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

# 4. WITHIN SUBCLADES
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

# ---- X) Save output
write_xlsx(
  list(
    pairwise_annotated = pairwise_annotated,
    between_clades = between_clades,
    within_clades = within_clades,
    between_subclades = between_subclades,
    within_subclades = within_subclades,
    unmatched_taxa = unmatched_taxa
  ),
  "clade_distance_summary.xlsx"
)
