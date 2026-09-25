# Retrieved from clade_plot_v6.R; full panel titles removed on 2026-09-24.
# ~ HEADER --------------------------------------------
#
# ~ Author:         Alfredo Marchiò
# ~ Email:          alfredo.marchio@research.uwa.edu.au
# ~ Organization:   Minderoo-UWA Deep-Sea Research Centre
#
# ~ Date:           2026-08-02
# ~ Version:        1.3
#
# ~ Script Name:    clade_plot_v6.R
#
# ~ Script Description:
# Create a three-panel boxplot figure showing all
# single uncorrected p-distance values within and
# between:
#   A) the two main Pheronematidae clades;
#   B) Poliopogon, Semperella, and Semperella-like;
#   C) Semperella versus Semperella-like.
#
# White boxplots represent within-lineage comparisons.
# Black boxplots represent between-lineage comparisons.
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
library(ggplot2)
library(patchwork)

# ---- Parameters ---------------------------------------------------------
dist_file <- "distances_samples.xlsx"
meta_file <- "clade_metadata.csv"

plot_dir <- "plots/genetic_distances"
plot_png <- file.path(plot_dir, "clade_genetic_distances_boxplot_ABC.png")
plot_pdf <- file.path(plot_dir, "clade_genetic_distances_boxplot_ABC.pdf")
plot_data_file <- file.path(plot_dir, "clade_genetic_distances_boxplot_ABC_data.csv")

marker_order <- c("16S", "28S", "COI")
target_subclades <- c("Poliopogon", "Semperella", "Semperella-like")

figure_width_mm <- 180
figure_height_mm <- 190
figure_dpi <- 600

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Helper -------------------------------------------------------------
clean_taxon <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("_", " ") %>%
    str_squish()
}

# ---- 1) Read metadata ---------------------------------------------------
metadata <- read.csv(meta_file, stringsAsFactors = FALSE) %>%
  select(Taxon, Clade, Subclade) %>%
  filter(!is.na(Taxon), Taxon != "") %>%
  mutate(
    Taxon = clean_taxon(Taxon),
    Clade = str_squish(Clade),
    Subclade = str_squish(Subclade)
  )

# ---- 2) Read pairwise matrices -----------------------------------------
sheets <- excel_sheets(dist_file)

pairwise_all <- map_dfr(sheets, function(sh) {
  
  marker_raw <- str_extract(sh, "^[^_]+")
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
      marker = case_when(
        str_to_lower(marker_raw) == "16s" ~ "16S",
        str_to_lower(marker_raw) == "28s" ~ "28S",
        str_to_lower(marker_raw) == "coi" ~ "COI",
        TRUE ~ marker_raw
      ),
      metric = metric,
      taxon1 = clean_taxon(taxon1),
      taxon2 = clean_taxon(taxon2),
      value = as.numeric(value)
    )
})

# ---- 3) Annotate pairwise data -----------------------------------------
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

pairwise_ingroup <- pairwise_annotated %>%
  filter(
    metric == "p-distance",
    !is.na(clade1),
    !is.na(clade2),
    !is.na(marker)
  )

# ---- 4) Panel A: main clades -------------------------------------------
panel_a_within <- pairwise_ingroup %>%
  filter(clade1 == clade2) %>%
  transmute(
    panel = "A",
    marker,
    distance_class = "Within lineages",
    comparison = paste("Within", clade1),
    distance_percent = value * 100
  )

panel_a_between <- pairwise_ingroup %>%
  filter(clade1 != clade2) %>%
  mutate(
    clade_a = pmin(clade1, clade2),
    clade_b = pmax(clade1, clade2)
  ) %>%
  transmute(
    panel = "A",
    marker,
    distance_class = "Between lineages",
    comparison = paste(clade_a, "vs", clade_b),
    distance_percent = value * 100
  )

panel_a_data <- bind_rows(panel_a_within, panel_a_between)

# ---- 5) Identify focal main clade --------------------------------------
subclade_records <- bind_rows(
  pairwise_ingroup %>% transmute(clade = clade1, subclade = subclade1),
  pairwise_ingroup %>% transmute(clade = clade2, subclade = subclade2)
) %>%
  filter(subclade %in% target_subclades) %>%
  distinct(clade, subclade)

target_main_clade <- subclade_records %>%
  count(clade, name = "n_target_subclades") %>%
  arrange(desc(n_target_subclades)) %>%
  slice(1) %>%
  pull(clade)

if (length(target_main_clade) == 0) {
  stop(
    "Could not identify the main clade containing ",
    paste(target_subclades, collapse = ", "),
    "."
  )
}

# ---- 6) Panel B: focal subclades ---------------------------------------
panel_b_within <- pairwise_ingroup %>%
  filter(
    clade1 == target_main_clade,
    clade2 == target_main_clade,
    subclade1 == subclade2,
    subclade1 %in% target_subclades
  ) %>%
  transmute(
    panel = "B",
    marker,
    distance_class = "Within lineages",
    comparison = paste("Within", subclade1),
    distance_percent = value * 100
  )

panel_b_between <- pairwise_ingroup %>%
  filter(
    clade1 == target_main_clade,
    clade2 == target_main_clade,
    subclade1 %in% target_subclades,
    subclade2 %in% target_subclades,
    subclade1 != subclade2
  ) %>%
  mutate(
    subclade_a = pmin(subclade1, subclade2),
    subclade_b = pmax(subclade1, subclade2)
  ) %>%
  transmute(
    panel = "B",
    marker,
    distance_class = "Between lineages",
    comparison = paste(subclade_a, "vs", subclade_b),
    distance_percent = value * 100
  )

panel_b_data <- bind_rows(panel_b_within, panel_b_between)

# ---- 7) Panel C: Semperella vs Semperella-like -------------------------
panel_c_within <- pairwise_ingroup %>%
  filter(
    clade1 == target_main_clade,
    clade2 == target_main_clade,
    subclade1 %in% c("Semperella", "Semperella-like"),
    subclade2 %in% c("Semperella", "Semperella-like"),
    subclade1 == subclade2
  ) %>%
  transmute(
    panel = "C",
    marker,
    distance_class = "Within lineages",
    comparison = case_when(
      subclade1 == "Semperella" ~ "Within Semperella",
      subclade1 == "Semperella-like" ~ "Within Semperella-like",
      TRUE ~ "Within lineages"
    ),
    distance_percent = value * 100
  )

panel_c_between <- pairwise_ingroup %>%
  filter(
    clade1 == target_main_clade,
    clade2 == target_main_clade,
    subclade1 %in% c("Semperella", "Semperella-like"),
    subclade2 %in% c("Semperella", "Semperella-like"),
    subclade1 != subclade2
  ) %>%
  transmute(
    panel = "C",
    marker,
    distance_class = "Between lineages",
    comparison = "Semperella vs Semperella-like",
    distance_percent = value * 100
  )

panel_c_data <- bind_rows(panel_c_within, panel_c_between)

# ---- 8) Combine and save plotting data ---------------------------------
plot_data <- bind_rows(panel_a_data, panel_b_data, panel_c_data) %>%
  mutate(
    marker = factor(marker, levels = marker_order),
    distance_class = factor(
      distance_class,
      levels = c("Within lineages", "Between lineages")
    )
  ) %>%
  arrange(panel, marker, distance_class, comparison)

write.csv(plot_data, plot_data_file, row.names = FALSE)

# ---- 9) Shared plot settings -------------------------------------------
fill_values <- c(
  "Within lineages" = "white",
  "Between lineages" = "black"
)

base_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    plot.margin = margin(t = 5, r = 8, b = 4, l = 5)
  )

make_panel_plot <- function(dat, panel_title, x_title = NULL, show_legend = FALSE) {
  
  box_offset <- 0.18
  box_width <- 0.30
  
  ggplot() +
    
    # Within-lineage boxplots: left side, white fill, black median
    geom_boxplot(
      data = dat %>% filter(distance_class == "Within lineages"),
      aes(
        x = as.numeric(marker) - box_offset,
        y = distance_percent,
        group = marker,
        fill = distance_class
      ),
      width = box_width,
      colour = "black",
      median.colour = "black",
      outlier.size = 1.4,
      outlier.stroke = 0.4
    ) +
    
    # Between-lineage boxplots: right side, black fill, white median
    geom_boxplot(
      data = dat %>% filter(distance_class == "Between lineages"),
      aes(
        x = as.numeric(marker) + box_offset,
        y = distance_percent,
        group = marker,
        fill = distance_class
      ),
      width = box_width,
      colour = "black",
      median.colour = "white",
      outlier.size = 1.4,
      outlier.stroke = 0.4
    ) +
    
    scale_x_continuous(
      breaks = seq_along(marker_order),
      labels = marker_order,
      limits = c(0.55, length(marker_order) + 0.45)
    ) +
    
    scale_fill_manual(
      values = fill_values,
      breaks = c("Within lineages", "Between lineages")
    ) +
    
    scale_y_continuous(
      labels = scales::label_number(accuracy = 0.1),
      expand = expansion(mult = c(0.03, 0.10))
    ) +
    
    labs(
      title = panel_title,
      x = x_title,
      y = "Uncorrected p-distance (%)"
    ) +
    
    base_theme +
    
    guides(
      fill = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    
    theme(
      legend.position = if (show_legend) "bottom" else "none"
    )
}

# ---- 10) Build panels ---------------------------------------------------
plot_a <- make_panel_plot(
  dat = panel_a_data %>% mutate(marker = factor(marker, levels = marker_order),
                                distance_class = factor(distance_class, levels = c("Within lineages", "Between lineages"))),
  panel_title = NULL,
  x_title = NULL,
  show_legend = FALSE
)

plot_b <- make_panel_plot(
  dat = panel_b_data %>% mutate(marker = factor(marker, levels = marker_order),
                                distance_class = factor(distance_class, levels = c("Within lineages", "Between lineages"))),
  panel_title = NULL,
  x_title = NULL,
  show_legend = FALSE
)

plot_c <- make_panel_plot(
  dat = panel_c_data %>% mutate(marker = factor(marker, levels = marker_order),
                                distance_class = factor(distance_class, levels = c("Within lineages", "Between lineages"))),
  panel_title = NULL,
  x_title = "Molecular marker",
  show_legend = TRUE
)

# ---- 11) Combine and export --------------------------------------------
combined_plot <- plot_a / plot_b / plot_c +
  plot_layout(heights = c(1, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )

ggsave(
  filename = plot_png,
  plot = combined_plot,
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  dpi = figure_dpi,
  bg = "white"
)

ggsave(
  filename = plot_pdf,
  plot = combined_plot,
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  device = grDevices::cairo_pdf,
  bg = "white"
)

cat(
  "\nBoxplot figure saved successfully:\n",
  "- ", plot_png, "\n",
  "- ", plot_pdf, "\n",
  "- ", plot_data_file, "\n",
  sep = ""
)
