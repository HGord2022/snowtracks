# === Load required libraries ===
library(tidyverse)
library(shadowtext)

setwd("~/github_repos/snowtracks/popgen/cow/pca")

# === File paths ===
evec_file <- "cow_proj.evec"
eigval_file <- "cow_proj.eigenvalues"
pop_file    <- "ena_metadata.tsv"

# === Control label size ===
label_text_size <- 4

# === Read eigenvectors (.evec) ===
evec <- read_table2(evec_file, comment = "#", col_names = FALSE)
evec <- evec[-1,]

num_pcs <- ncol(evec) - 2
pc_names <- paste0("PC", 1:num_pcs)
colnames(evec) <- c("sample", pc_names, "sample_dup")
evec <- evec %>% select(-sample_dup)

# === Read eigenvalues ===
eigvals <- read_table2(eigval_file, col_names = FALSE) %>% pull(X1)
var_explained <- eigvals / sum(eigvals) * 100

# === Read population info ===
pop <- read_tsv(pop_file, col_names = FALSE, trim_ws = TRUE)
colnames(pop)[1:2] <- c("sample", "population")

# === Merge data ===
merged <- left_join(evec, pop, by = "sample")

# === Highlighted samples ===
highlight_colors <- c(CX1138 = "#E41A1C")
highlight_labels <- c(CX1138 = "CX1138 - Slovenia, 0.45x")

merged <- merged %>%
  mutate(highlight = sample %in% names(highlight_colors))

# === Build population colour palette (override only one) ===
pop_levels <- sort(unique(merged$population))

pop_palette <- setNames(
  scales::hue_pal()(length(pop_levels)),
  pop_levels
)

pop_palette["Holstein x Freisian"] <- "#F2A900"  # yellow–orange

# === Populations to label ===
pop_to_label <- c(
  "Simmental", "Brown Swiss", "Holstein", "Zebu",
  "Tuxer", "Tyrolean Grey", "Scottish Highland",
  "Hereford", "Gelbvieh"
)

# === Population centroids ===
pop_labels <- merged %>%
  group_by(population) %>%
  summarise(across(starts_with("PC"), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    label = population,
    across(
      starts_with("PC"),
      ~ if_else(population == "Simmental" & cur_column() == "PC2", .x, .x + 0.005)
    )
  ) %>%
  filter(population %in% pop_to_label)

# === Theme ===
base_theme <- theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# === PCA plotting function ===
plot_pca <- function(merged, var_explained, pc_x, pc_y, file_name) {
  
  p <- ggplot() +
    geom_point(
      data = merged %>% filter(!highlight),
      aes(x = .data[[pc_x]], y = .data[[pc_y]], color = population),
      size = 3,
      show.legend = FALSE
    ) +
    geom_point(
      data = merged %>% filter(highlight),
      aes(x = .data[[pc_x]], y = .data[[pc_y]], fill = sample),
      shape = 24, color = "black", size = 4, stroke = 1.2
    ) +
    scale_fill_manual(
      name = "Highlighted samples",
      values = highlight_colors,
      labels = highlight_labels,
      guide = guide_legend(override.aes = list(shape = 24, size = 5))
    ) +
    scale_color_manual(values = pop_palette) +
    geom_shadowtext(
      data = pop_labels,
      aes(x = .data[[pc_x]], y = .data[[pc_y]], label = label, colour = population),
      bg.colour = "black",
      bg.r = 0.15,
      fontface = "bold",
      size = label_text_size,
      hjust = -0.1,
      show.legend = FALSE
    ) +
    labs(
      x = paste0(pc_x, " (", round(var_explained[as.numeric(sub("PC", "", pc_x))], 1), "%)"),
      y = paste0(pc_y, " (", round(var_explained[as.numeric(sub("PC", "", pc_y))], 1), "%)")
    ) +
    base_theme
  
  ggsave(file_name, p, width = 11, height = 5, dpi = 300)
}

# === PC pairs ===
pc_pairs <- list(
  c("PC1", "PC2"),
  c("PC2", "PC3"),
  c("PC1", "PC4"),
  c("PC1", "PC5")
)

# === Plot ===
for (pair in pc_pairs) {
  plot_pca(
    merged,
    var_explained,
    pair[1],
    pair[2],
    paste0("cow_pca_", pair[1], "_", pair[2], ".png")
  )
}
