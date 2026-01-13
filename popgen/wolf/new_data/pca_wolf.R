# === Load required libraries ===
library(tidyverse)
library(shadowtext)

setwd("~/github_repos/snowtracks/popgen/wolf/new_data")

# === File paths ===
evec_file <- "wolf_proj.evec"
eigval_file <- "wolf_proj.eigenvalues"
pop_file  <- "populations.txt"

# === Control label size ===
label_text_size <- 4  # adjust as needed

# === Read eigenvectors (.evec) ===
evec <- read_table2(
  evec_file,
  comment = "#",
  col_names = FALSE
)
evec <- evec[-1,]

num_pcs <- ncol(evec) - 2
pc_names <- paste0("PC", 1:num_pcs)
colnames(evec) <- c("sample", pc_names, "sample_dup")
evec <- evec %>% select(-sample_dup)

# === Read eigenvalues (.eigenvalues) ===
eigvals <- read_table2(eigval_file, col_names = FALSE) %>% pull(X1)
var_explained <- eigvals / sum(eigvals) * 100

# === Read population info ===
pop <- read_tsv(pop_file, col_names = FALSE, trim_ws = TRUE)
colnames(pop)[1:2] <- c("sample", "population")

# === Merge data ===
merged <- left_join(evec, pop, by = "sample")

# === Highlighted samples ===
highlight_colors <- c(
  CX113F    = "#E41A1C",  # red
  Neige_2_3 = "#377EB8",  # blue
  Neige_3_2 = "#4DAF4A",  # green
  Neige_4   = "#984EA3"   # purple
)

highlight_labels <- c(
  CX113F    = "CX113F - Slovenia, 0.08x",
  Neige_2_3 = "Neige_2_3 - France, 27x",
  Neige_3_2 = "Neige_3_2 - France, 0.04x",
  Neige_4   = "Neige_4 - France, 0.01x"
)

highlight_samples <- names(highlight_colors)

merged <- merged %>%
  mutate(
    highlight = sample %in% highlight_samples,
    highlight_color = if_else(highlight, highlight_colors[sample], NA_character_)
  )

# === Define which populations to label ===
pop_to_label <- c("Italy", "Iberia", "Spain", "Croatia", "Scandinavia", "Bulgaria", "Slovakia", "Greece", "Central Asia", "Yakutia", "Azerbaijan")  # replace with your desired population names

# === Filter population centroids for labeling ===
pop_labels <- merged %>%
  group_by(population) %>%
  summarise(
    across(starts_with("PC"), ~mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(label = population) %>%
  filter(population %in% pop_to_label)

# === Common theme ===
base_theme <- theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# === Helper function for PCA plotting with labels ===
plot_pca <- function(merged, var_explained, pc_x, pc_y, file_name, highlight_colors, highlight_labels) {
  
  p <- ggplot() +
    # Regular samples
    geom_point(
      data = merged %>% filter(!highlight),
      aes(x = .data[[pc_x]], y = .data[[pc_y]], color = population),
      size = 3
    ) +
    # Highlighted samples
    geom_point(
      data = merged %>% filter(highlight),
      aes(x = .data[[pc_x]], y = .data[[pc_y]], fill = sample),
      shape = 24, color = "black", size = 4, stroke = 1.2
    ) +
    # Highlighted samples legend
    scale_color_discrete(name = "Population") +
    scale_fill_manual(
      name = "Highlighted samples",
      values = highlight_colors,
      labels = highlight_labels,
      guide = guide_legend(override.aes = list(shape = 24, size = 5, color = "black"))
    ) +
    # Population labels at centroids
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
    # Axis labels
    labs(
      x = paste0(pc_x, " (", round(var_explained[as.numeric(sub("PC", "", pc_x))], 1), "%)"),
      y = paste0(pc_y, " (", round(var_explained[as.numeric(sub("PC", "", pc_y))], 1), "%)")
    ) +
    base_theme
  
  show(p)
  ggsave(file_name, plot = p, width = 11, height = 5, dpi = 300)
}

# === List of PC pairs to plot ===
pc_pairs <- list(
  c("PC1", "PC2"),
  c("PC2", "PC3"),
  c("PC1", "PC3"),
  c("PC1", "PC4"),
  c("PC1", "PC5")
)

# === Loop through PC pairs and plot ===
for (pair in pc_pairs) {
  pc_x <- pair[1]
  pc_y <- pair[2]
  file_name <- paste0("wolf_pca_", pc_x, "_", pc_y, ".png")
  plot_pca(merged, var_explained, pc_x, pc_y, file_name, highlight_colors, highlight_labels)
}
