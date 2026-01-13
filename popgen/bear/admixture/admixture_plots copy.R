setwd("~/github_repos/snowtracks/popgen/bear/admixture")

library(tidyverse)
library(scales)

# === Load static data ===
fam  <- read.table("bears.fam", header = FALSE)
pops <- read.csv("populations.csv", header = FALSE)

# --- Extract sample names in ADMIXTURE order ---
samples <- fam$V2

# --- Build population lookup table ---
pop_df <- pops %>%
  transmute(
    sample = V1,
    population = V2
  )

# --- Population order for plotting ---
pop_order <- readLines("pop_order.txt")

# === Helper function for plotting ADMIXTURE results ===
plot_admixture <- function(K) {
  
  # --- Load Q file ---
  Q <- read.table(paste0("bears.", K, ".Q"), header = FALSE)
  
  n_clusters <- ncol(Q)
  
  # --- Manually specify default ggplot colours ---
  cluster_cols <- hue_pal()(n_clusters)
  names(cluster_cols) <- paste0("V", seq_len(n_clusters))
  
  cluster_cols["V10"] <- "#D73027"
  
  # --- Combine Q + metadata safely ---
  tbl <- Q %>%
    mutate(sample = samples) %>%
    left_join(pop_df, by = "sample")
  
  # --- Prepare data for plotting ---
  plot_data <- tbl %>%
    pivot_longer(
      cols = starts_with("V"),
      names_to = "cluster",
      values_to = "prob"
    ) %>%
    group_by(sample) %>%
    mutate(
      likely_assignment = cluster[which.max(prob)],
      assignment_prob   = max(prob)
    ) %>%
    ungroup() %>%
    arrange(
      factor(population, levels = pop_order),
      desc(assignment_prob)
    ) %>%
    mutate(
      sample          = factor(sample, levels = unique(sample)),
      cluster         = factor(cluster, levels = paste0("V", 1:n_clusters)),
      population_plot = factor(population, levels = pop_order)
    )
  
  # --- Plot ---
  p <- ggplot(plot_data, aes(x = sample, y = prob, fill = cluster)) +
    geom_col(width = 1) +
    facet_grid(~population_plot, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = cluster_cols) +
    theme_minimal() +
    theme(
      axis.text.x   = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      axis.ticks.x  = element_blank(),
      panel.spacing = unit(0.2, "lines"),
      strip.text    = element_text(size = 9, angle = 70),
      legend.position = "none",
      plot.margin   = margin(10, 10, 10, 10),
      strip.clip    = "off"
    ) +
    labs(
      x = "",
      y = "Ancestry proportion",
      title = paste("K =", K)
    )
  
  print(p)
  
  ggsave(
    filename = paste0("bear_plot_K", K, ".png"),
    plot = p,
    width = 12, height = 5, units = "in", dpi = 300
  )
}

# === Run ONLY for K = 10 ===
plot_admixture(10)
