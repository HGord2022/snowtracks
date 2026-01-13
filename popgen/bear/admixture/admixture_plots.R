setwd("~/github_repos/snowtracks/popgen/bear/admixture")

library(tidyverse)

# === Load static data ===
fam  <- read.table("bears.fam", header = FALSE)
pops <- read.csv("populations.csv", header = FALSE)

# --- Extract sample names in ADMIXTURE order ---
samples <- fam$V2

# --- Build population lookup table ---
pop_df <- pops %>%
  transmute(
    sample = V1,      # sample name
    population = V2   # population label
  )

# --- Population order for plotting ---
# --- Population order for plotting (from file) ---
pop_order <- readLines("pop_order.txt")

# === Helper function for plotting ADMIXTURE results ===
plot_admixture <- function(K) {
  
  # --- Load Q file ---
  Q <- read.table(paste0("bears.", K, ".Q"), header = FALSE)
  
  # --- Combine Q + metadata safely ---
  tbl <- Q %>%
    mutate(sample = samples) %>%
    left_join(pop_df, by = "sample")
  
  n_clusters <- ncol(Q)
  
  # --- Prepare data for plotting ---
  plot_data <- tbl %>%
    tidyr::pivot_longer(
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
  
  # --- Save plot ---
  ggsave(
    filename = paste0("bear_plot_K", K, ".png"),
    plot = p,
    width = 12, height = 5, units = "in", dpi = 300
  )
}

# === Run for selected K values ===
Ks <- c(2, 3, 4, 5, 10)
purrr::walk(Ks, plot_admixture)
