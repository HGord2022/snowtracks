setwd("~/github_repos/snowtracks/popgen/wolf/admixture")

library(tidyverse)

# --- Load data ---
Q <- read.table("wolf_final.10.Q", header = FALSE)
fam <- read.table("wolf_final.fam", header = FALSE)
pops <- read.table("populations.txt", header = FALSE, sep = "\t")

samples <- fam$V2
pop_labels <- pops$V3  # population names

tbl <- Q %>%
  mutate(sample = samples,
         population = pop_labels)

# --- Population order for plotting ---
pop_order <- c("Yakutia","Mongolia","Altai","Azerbaijan","Georgia","Bulgaria",
               "Greece","Croatia","Slovakia","Ukraine","Poland","Russia","Latvia", "Belarus",
               "Finland","Europe","France","Slovenia")

# --- Prepare data ---
n_clusters <- ncol(Q)

plot_data <- tbl %>%
  pivot_longer(cols = starts_with("V"), names_to = "cluster", values_to = "prob") %>%
  group_by(sample) %>%
  mutate(
    likely_assignment = cluster[which.max(prob)],
    assignment_prob = max(prob)
  ) %>%
  ungroup() %>%
  # Arrange by population order only
  arrange(factor(population, levels = pop_order), desc(assignment_prob)) %>%
  mutate(
    sample = factor(sample, levels = unique(sample)),
    cluster = factor(cluster, levels = paste0("V", 1:n_clusters)),
    population_plot = factor(population, levels = pop_order)
  )

# --- Final plot ---
p_final <- ggplot(plot_data, aes(x = sample, y = prob, fill = cluster)) +
  geom_col(width = 1, position = "stack") +
  facet_grid(~population_plot, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.2, "lines"),
    strip.text = element_text(size = 9, angle = 70),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10),
    strip.clip = "off"
  ) +
  labs(
    x = "",
    y = "Ancestry proportion",
    title = paste("K =", n_clusters)
  )

# Show plot
p_final

# --- Dynamic file name and save ---
file_name <- paste0("wolf_plot_K", n_clusters, ".png")
ggsave(filename = file_name,
       plot = p_final,
       width = 12, height = 5, units = "in", dpi = 300)

