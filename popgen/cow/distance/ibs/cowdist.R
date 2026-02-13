# === Load libraries ===
library(tidyverse)
setwd("~/github_repos/snowtracks/popgen/cow/distance")

# === File paths ===
ibs_file <- "ibs_distances.csv"    # columns: sample, ibs_distance
pop_file <- "populations.csv"      # columns: sample, population

# === Read data ===
ibs <- read_csv(ibs_file, col_names = c("sample", "ibs_distance"))
pop <- read_csv(pop_file, col_names = c("sample", "population"))

# === Merge IBS distances with population info ===
merged <- left_join(ibs, pop, by = "sample")

# === Compute mean IBS distance per population ===
pop_avg <- merged %>%
  group_by(population) %>%
  summarise(mean_ibs = mean(ibs_distance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_ibs)) %>%
  slice_head(n = 10)

# === Reorder populations by mean_ibs for plotting ===
pop_avg <- pop_avg %>%
  mutate(population = fct_reorder(population, mean_ibs))

# === Create plot ===
p <- ggplot(pop_avg, aes(x = mean_ibs, y = population)) +
  geom_point(size = 5, color = "#00008B") +  # dark blue
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +  # 2 decimals
  labs(x = "Average IBS distance", y = "Population") +
  theme_classic(base_size = 14) +  # keeps axes lines
  theme(
    panel.grid.major = element_blank(),  # remove any major grid lines
    panel.grid.minor = element_blank(),  # remove any minor grid lines
    plot.title = element_blank()         # no title
  )

# === Show plot in window ===
print(p)

# === Save plot to file ===
ggsave("top10_populations_ibs.png", plot = p, width = 8, height = 6, dpi = 300)
