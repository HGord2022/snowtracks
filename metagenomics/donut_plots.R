# === Load libraries ===
library(tidyverse)
setwd("~/github_repos/snowtracks/metagenomics")

# === Read CSV ===
df <- read.csv("deep_stats.csv", stringsAsFactors = FALSE)

# Remove % signs, convert to numeric, filter categories, NORMALISE to 100
df_long <- df %>%
  pivot_longer(-Name, names_to = "Category", values_to = "Percent") %>%
  mutate(Percent = as.numeric(sub("%", "", Percent))) %>%
  filter(Category %in% c("Chordate_reads","Unclassified_reads","Microbial_reads")) %>%
  group_by(Name) %>%
  mutate(Percent = Percent / sum(Percent) * 100) %>%
  ungroup()

# Set category order (BOTTOM → TOP of bars)
df_long$Category <- factor(
  df_long$Category,
  levels = c("Microbial_reads","Chordate_reads","Unclassified_reads")
)

# Color palette (subset)
cols <- c(
  "Chordate_reads" = "#1b9e77",
  "Unclassified_reads" = "#7570b3",
  "Microbial_reads" = "#e7298a"
)

# === Stacked vertical barplot ===
final_plot <- ggplot(df_long, aes(x = Name, y = Percent, fill = Category)) +
  geom_col(width = 0.6, color = NA) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Percent of reads", x = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

# === Print figure ===
final_plot

# === Save figure ===
ggsave(
  filename = "stacked_barplots.png",
  plot = final_plot,
  width = 12,
  height = 5,
  units = "in",
  dpi = 300,
  bg = "white"
)
