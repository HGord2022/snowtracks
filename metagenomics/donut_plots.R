# === Load libraries ===
library(tidyverse)
library(patchwork)
setwd("~/github_repos/snowtracks/metagenomics")

# === Read CSV ===
df <- read.csv("deep_stats.csv", stringsAsFactors = FALSE)

# Remove % signs, convert to numeric, NORMALISE to 100
df_long <- df %>%
  pivot_longer(-Name, names_to = "Category", values_to = "Percent") %>%
  mutate(Percent = as.numeric(sub("%","",Percent))) %>%
  group_by(Name) %>%
  mutate(Percent = Percent / sum(Percent) * 100) %>%
  ungroup()

# Set a consistent order of categories
df_long$Category <- factor(df_long$Category, levels = c(
  "Chordate_reads","Artificial_reads","Unclassified_reads","Microbial_reads",
  "Bacterial_reads","Viral_reads","Fungal_reads","Protozoan_reads"))

# Color palette
cols <- c(
  "Chordate_reads" = "#1b9e77",
  "Artificial_reads" = "#d95f02",
  "Unclassified_reads" = "#7570b3",
  "Microbial_reads" = "#e7298a",
  "Bacterial_reads" = "#66a61e",
  "Viral_reads" = "#e6ab02",
  "Fungal_reads" = "#a6761d",
  "Protozoan_reads" = "#666666"
)

# Helper to make a single donut chart
make_donut <- function(df_subset, show_legend = FALSE){
  ggplot(df_subset, aes(x = 2, y = Percent, fill = Category)) +
    geom_col(width = 1, color = NA) +     # solid fills, no borders
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    scale_y_continuous(expand = c(0, 0)) + # CRITICAL: removes gaps
    theme_void() +
    theme(legend.position = if(show_legend) "right" else "none") +
    scale_fill_manual(values = cols) +
    ggtitle(unique(df_subset$Name))
}

# === Split data for plotting ===
plots <- df_long %>%
  split(.$Name) %>%
  imap(~ make_donut(.x, show_legend = .y == "CX113F"))

# === Layout with extra vertical spacing ===
final_plot <-
  wrap_plots(plots[1:4], nrow = 1) /
  plot_spacer() /
  wrap_plots(plots[5:9], nrow = 1) +
  plot_layout(
    heights = c(1, 0.25, 1),   # controls gap size
    guides = "collect"
  )

# === Print figure ===
final_plot

ggsave(
  filename = "donut_plots.png",
  plot = final_plot,
  width = 14,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

