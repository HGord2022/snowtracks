setwd("~/github_repos/snowtracks")

# Clean environment
rm(list = ls())

# Load libraries
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# Your data
df <- data.frame(
  lon = c(
    rep(14.446978287833897, 2),   # Bear x2
    14.46256570314588,            # Bear
    14.471102965596586,           # Bear
    rep(6.107745516021383, 4)     # Wolf x4
  ),
  lat = c(
    rep(45.58896851320194, 2),
    45.62552097267587,
    45.68457290520403,
    rep(45.303352202944104, 4)
  ),
  species = c(
    rep("Bear", 4),
    rep("Wolf", 4)
  )
)

# Keep only one point per species
df_adj <- df %>%
  group_by(species) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(adj_lon = lon, adj_lat = lat)

# Map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Country labels
labels <- data.frame(
  lon = c(4.606217081716726, 14.467039197910825),
  lat = c(45.74666562072359, 46.12839859549004),
  name = c("France", "Slovenia")
)

# Compute map extent based on both points
clon <- mean(df_adj$lon)
clat <- mean(df_adj$lat)
zoom_km <- 1200   # adjust zoom
deg_lat <- (zoom_km/111)/2
deg_lon <- (zoom_km/85)/2
lon_range <- c(clon - deg_lon, clon + deg_lon)
lat_range <- c(clat - deg_lat, clat + deg_lat)

# Plot single map with both species
p <- ggplot(data = world) +
  geom_sf(fill = "antiquewhite", color = "grey40", size = 0.3) +
  geom_point(
    data = df_adj,
    aes(x = adj_lon, y = adj_lat, color = species),
    size = 3, alpha = 0.9
  ) +
  geom_text(
    data = labels,
    aes(x = lon, y = lat, label = name),
    size = 3.5, family = "Arial", color = "black"
  ) +
  scale_color_manual(values = c("Bear" = "red", "Wolf" = "blue")) +
  coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, hjust = 0.5, family = "Arial"),
    legend.position = "right"
  )

# Save
ggsave("sample_map.png", plot = p, width = 8, height = 6, dpi = 300)

# Show plot
p
