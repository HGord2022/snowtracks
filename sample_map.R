setwd("~/github_repos/snowtracks")

# Load libraries
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(patchwork)

# Your data
df <- data.frame(
  lon = c(
    rep(14.446978287833897),   # Bear x2
    14.46256570314588,            # Bear
    14.471102965596586,           # Bear
    rep(6.107745516021383, 4)     # Wolf x4
  ),
  lat = c(
    rep(45.58896851320194),
    45.62552097267587,
    45.68457290520403,
    rep(45.303352202944104, 4)
  ),
  species = c(
    rep("Bear", 3),
    rep("Wolf", 4)
  )
)

# Function to spread duplicates radially (both species)
spread_duplicates <- function(df, radius = 0.05) {
  df %>%
    group_by(lon, lat, species) %>%
    mutate(
      n_dup = n(),
      angle = ifelse(n_dup > 1,
                     seq(0, 2*pi, length.out = n_dup + 1)[1:n_dup],
                     0),
      adj_lon = lon + ifelse(n_dup > 1, radius * cos(angle), 0),
      adj_lat = lat + ifelse(n_dup > 1, radius * sin(angle), 0)
    ) %>%
    ungroup()
}

df_adj <- spread_duplicates(df)

# Map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Manual country labels
wolf_labels <- data.frame(
  lon = c(4.606217081716726, 7.20962407728195, 7.473695653694315),
  lat = c(45.74666562072359, 46.677299493195406, 45.25969436424621),
  name = c("France", "Switzerland", "Italy")
)

bear_labels <- data.frame(
  lon = c(14.467039197910825, 13.00098665300081, 14.157438039358402, 14.715808487720494),
  lat = c(46.11839859549004, 46.15625728965163, 46.67729948651543, 45.35878285476793),
  name = c("Slovenia", "Italy", "Austria", "Croatia")
)

# Function to plot species map
plot_species_map <- function(species_name, df, world, labels, zoom_km = 400) {
  subdf <- df %>% filter(species == species_name)
  
  clat <- mean(subdf$lat)
  clon <- mean(subdf$lon)
  
  deg_lat <- (zoom_km/111)/2
  deg_lon <- (zoom_km/85)/2
  
  lon_range <- c(clon - deg_lon, clon + deg_lon)
  lat_range <- c(clat - deg_lat, clat + deg_lat)
  
  ggplot(data = world) +
    geom_sf(fill = "antiquewhite", color = "grey40", size = 0.3) +
    geom_point(
      data = subdf,
      aes(x = adj_lon, y = adj_lat, color = species),
      size = 2, alpha = 0.9
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
      legend.position = "none"
    ) +
    ggtitle(species_name)
}

# Wolf map: regional zoom (400 km), Bear map: slightly more zoomed in (200 km)
map_wolf <- plot_species_map("Wolf", df_adj, world, wolf_labels, zoom_km = 400)
map_bear <- plot_species_map("Bear", df_adj, world, bear_labels, zoom_km = 200)

# Combine side by side
p <- map_wolf + map_bear

# Save
ggsave("sample_map.png", plot = p, width = 12, height = 6, dpi = 300)

p
