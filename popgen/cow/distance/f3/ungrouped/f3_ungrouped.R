#!/usr/bin/env Rscript
setwd("~/github_repos/snowtracks/popgen/cow/distance/f3/ungrouped")
args <- commandArgs(trailingOnly = TRUE)
library(ggplot2)

# parse args
file_path <- "out_f3.txt"
outdir <- "plots"

if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# read data
f3_data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

samples <- unique(f3_data$pop2)

for(focal in samples) {
  
  # subset for this focal
  f3_values <- f3_data[f3_data$pop2 == focal & !is.na(f3_data$est), ]
  
  # order pop3 by decreasing est
  f3_values$pop3_label <- paste0(f3_values$pop3, " (", f3_values$pop3, ")")
  f3_values <- f3_values[order(f3_values$est, decreasing = TRUE), ]
  
  # take top 50
  top50 <- head(f3_values, 50)
  
  # make pop3 a factor with levels in order for plotting
  top50$pop3_label <- factor(top50$pop3_label, levels = rev(top50$pop3_label))
  
  # plot
  p <- ggplot(top50, aes(x = est, y = pop3_label)) +
    geom_errorbarh(aes(xmin = est - se, xmax = est + se), height = 0.2, color = "darkgray") +
    geom_point(size = 3, color = "steelblue") +
    labs(
      x = paste0("f3 (", focal, "; population)"),
      y = "Population"
    ) +
    theme_minimal() +
    theme(axis.text = element_text(size = 12))
  
  # save
  ggsave(
    filename = file.path(outdir, paste0("f3_", focal, ".png")),
    plot = p,
    width = 6,
    height = 9,
    dpi = 300
  )
}
