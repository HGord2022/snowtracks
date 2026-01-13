

# ----------------------------------------
# Script: plot_pca.R
# Purpose: Plot PLINK PCA and highlight a sample
# Usage: Rscript plot_pca.R
# ----------------------------------------
setwd("~/github_repos/snowtracks/popgen")

# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)

# ----- Step 1: Read PLINK PCA files -----
pca_file <- "merged_plink.europe.pca.eigenvec"       # PLINK eigenvectors
eigval_file <- "merged_plink.europe.pca.eigenval"    # PLINK eigenvalues

pca <- read.table(pca_file, header = FALSE, stringsAsFactors = FALSE)
eigvals <- scan(eigval_file)

# PLINK .eigenvec format: first column is FID, second is IID, remaining are PCs
colnames(pca) <- c("FID", "Sample", paste0("PC", 1:(ncol(pca)-2)))

# Compute % variance explained
var_explained <- eigvals / sum(eigvals) * 100

# ----- Step 2: Read population dictionary (only first 2 columns) -----
pop_file <- "eu_pops.txt"   # Replace with your file path
pop_dict <- read_tsv(pop_file, col_names = FALSE) %>%
  select(1,2) %>%
  rename(Sample = X1, Population = X2)

# Merge PCA with population info
pca_pop <- left_join(pca, pop_dict, by = "Sample")

# ----- Step 3: Assign colors to populations dynamically -----
n_pops <- length(unique(pca_pop$Population))
pop_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_pops),
  unique(pca_pop$Population)
)

# ----- Step 4: Plot PCA -----
p <- ggplot(pca_pop, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Population), size = 2, alpha = 0.7) +                   # all points
  geom_point(data = subset(pca_pop, Sample == "Neige_2_3"),                      # highlight Neige_2_3
             aes(color = Population), size = 4, shape = 21, fill = "red") +
  geom_text(data = subset(pca_pop, Sample == "Wlf_22799"), aes(label = Sample), vjust = 1.5, size = 3) +  scale_color_manual(values = pop_colors) +
  theme_minimal(base_size = 14) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 2), "%)")
  )

# ----- Step 5: Display plot in RStudio -----
print(p)

# ----- Step 6: Save plot -----
ggsave("PCA_plot_europe.pdf", p, width = 8, height = 6)
