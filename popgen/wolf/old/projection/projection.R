
setwd("~/github_repos/snowtracks/popgen/projection")
# ---- R script: project_samples_to_refPCA.R ----
# Run in RStudio: plots will appear in the plot window

library(data.table)
library(dplyr)
library(ggplot2)

# ---- parameters ----
ref_raw <- "ref_pruned.raw"                # reference panel
query_raws <- c("query_pruned.raw")       # your query samples
pops_file <- "populations.txt"            # Sample Population table (first two columns)
out_plot <- "PCA_projected.png"
nPCs <- 10

# ---- read reference ----
ref <- fread(ref_raw)
geno_cols <- colnames(ref)[-(1:6)]
ref_ids <- ref[, .(FID = get(colnames(ref)[1]), IID = get(colnames(ref)[2]))]
ref_geno <- as.matrix(ref[, ..geno_cols])
mode(ref_geno) <- "double"

# ---- read populations ----
pops <- fread(pops_file, header = FALSE)[, .(Sample = as.character(V1), Population = as.character(V2))]

# ---- impute missing in reference ----
ref_means <- colMeans(ref_geno, na.rm = TRUE)
for(j in seq_len(ncol(ref_geno))) {
  nas <- is.na(ref_geno[, j])
  if(any(nas)) ref_geno[nas, j] <- ref_means[j]
}

# ---- center reference ----
ref_centered <- sweep(ref_geno, 2, ref_means, "-")

# ---- PCA on reference ----
pca_ref <- prcomp(ref_centered, center = FALSE, scale. = FALSE, rank. = nPCs)
var_explained <- (pca_ref$sdev^2) / sum(pca_ref$sdev^2) * 100
ref_scores <- as.data.frame(pca_ref$x)
ref_scores$Sample <- ref_ids$IID
ref_scores <- left_join(ref_scores, pops, by = c("Sample"))

# ---- function to project a query raw ----
project_sample <- function(rawfile) {
  s <- fread(rawfile)
  s_ids <- s[, .(FID = get(colnames(s)[1]), IID = get(colnames(s)[2]))]
  s_geno <- as.matrix(s[, -(1:6), with = FALSE])
  mode(s_geno) <- "double"
  
  # --- align SNPs with reference ---
  common_snps <- intersect(colnames(ref_geno), colnames(s_geno))
  cat("Query", basename(rawfile), "has", length(common_snps), "SNPs in common with reference\n")
  if(length(common_snps) == 0) {
    warning("No overlapping SNPs for query ", basename(rawfile))
    return(NULL)
  }
  
  s_geno_sub <- s_geno[, common_snps, drop = FALSE]
  ref_means_sub <- ref_means[common_snps]
  
  # impute missing in query with reference means
  nas <- is.na(s_geno_sub)
  s_geno_sub[nas] <- ref_means_sub[nas]
  
  # center
  s_centered <- sweep(s_geno_sub, 2, ref_means_sub, "-")
  
  # project
  proj_scores <- s_centered %*% pca_ref$rotation[common_snps, 1:nPCs, drop = FALSE]
  out <- as.data.frame(proj_scores)
  colnames(out) <- paste0("PC", 1:nPCs)
  out$Sample <- s_ids$IID
  return(out)
}

# ---- project all queries ----
proj_list <- lapply(query_raws, project_sample)
proj_df <- do.call(rbind, proj_list)
proj_df <- left_join(proj_df, pops, by = c("Sample"))

# ---- combine for plotting ----
plot_df <- bind_rows(
  ref_scores %>% mutate(Type = "Reference") %>% select(Sample, PC1, PC2, Population, Type),
  proj_df %>% mutate(Type = "Projected") %>% select(Sample, PC1, PC2, Population, Type)
)

# ---- plot PCA ----
plot_df_complete <- plot_df %>% filter(!is.na(PC1) & !is.na(PC2))

p <- ggplot(plot_df_complete, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(data = subset(plot_df_complete, Type == "Reference"), size = 2, alpha = 0.7) +
  geom_point(data = subset(plot_df_complete, Type == "Projected"), size = 4) +
  geom_text(data = subset(plot_df_complete, Type == "Projected"), aes(label = Sample), vjust = -1) +
  xlab(paste0("PC1 (", round(var_explained[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(var_explained[2], 2), "%)")) +
  theme_minimal()


print(p)   # shows in RStudio
ggsave(out_plot, p, width = 8, height = 6, dpi = 300)

# ---- save coordinates ----
write.table(
  plot_df,
  file = "pca_combined_coordinates.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
