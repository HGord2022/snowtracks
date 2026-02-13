setwd("~/github_repos/snowtracks/popgen/cow/distance/f3/ungrouped")

library(admixtools)
library(tidyverse)
library(dplyr)

prefix = 'data/cows'

inds <- read.table(paste0(prefix, ".ind"), stringsAsFactors = FALSE, header = FALSE)
colnames(inds) <- c("ind", "sex", "pop")
# Strip whitespace first
inds$ind <- trimws(inds$ind)
# Now remove leading 0:
inds$ind <- sub("^0:", "", inds$ind)
# Optional: make pop match ind
inds$pop <- inds$ind
all_samples <- inds$ind

# Read populations file
pop_df <- read.delim("ena_metadata.txt", header = FALSE, stringsAsFactors = FALSE)
# Check the structure
head(pop_df)
# Example columns: V1 = sample, V2 = population
# Match samples in inds with populations and replace pop column
inds$pop <- pop_df$V2[match(inds$ind, pop_df$V1)]



# extract f2 blocks
f2_dir = 'f2data/'
pops <- inds$ind


extract_f2(prefix, f2_dir, 
           auto_only = FALSE,
           maxmiss = 1,
           pops=pops,
           maxmem = 20000,
           n_cores = 4,
           adjust_pseudohaploid = FALSE)

f2_blocks = f2_from_precomp(f2_dir, remove_na = FALSE)


# f3 stats

samples <- c("CX1138.fixed.fa.gz")

for (focal in samples) {
  
  outgroup <- "SAMN08242229"
  pop1 <- c(focal)
  pop2 <- inds$ind[!inds$ind %in% focal]
  
  f3_values <- f3(f2_blocks, outgroup, pop1, pop2) %>%
    left_join(
      inds %>% select(ind, pop),
      by = c("pop3" = "ind")
    ) %>%
    mutate(pop3_label = paste0(pop3, " (", pop, ")"))
  
  
  # Take top 50 by est
  top50 <- f3_values %>% 
    filter(!is.na(est)) %>% 
    arrange(desc(est)) %>% 
    slice(1:50)
  
  # Horizontal point plot with error bars
  p <- ggplot(top50, aes(x = est, y = fct_reorder(pop3_label, est))) +
    geom_errorbarh(
      aes(xmin = est - se, xmax = est + se),
      height = 0.2,
      color = "darkgray"
    ) +
    geom_point(size = 3, color = "steelblue") +
    labs(
      x = paste0("f3 (Zebu, CX1138; population)"),
      y = "Population"
    ) +
    theme_minimal() +
    theme(axis.text = element_text(size = 12))
  show(p)
  # save
  ggsave(
    filename = paste0("f3_", focal, "_ind.png"),
    plot = p,
    width = 6,
    height = 9,
    dpi = 300
  )
}
