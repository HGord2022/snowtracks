setwd("~/github_repos/snowtracks/popgen/bear/distance/f3/grouped")

library(admixtools)
library(tidyverse)
library(dplyr)

## create grouped populations in the ind file
ind <- read.table("data/bears_old.ind", stringsAsFactors = FALSE, col.names = c("ind", "sex", "pop"))
colnames(ind) <- c("ind", "sex", "pop")
# Strip whitespace first
ind$ind <- trimws(ind$ind)
# Now remove leading 0:
ind$ind <- sub("^0:", "", ind$ind)

map <- read.table("pop_map.txt", stringsAsFactors = FALSE, col.names = c("ind", "new_pop"))
ind <- ind %>% left_join(map, by = "ind") %>% mutate(pop = new_pop) %>% select(ind, sex, pop)
write.table(ind, "data/bears.ind", quote = FALSE, row.names = FALSE, col.names = FALSE)

prefix <- "data/bears"



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
pop_df <- read.delim("populations.txt", header = FALSE, stringsAsFactors = FALSE)
# Check the structure
head(pop_df)
# Example columns: V2 = sample, V4 = population
# Match samples in inds with populations and replace pop column
inds$pop <- pop_df$V4[match(inds$ind, pop_df$V2)]



# extract f2 blocks
f2_dir = 'f2data/'
pops <- inds$pop


extract_f2(prefix, f2_dir, 
           auto_only = FALSE,
           maxmiss = 1,
           pops=pops,
           maxmem = 20000,
           n_cores = 4,
           adjust_pseudohaploid = FALSE)

f2_blocks = f2_from_precomp(f2_dir, remove_na = FALSE)


# f3 stats

samples <- c("CX115H", "CX1138", "CX113C", "CX113E")

for (focal in samples) {
  
  outgroup <- "Polar_Bear"
  pop1 <- c(focal)
  pop2 <- inds$pop[!inds$pop %in% focal]
  
  f3_values <- f3(f2_blocks, outgroup, pop1, pop2)
  
  
  # Take top 15 by est
  top15 <- f3_values %>% 
    filter(!is.na(est)) %>% 
    arrange(desc(est)) %>% 
    slice(1:15)
  
  # Horizontal point plot with error bars
  p <- ggplot(top15, aes(x = est, y = fct_reorder(pop3, est))) +
    geom_errorbarh(
      aes(xmin = est - se, xmax = est + se),
      height = 0.2,
      color = "darkgray"
    ) +
    geom_point(size = 3, color = "steelblue") +
    labs(
      x = paste0("f3 (Polar Bear, ", focal, "; population)"),
      y = "Population"
    ) +
    theme_minimal() +
    theme(axis.text = element_text(size = 12))
  
  show(p)
  
  # save
  ggsave(
    filename = paste0("f3_", focal, ".png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}
