setwd("~/github_repos/snowtracks/popgen/wolf/new_data/distance/f3")

library(dplyr)
ind <- read.table("data/wolf_old.ind", stringsAsFactors = FALSE, col.names = c("ind", "sex", "pop"))
colnames(ind) <- c("ind", "sex", "pop")
# Strip whitespace first
ind$ind <- trimws(ind$ind)
# Now remove leading 0:
ind$ind <- sub("^0:", "", ind$ind)

map <- read.table("pop_map.txt", stringsAsFactors = FALSE, col.names = c("ind", "new_pop"))
ind <- ind %>% left_join(map, by = "ind") %>% mutate(pop = new_pop) %>% select(ind, sex, pop)
write.table(ind, "wolf_updated.ind", quote = FALSE, row.names = FALSE, col.names = FALSE)


library(admixtools)
library(tidyverse)
library(dplyr)


prefix <- "data/wolf"


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
# Example columns: V1 = sample, V3 = population
# Match samples in inds with populations and replace pop column
inds$pop <- pop_df$V3[match(inds$ind, pop_df$V1)]
# Sanity check
head(inds)
all(!is.na(inds$pop))  # Should be TRUE


outgroup <- "Wolf02"

# extract f2 blocks
f2_dir = 'f2data_n23/'
pops <- inds$pop[!inds$pop %in% c("CX113H", "CX113F", "Neige_3_2", "Neige_4", "Europe")]


extract_f2(prefix, f2_dir, 
           auto_only = FALSE, 
           overwrite = TRUE,
           pops=pops,
           maxmem = 20000,
           n_cores = 4)

f2_blocks = f2_from_precomp(f2_dir)


# f3 stats

outgroup <- "Chukotka"
pop1 = c('Neige_2_3')
pop2 <- inds$pop[!inds$pop %in% c("Neige_2_3", "CX113H", "CX113F", "Neige_3_2", "Neige_4", "Europe")]

f3_values <- f3(f2_blocks, outgroup, pop1, pop2)


# Take top 15 by est
top15 <- f3_values %>% 
  filter(!is.na(est)) %>% 
  arrange(desc(est)) %>% 
  slice(1:15)

# Horizontal point plot with error bars
ggplot(top15, aes(x = est, y = fct_reorder(pop3, est))) +
  geom_errorbarh(aes(xmin = est - se, xmax = est + se), height = 0.2, color = "darkgray") +
  geom_point(size = 3, color = "steelblue") +
  labs(x = "f3 (Chukotka, Neige_2_3; population)", y = "Population") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12))


# save
ggsave("f3_n23.png", width = 8, height = 6, dpi = 300)

