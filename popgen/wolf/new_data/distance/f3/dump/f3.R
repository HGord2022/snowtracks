setwd("~/github_repos/snowtracks/popgen/wolf/new_data/distance/f3")
# ============================================================
# Outgroup f3 analysis with admixtools
# Produces a dataframe:
#   rows   = sample B (all samples)
#   cols   = test samples A
#   values = f3(outgroup; A, B)
# ============================================================

library(admixtools)
library(tidyverse)
library(dplyr)

# ------------------------------------------------------------
# 1) DATA PREFIX (EIGENSTRAT format)
# ------------------------------------------------------------
# expects:
#   prefix.geno
#   prefix.snp
#   prefix.ind
prefix <- "data/wolf"

# ------------------------------------------------------------
# 2) READ SAMPLE LIST
# ------------------------------------------------------------
inds <- read.table(paste0(prefix, ".ind"), stringsAsFactors = FALSE, header = FALSE)
colnames(inds) <- c("ind", "sex", "pop")
# Strip whitespace first
inds$ind <- trimws(inds$ind)
# Now remove leading 0:
inds$ind <- sub("^0:", "", inds$ind)
# Optional: make pop match ind
inds$pop <- inds$ind
all_samples <- inds$ind


# Read your populations file
pop_df <- read.delim("populations.txt", header = FALSE, stringsAsFactors = FALSE)
# Check the structure
head(pop_df)
# Example columns: V1 = sample, V2 = population
# Match samples in inds with populations and replace pop column
inds$pop <- pop_df$V2[match(inds$ind, pop_df$V1)]
# Sanity check
head(inds)
all(!is.na(inds$pop))  # Should be TRUE

# ------------------------------------------------------------
# 3) DEFINE TEST SAMPLES (A) AND OUTGROUP (O)
# ------------------------------------------------------------
test_samples <- c(
  "Neige_2_3",
  "Neige_3_2",
  "Neige_4",
  "CX113H",
  "CX113F"
)

outgroup <- "Wolf02"

# ------------------------------------------------------------
# 4) DEFINE SAMPLE B POOL
# ------------------------------------------------------------
# Includes *all* samples, including test samples
B_samples <- all_samples

# extract f2 blocks
f2_dir = 'f2data/'
pops <- inds$ind[!inds$ind %in% c("CX113H", "CX113F", "Neige_3_2", "Neige_4", "Europe")]


extract_f2(prefix, f2_dir, 
           auto_only = FALSE, 
           overwrite = TRUE,
           pops=pops,
           maxmem = 20000,
           n_cores = 4)

f2_blocks = f2_from_precomp(f2_dir)


# f3 stats

outgroup <- "Wolf02"
pop1 = c('Neige_2_3')
pop2 <- inds$ind[!inds$ind %in% c("Neige_2_3", "CX113H", "CX113F", "Neige_3_2", "Neige_4", "Europe")]

f3_values <- f3(f2_blocks, outgroup, pop1, pop2)

f3_values$pop3 <- pop_df$V2[match(inds$ind, pop_df$V1)]

# keep only sample ID + population (2nd column)
pops <- pop_df[, c(1, 2)]
colnames(pops) <- c("pop3", "population")

# join + select
f3_pop <- f3_values %>%
  left_join(pops, by = "pop3") %>%
  select(
    population,
    est,
    se
  )

f3_pop %>%
  group_by(population) %>%                             # average across population
  summarise(est_mean = mean(est, na.rm = TRUE),
            se_mean = mean(se, na.rm = TRUE)) %>%
  arrange(desc(est_mean)) %>%                          # sort by mean
  slice_head(n = 15) %>%                              # take top 15
  ggplot(aes(x = est_mean, y = fct_reorder(population, est_mean))) +
  geom_point(size = 3, color = "steelblue") +         # points instead of bars
  geom_errorbarh(aes(xmin = est_mean - se_mean, xmax = est_mean + se_mean),
                 height = 0.2, color = "gray50") +    # optional: error bars
  labs(x = "Mean f3 (East Russia, Neige_2_3; population)", 
       y = "") +
  theme_minimal(base_size = 14)

ggsave("f3_n23.png", width = 8, height = 6, dpi = 300)


# ------------------------------------------------------------
# 5) HELPER FUNCTION: outgroup f3
# ------------------------------------------------------------
# Computes f3(O; A, B)
run_f3 <- function(A, B, O, prefix) {
  f3(
    A = A,
    B = B,
    C = O,
    data = prefix
  )$f3
}

# ------------------------------------------------------------
# 6) RUN ALL A × B COMBINATIONS
# ------------------------------------------------------------
f3_df <- map_dfc(
  test_samples,
  function(A) {
    tibble(
      sample_B = B_samples,
      !!A := map_dbl(
        B_samples,
        ~ run_f3(A, .x, outgroup, prefix)
      )
    ) %>% select(-sample_B)
  }
)

f3_df <- bind_cols(
  tibble(sample_B = B_samples),
  f3_df
)

# ------------------------------------------------------------
# 7) VIEW RESULT
# ------------------------------------------------------------
print(f3_df)

# ------------------------------------------------------------
# 8) SAVE TO FILE (OPTIONAL)
# ------------------------------------------------------------
write_tsv(f3_df, "outgroup_f3_results.tsv")
