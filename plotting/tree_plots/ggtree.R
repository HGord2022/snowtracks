library(ape)
library(ggtree)
setwd("~/github_repos/snowtracks/plotting/tree_plots")


tree <- read.tree("core_tree.treefile")

# Root the tree on HPA21
tree <- root(tree, outgroup = "HPA21_ref", resolve.root = TRUE)

p <- ggtree(tree) +
  geom_tiplab(size = 2)

ggsave("tree.jpeg", p, width = 10, height = 15)
