library(ape)
library(ggtree)
setwd("~/github_repos/snowtracks/plotting/tree_plots")

tree <- read.tree("core_tree.treefile")

ggtree(tree) +
  geom_tiplab(size = 3)

library(ape)
library(ggtree)

tree <- read.tree("mytree.treefile")

p <- ggtree(tree) +
  geom_tiplab(size = 2)

ggsave("tree.jpeg", p, width = 10, height = 15)
