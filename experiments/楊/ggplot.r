if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ggtree")

library(ggtree)
library(phyloseq)


aldex_res <- read.table("/home/patwuch/projects/microbiome/experiments/楊/ALDEx2_results/sig_taxa_comp4.csv", header=TRUE, sep=",")
taxonomy <- read.table("experiments/楊/exported/exported-taxonomy/taxonomy.tsv", header=TRUE, sep="\t")
merged <- merge(aldex_res, taxonomy, by="Feature")

# Load tree and attach taxonomic data
ps <- phyloseq(otu_table(otu_table), tax_table(tax_table), phy_tree(tree))

# Extract tree
tree <- phy_tree(ps)

# Merge ALDEx2 results
tree_data <- data.frame(tax_table(ps))
tree_data$Feature <- rownames(tree_data)
tree_data <- merge(tree_data, aldex_res, by="Feature", all.x=TRUE)

# Plot
ggtree(tree) %<+% tree_data +
  geom_tippoint(aes(color = -log10(we.ep), size = abs(diff.btw))) +
  scale_color_viridis_c() +
  theme_tree2()
