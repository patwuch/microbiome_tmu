library(phyloseq)
library(data.table)
library(biomformat)
library(tidytree)
library(treeio)
library(ALDEx2)
library(ape)

# --- Edit this if needed ---
base <- "/home/patwuch/projects/microbiome/experiments/楊"
otu_fp <- file.path(base, "exported/exported-feature-table/table.from_biom.tsv")
tax_fp <- file.path(base, "exported/exported-taxonomy/taxonomy.tsv")
meta_fp <- file.path(base, "metadata.tsv")
# ----------------------------

# Read OTU count table
otu_counts <- read.delim(otu_fp,
                         header = TRUE,       # first non-skipped line is header
                         row.names = 1,
                         check.names = FALSE,
                         skip = 1)  # skip the first line starting with '#'

# Read taxonomy table
tax_mat = read.csv(tax_fp, sep="\t", row.names=1)
tax_mat = as.matrix(tax_mat)  # convert to matrix

# Read metadata
meta = read.csv(meta_fp, sep="\t", row.names=1)

# Make phyloseq components
OTU = otu_table(otu_counts, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
SAMP = sample_data(meta)

setdiff(colnames(OTU), rownames(meta))  # in OTU but not metadata
setdiff(rownames(meta), colnames(OTU))  # in metadata but not OTU

# Combine into phyloseq object
ps = phyloseq(OTU, TAX, SAMP)

head(otu_table(ps))
head(tax_table(ps))
head(sample_data(ps))

# --- ALDEx2 Analysis ---

# --- Extract data from phyloseq ---
counts <- as.data.frame(otu_table(ps))
conds <- sample_data(ps)$Group
tax <- as.data.frame(tax_table(ps))

# --- Filter rare taxa ---
# MODIFICATION 1: Decreased the minimum prevalence from 0.2 to 0.05 (5%)
# This keeps more rare taxa in the analysis.
min_prev <- 0.05
taxa_keep <- rowSums(counts > 0) / ncol(counts) >= min_prev
counts <- counts[taxa_keep, ]
tax <- tax[taxa_keep, ]

# --- Define comparisons ---
comparisons <- list(
  # comp1 = c("E","EH","EL","C","CH","CL"),      # multi-group
  # comp2 = c("E","EH","EL"),                   # multi-group
  # comp3 = c("C","CH","CL"),                   # multi-group
  comp4 = list(group1 = c("C","CH","CL"), group2 = c("E","EH","EL"))  # two-group
)

# MODIFICATION 2: Set a less strict p-value threshold (alpha)
# Increased from 0.05 to 0.10 for the corrected p-value (eBH)
p_threshold <- 0.10

for (i in names(comparisons)) {

  cat("Processing", i, "...\n")
  comp <- comparisons[[i]]

  if (i == "comp4") {
    # Two-group comparison
    conds_sub <- ifelse(conds %in% comp$group1, "Cgroup",
                        ifelse(conds %in% comp$group2, "Egroup", NA))
    keep <- !is.na(conds_sub)
    counts_sub <- counts[, keep]
    conds_sub <- conds_sub[keep]

    x <- aldex.clr(counts_sub, conds_sub, mc.samples = 128)
    res <- aldex.ttest(x)
    eff <- aldex.effect(x)
    res <- data.frame(res, eff)

    # MODIFICATION 2a: Use the less strict threshold for two-group t-test (we.eBH)
    sig <- res[res$we.eBH < p_threshold, ]

    if (nrow(sig) > 0) {
      sig$Size <- abs(sig$effect)
      sig$Color <- ifelse(sig$effect > 0, "red", "blue")
    }

  } else {
    # Multi-group comparison
    keep <- conds %in% comp
    counts_sub <- counts[, keep]
    conds_sub <- conds[keep]

    x <- aldex.clr(counts_sub, conds_sub, mc.samples = 128)
    res <- aldex.kw(x)
    # MODIFICATION 2b: Use the less strict threshold for multi-group Kruskal-Wallis (kw.eBH)
    sig <- res[res$kw.eBH < p_threshold, ]

    if (nrow(sig) > 0) {
      sig$Size <- 1        # placeholder for GraPhlAn
      sig$Color <- "grey"
    }
  }
  print(head(sig))
  print("Finished ALDEx2 analysis.")

    if (nrow(sig) > 0) {
    # Merge with taxonomy
    sig$taxon <- rownames(sig)
    sig <- merge(sig, tax, by.x = "taxon", by.y = "row.names", all.x = TRUE)
    
    # Ensure Rank columns exist
    rank_cols <- paste0("Rank", 1:7)
    missing_cols <- setdiff(rank_cols, colnames(sig))
    if(length(missing_cols) > 0){
      sig[missing_cols] <- NA
    }

    # Prepare sig_tax for tree
    sig_tax <- sig[, rank_cols]
    rownames(sig_tax) <- sig$taxon

    # Replace NAs with 'Unclassified' and convert to factors
    sig_tax[, rank_cols] <- lapply(sig_tax[, rank_cols], function(x) {
      x <- as.character(x)
      x[is.na(x)] <- "Unclassified"
      factor(x)
    })

    # Build phylo tree
    tree <- as.phylo(~ Rank1/Rank2/Rank3/Rank4/Rank5/Rank6/Rank7, data = sig_tax)
    write.tree(tree, paste0("graphlan_tree_", i, ".txt"))

    # Create Node column for GraPhlAn
    sig$Node <- sig$taxon  # Node = taxon name
    # Save GraPhlAn annotation
    anno <- sig[, c("Node","Color","Size")]
    fwrite(anno, paste0("graphlan_annotation_", i, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

    # Save full significant taxa table
    fwrite(sig, paste0("sig_taxa_", i, ".csv"), row.names = FALSE)
    
  } else {
    cat("No significant taxa for", i, "\n")
}

}

