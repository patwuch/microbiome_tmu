library(phyloseq)
library(data.table)
library(biomformat)
library(tidytree)
library(treeio)
library(ALDEx2)
library(ape)

# --- Helper function for safe writing ---
write_safe <- function(obj, filename, write_func) {
  filepath <- file.path(outdir, filename)
  tryCatch({
    write_func(obj, filepath)
  }, error = function(e) {
    message("⚠️ Could not write file: ", filepath)
    message("Error: ", e$message)
  })
}

# --- Edit this if needed ---
base <- "/home/patwuch/projects/microbiome/experiments/yang"
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
# ============================================================
#        ALDEx2 Pairwise Analysis Within Each Metadata Level
# ============================================================

library(ALDEx2)
library(data.table)
library(ape)

# --- Extract data from phyloseq ---
counts <- as.data.frame(otu_table(ps))
meta <- as.data.frame(sample_data(ps))
tax <- as.data.frame(tax_table(ps))

# --- Filter rare taxa ---
min_prev <- 0.05
taxa_keep <- rowSums(counts > 0) / ncol(counts) >= min_prev
counts <- counts[taxa_keep, ]
tax <- tax[taxa_keep, ]

# --- Less strict significance threshold ---
p_threshold <- 0.10

# --- Metadata variables to analyze ---
meta_vars <- c("Group", "Disease", "Butyrate")

# --- Create output folder ---
outdir <- "~/microbiome_ALDEx2_results"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = TRUE)

if (file.access(outdir, 2) != 0) stop("Cannot write to folder ", outdir)

cat("All outputs will be saved in:", outdir, "\n")


# ============================================================
#              Loop over each metadata variable
# ============================================================

for (var in meta_vars) {
  cat("\n==== Running pairwise ALDEx2 for:", var, "====\n")
  
  conds <- as.factor(meta[[var]])
  levels_var <- levels(conds)
  
  if (length(levels_var) < 2) {
    cat("Skipping", var, "- only one level found.\n")
    next
  }
  
  # --- All pairwise combinations of levels ---
  pairwise_comps <- combn(levels_var, 2, simplify = FALSE)
  
  for (pair in pairwise_comps) {
    group1 <- pair[1]
    group2 <- pair[2]
    comp_name <- paste(var, group1, "vs", group2, sep = "_")
    cat("Processing:", comp_name, "\n")
    
    # Subset counts and conditions to just these two levels
    keep <- conds %in% c(group1, group2)
    counts_sub <- counts[, keep]
    conds_sub <- as.character(conds[keep])  # ALDEx2 requires character vector
    
    # --- ALDEx2 analysis ---
    x <- aldex.clr(counts_sub, conds_sub, mc.samples = 128)
    res <- aldex.ttest(x)
    eff <- aldex.effect(x)
    res <- data.frame(res, eff)
    
    # --- Filter significant taxa ---
    sig <- res[res$we.eBH < p_threshold, ]
    
    if (nrow(sig) > 0) {
      sig$Size <- abs(sig$effect)
      sig$Color <- ifelse(sig$effect > 0, "red", "blue")
      sig$taxon <- rownames(sig)
      
      # Merge with taxonomy
      sig <- merge(sig, tax, by.x = "taxon", by.y = "row.names", all.x = TRUE)
      
      # Ensure Rank1-Rank7 exist
      rank_cols <- paste0("Rank", 1:7)
      missing_cols <- setdiff(rank_cols, colnames(sig))
      if (length(missing_cols) > 0) sig[missing_cols] <- NA
      
      # --- Prepare phylo tree ---
      sig_tax <- sig[, rank_cols]
      rownames(sig_tax) <- sig$taxon
      sig_tax[, rank_cols] <- lapply(sig_tax[, rank_cols], function(x) {
        x <- as.character(x)
        x[is.na(x)] <- "Unclassified"
        factor(x)
      })
      # Ensure unique tip labels
      sig_tax$Rank7 <- make.unique(as.character(sig_tax$Rank7))

      # Convert missing ranks to Unclassified (you already do this)
      sig_tax[, rank_cols] <- lapply(sig_tax[, rank_cols], function(x) {
        x <- as.character(x)
        x[is.na(x)] <- "Unclassified"
        factor(x)
      })

      # --- Save GraPhlAn annotation ---
      sig$Node <- sig$taxon
      anno <- sig[, c("Node", "Color", "Size")]
      write_safe(anno, paste0("graphlan_annotation_", comp_name, ".txt"), function(x, f) {
        fwrite(x, f, sep = "\t", quote = FALSE, row.names = FALSE)
      })
      
      # --- Save full significant taxa table ---
      write_safe(sig, paste0("sig_taxa_", comp_name, ".csv"), function(x, f) {
        fwrite(x, f, row.names = FALSE)
      })
      
      cat("→ Significant taxa saved for", comp_name, "\n")
      
    } else {
      cat("No significant taxa for", comp_name, "\n")
    }
  }
}
