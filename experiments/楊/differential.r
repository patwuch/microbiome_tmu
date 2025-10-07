library(phyloseq)
library(data.table)
library(biomformat)

# --- Edit this if needed ---
base <- "/home/patwuch/projects/microbiome/experiments/楊"
biom_fp <- file.path(base, "exported/exported-feature-table/feature-table.biom")
tax_fp  <- file.path(base, "exported/exported-taxonomy/taxonomy.tsv")
meta_fp <- file.path(base, "metadata.tsv")
# ----------------------------

clean_id <- function(x){
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub('^"|"$', '', x)     # remove surrounding quotes if present
  x
}

cat("BASE:", base, "\n\n")

# 1) Read BIOM (OTU matrix only)
cat("1) Reading BIOM:\n")
cat("BIOM path:", biom_fp, "\n")
biom_obj <- read_biom(biom_fp)
cat("biom object class:", class(biom_obj), "\n")
otu_mat_raw <- as(biom_data(biom_obj), "matrix")
cat("Raw OTU matrix dims (rows, cols):", dim(otu_mat_raw), "\n")
cat("Rowname (taxa) sample (first 6):\n"); print(head(rownames(otu_mat_raw), 6))
cat("Colname (samples) sample (first 6):\n"); print(head(colnames(otu_mat_raw), 6))

# sanitize OTU ids and sample ids in the matrix
if(!is.null(rownames(otu_mat_raw))) rownames(otu_mat_raw) <- clean_id(rownames(otu_mat_raw))
if(!is.null(colnames(otu_mat_raw))) colnames(otu_mat_raw) <- clean_id(colnames(otu_mat_raw))

cat("\nAfter cleaning names:\n")
cat("Taxa (first 6):\n"); print(head(rownames(otu_mat_raw), 6))
cat("Samples (first 6):\n"); print(head(colnames(otu_mat_raw), 6))

# 2) Read metadata
cat("\n2) Reading metadata:\n")
cat("Metadata path:", meta_fp, "\n")
meta_dt <- fread(meta_fp, na.strings = c("", "NA"))
cat("Metadata dims:", dim(meta_dt), "\n")
cat("Metadata columns:\n"); print(colnames(meta_dt))
cat("Metadata head (first 6 rows):\n"); print(head(meta_dt, 6))

# detect sample ID column (try common names)
candidate_cols <- c("#SampleID","SampleID","sample-id","Sample ID","sampleID","ID","id","sample")
found <- intersect(candidate_cols, colnames(meta_dt))
if(length(found) > 0){
  sid_col <- found[1]
  cat("Detected sample-ID column by name:", sid_col, "\n")
} else {
  # fall back to detect column whose values match sample names in biom
  match_cols <- sapply(colnames(meta_dt), function(col) {
    sum(clean_id(as.character(meta_dt[[col]])) %in% colnames(otu_mat_raw))
  })
  if(max(match_cols, na.rm=TRUE) > 0){
    sid_col <- names(match_cols)[which.max(match_cols)]
    cat("No standard column name found; using column", sid_col, "because it shares values with OTU sample IDs.\n")
  } else {
    sid_col <- colnames(meta_dt)[1]
    cat("No obvious sample-ID column found; defaulting to first column:", sid_col, "\n")
  }
}

# set rownames and convert to data.frame
meta_df <- as.data.frame(meta_dt)
rownames(meta_df) <- clean_id(as.character(meta_df[[sid_col]]))
cat("After setting rownames, metadata rownames (first 8):\n"); print(head(rownames(meta_df), 8))

# also clean any character columns inside metadata
meta_df[] <- lapply(meta_df, function(col) if(is.character(col)) clean_id(col) else col)

# 3) Determine OTU orientation: are samples columns or rows?
cat("\n3) Checking orientation / matching sample IDs between OTU and metadata\n")
sample_ids_meta <- rownames(meta_df)
matches_col <- sum(sample_ids_meta %in% colnames(otu_mat_raw))
matches_row <- sum(sample_ids_meta %in% rownames(otu_mat_raw))
cat("Metadata sample IDs matching OTU colnames:", matches_col, "\n")
cat("Metadata sample IDs matching OTU rownames:", matches_row, "\n")

if(matches_col >= matches_row){
  cat("Assuming samples are columns in OTU matrix (taxa are rows). No transpose.\n")
  otu_mat <- otu_mat_raw
} else {
  cat("Samples appear to be rows in OTU matrix. Transposing so taxa are rows.\n")
  otu_mat <- t(otu_mat_raw)
}
cat("OTU matrix dims after orientation (taxa rows):", dim(otu_mat), "\n")
cat("OTU sample names (first 6):\n"); print(head(colnames(otu_mat), 6))
cat("OTU taxa names (first 6):\n"); print(head(rownames(otu_mat), 6))

# 4) Read taxonomy file
cat("\n4) Reading taxonomy file:\n")
cat("Taxonomy path:", tax_fp, "\n")
tax_dt <- fread(tax_fp, na.strings = c("", "NA"))
cat("Taxonomy dims:", dim(tax_dt), "\n")
cat("Taxonomy columns:\n"); print(colnames(tax_dt))
cat("Taxonomy head (first 6 rows):\n"); print(head(tax_dt, 6))

# detect feature id and taxonomy columns
feat_candidates <- c("Feature ID","FeatureID","Feature.Id","Feature.ID","OTU ID","OTU_ID","#OTU ID","feature_id","Feature")
feat_col <- intersect(feat_candidates, colnames(tax_dt))
if(length(feat_col) > 0) feat_col <- feat_col[1] else feat_col <- colnames(tax_dt)[1]
cat("Using feature-id column:", feat_col, "\n")

tax_candidates <- c("Taxon","taxonomy","Taxonomy","consensus","Consensus","Classification")
tax_col <- intersect(tax_candidates, colnames(tax_dt))
if(length(tax_col) > 0) tax_col <- tax_col[1] else tax_col <- NA

if(!is.na(tax_col)){
  cat("Found taxonomy string column:", tax_col, "\n")
  tax_strings <- as.character(tax_dt[[tax_col]])
  # split on semicolon and clean
  tax_split <- strsplit(tax_strings, ";", fixed = TRUE)
  clean_token <- function(x){
    x <- gsub("^\\s+|\\s+$", "", x)
    x <- gsub("^[a-zA-Z]__+", "", x)   # remove k__ p__ etc
    if(nchar(x) == 0) x <- NA
    x
  }
  tax_split <- lapply(tax_split, function(v) sapply(v, clean_token, USE.NAMES = FALSE))
  max_ranks <- max(lengths(tax_split), na.rm = TRUE)
  cat("Detected max ranks from taxonomy column:", max_ranks, "\n")
  tax_mat <- do.call(rbind, lapply(tax_split, function(x){ length(x) <- max_ranks; x }))
  rownames(tax_mat) <- clean_id(as.character(tax_dt[[feat_col]]))
  colnames(tax_mat) <- paste0("Rank", seq_len(ncol(tax_mat)))
} else {
  cat("No single taxonomy string column found -> assuming the tax file already has multiple rank columns.\n")
  rank_cols <- setdiff(colnames(tax_dt), feat_col)
  cat("Using rank columns:", paste(rank_cols, collapse = ", "), "\n")
  tax_mat <- as.matrix(tax_dt[, ..rank_cols])
  rownames(tax_mat) <- clean_id(as.character(tax_dt[[feat_col]]))
  colnames(tax_mat) <- make.names(colnames(tax_mat))
}
cat("Tax matrix dims:", dim(tax_mat), "\n")
cat("Tax matrix head (first 6 rows):\n"); print(head(tax_mat, 6))

# 5) Align OTU taxa IDs with taxonomy IDs
cat("\n5) Aligning taxa (OTU IDs) between OTU matrix and taxonomy\n")
otus_in_otu <- rownames(otu_mat)
otus_in_tax <- rownames(tax_mat)
cat("OTUs in OTU table:", length(otus_in_otu), "\n")
cat("OTUs in tax table:", length(otus_in_tax), "\n")
cat("Intersection size:", length(intersect(otus_in_otu, otus_in_tax)), "\n")
cat("First 10 OTUs only in OTU table:\n"); print(head(setdiff(otus_in_otu, otus_in_tax), 10))
cat("First 10 OTUs only in tax table:\n"); print(head(setdiff(otus_in_tax, otus_in_otu), 10))

common_otus <- intersect(otus_in_otu, otus_in_tax)
if(length(common_otus) == 0) stop("ERROR: No overlapping OTU IDs between OTU table and taxonomy table. Check IDs.")
# subset to common
otu_mat2 <- otu_mat[common_otus, , drop = FALSE]
tax_mat2 <- tax_mat[common_otus, , drop = FALSE]
cat("After subsetting to common OTUs - OTU dims:", dim(otu_mat2), " TAX dims:", dim(tax_mat2), "\n")

# 6) Align samples between OTU and metadata
cat("\n6) Aligning samples between OTU and metadata\n")
sample_ids_otu <- colnames(otu_mat2)
sample_ids_meta <- rownames(meta_df)
cat("Samples in OTU:", length(sample_ids_otu), "\n")
cat("Samples in metadata:", length(sample_ids_meta), "\n")
cat("Intersection samples:", length(intersect(sample_ids_otu, sample_ids_meta)), "\n")
cat("First 10 samples only in OTU:\n"); print(head(setdiff(sample_ids_otu, sample_ids_meta), 10))
cat("First 10 samples only in metadata:\n"); print(head(setdiff(sample_ids_meta, sample_ids_otu), 10))

common_samples <- intersect(sample_ids_otu, sample_ids_meta)
if(length(common_samples) == 0) stop("ERROR: No overlapping sample IDs between OTU table and metadata. Check sample naming.")
otu_mat3 <- otu_mat2[, common_samples, drop = FALSE]
meta_df2 <- meta_df[common_samples, , drop = FALSE]
cat("After subsetting to common samples - OTU dims:", dim(otu_mat3), " METADATA dims:", dim(meta_df2), "\n")

# 7) Build phyloseq object and print final diagnostics
cat("\n7) Building phyloseq object\n")
otu_phy <- otu_table(otu_mat3, taxa_are_rows = TRUE)
tax_phy <- tax_table(as.matrix(tax_mat2))
samp_phy <- sample_data(meta_df2)

cat("Sample names in OTU (first 6):\n"); print(head(sample_names(otu_phy), 6))
cat("Sample names in metadata (first 6):\n"); print(head(rownames(meta_df2), 6))
cat("Taxa names (first 6):\n"); print(head(taxa_names(otu_phy), 6))
cat("Tax table columns:", colnames(tax_phy), "\n")

ps <- tryCatch({
  phyloseq(otu_phy, tax_phy, samp_phy)
}, error = function(e){
  cat("ERROR when creating phyloseq object:\n")
  cat("Message:", e$message, "\n")
  cat("Length sample_names(otu_phy):", length(sample_names(otu_phy)), "\n")
  cat("Length rownames(meta_df2):", length(rownames(meta_df2)), "\n")
  cat("First 10 sample_names(otu_phy):\n"); print(head(sample_names(otu_phy),10))
  cat("First 10 rownames(meta_df2):\n"); print(head(rownames(meta_df2),10))
  stop(e)
})

cat("phyloseq created successfully!\n")
cat("ps summary:\n"); print(ps)
cat("Number of samples:", nsamples(ps), "\n")
cat("Number of taxa:", ntaxa(ps), "\n")
cat("Tax table ranks:", colnames(tax_table(ps)), "\n")
cat("Sample variables:", sample_variables(ps), "\n")

# 8) Print warnings (if any)
cat("\n8) Any warnings (first few):\n")
w <- warnings()
print(w)
library(tidytree)
library(treeio)

# Convert taxonomy data to hierarchical tree
tree <- as.phylo(~ Rank1/Rank2/Rank3/Rank4/Rank5/Rank6/Rank7, data = tax)

# Save Newick
write.tree(tree, "graphlan_tree.txt")

library(phyloseq)
library(ALDEx2)
library(data.table)
library(ape)

# --- Extract data from phyloseq ---
counts <- as.data.frame(otu_table(ps))
conds <- sample_data(ps)$Group
tax <- as.data.frame(tax_table(ps))

# --- Filter rare taxa ---
min_prev <- 0.2
taxa_keep <- rowSums(counts > 0)/ncol(counts) >= min_prev
counts <- counts[taxa_keep, ]
tax <- tax[taxa_keep, ]
colnames(tax)

# --- Define comparisons ---
comparisons <- list(
  comp1 = c("E","EH","EL","C","CH","CL"),      # multi-group
  comp2 = c("E","EH","EL"),                    # multi-group
  comp3 = c("C","CH","CL"),                    # multi-group
  comp4 = list(group1 = c("C","CH","CL"), group2 = c("E","EH","EL"))  # two-group
)
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
    
    x <- aldex.clr(counts_sub, conds_sub, mc.samples=128)
    res <- aldex.ttest(x)
    eff <- aldex.effect(x)
    res <- data.frame(res, eff)
    
    sig <- res[res$we.eBH < 0.05, ]
    
    if (nrow(sig) > 0) {
      sig$Size <- abs(sig$effect)
      sig$Color <- ifelse(sig$effect > 0, "red", "blue")
    }
    
  } else {
    # Multi-group comparison
    keep <- conds %in% comp
    counts_sub <- counts[, keep]
    conds_sub <- conds[keep]
    
    x <- aldex.clr(counts_sub, conds_sub, mc.samples=128)
    res <- aldex.kw(x)
    sig <- res[res$kw.eBH < 0.05, ]
    
    if (nrow(sig) > 0) {
      sig$Size <- 1        # placeholder for GraPhlAn
      sig$Color <- "grey"
    }
  }
  
  if (nrow(sig) > 0) {
    # Merge with taxonomy
    sig$taxon <- rownames(sig)
    sig <- merge(sig, tax, by.x="taxon", by.y="row.names", all.x=TRUE)
    
    # Create Node path
    sig$Node <- apply(sig[, c("Rank1","Rank2","Rank3","Rank4","Rank5","Rank6","Rank7")],
                      1, function(x) paste(na.omit(x), collapse="|"))
    
    anno <- sig[, c("Node","Color","Size")]
    
    # Save files
    fwrite(anno, paste0("graphlan_annotation_", i, ".txt"), sep="\t", quote=F, row.names=F)
    fwrite(sig, paste0("sig_taxa_", i, ".csv"), row.names=F)
  } else {
    cat("No significant taxa for", i, "\n")
  }
}


library(ape)

library(tidytree)

# # Convert tax table into a hierarchical tree
# paths <- apply(tax[, c("Rank1","Rank2","Rank3","Rank4","Rank5","Rank6","Rank7")],
#                1, function(x) paste(na.omit(x), collapse="|"))
# tree <- as.phylo(strsplit(paths, "\\|"))

# write.tree(tree, "graphlan_tree.txt")