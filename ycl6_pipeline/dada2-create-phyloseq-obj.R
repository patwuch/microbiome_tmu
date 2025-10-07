#

library("dada2")
library("data.table")
library("DECIPHER")
library("phangorn")
library("phyloseq")
library("ggplot2")
options(width=190)

# Get script location when running with Rscript
args <- commandArgs(trailingOnly = FALSE)
script.path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
script.dir <- dirname(script.path)

# Define raw, trimmed, filt, processed paths accordingly
fastq <- file.path(script.dir, "..", "data", "raw", "20250905")
processed.dir <- file.path(script.dir, "..", "data", "processed", "20250905")
ycl6.dir <- file.path(processed.dir, "ycl6")
trimmed <- file.path(ycl6.dir, "trimmed")
filt    <- file.path(ycl6.dir, "filt")

# Continued from DADA2 part 1 (dada2-per-run-processing.R or cutadapt-and-dada2-per-run-processing.R)

# Load saved workspace
load(file.path(ycl6.dir, "image1.RData"))

# Load sequence table from multiple runs if applicable
s01 = readRDS(file.path(ycl6.dir, "seqtab.rds"))
st.all = s01 # here we define st.all as s01 since there is only one run in this example
# s02 = readRDS("seqtab2.rds")
# st.all = mergeSequenceTables(s01, s02)

# If this script is used independently and the "sample.names" is not carried over from previous step/script
sample.names = gsub(".1.fastq.gz", "", rownames(st.all))


# To sum values of same sample from multiple sequence table (i.e. when a sample was re-sequenced due to low depth)
# st.all = mergeSequenceTables(s01, s02, repeats = "sum")

# Remove bimeras (two-parent chimeras)
seqtab.nochim = removeBimeraDenovo(st.all, method = 'consensus', verbose = TRUE, multithread = TRUE)
rownames(seqtab.nochim) = sample.names

dim(seqtab)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track reads through the pipeline
# If you load saved workspace from previous step/script, you will have the necessary objects to build the track data.frame
getN <- function(x) sum(getUniques(x))
track = cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) = c("Trimmed", "Filtered", "denoisedF", "denoisedR", "merged", "nonchim")
track = cbind(data.frame(SampleID = sample.names), track)
write.table(track, file = file.path(ycl6.dir, "track.txt"), sep = "\t", quote = F, row.names = F, col.names = T)

# Assign taxonomy
# Note to change the PATH to reference training datasets accordingly
dbpath = "/home/patwuch/projects/microbiome/reference/"
ref1 = paste0(dbpath, "silva_nr99_v138.2_toSpecies_trainset.fa.gz")
ref2 = paste0(dbpath, "silva_v138.2_assignSpecies.fa.gz")
ref3 = paste0(dbpath, "ncbi_refseq_16S_combined.fna") # Combination of NCBI RefSeq 16S Archaea and Bacteria

# Extracts the sequences
seqs = getSequences(seqtab.nochim)

# Classifies sequences against SILVA reference training dataset
set.seed(12345)
taxtab = assignTaxonomy(seqs, refFasta = ref1, minBoot = 80, tryRC = TRUE, 
	outputBootstraps = TRUE, verbose = TRUE, multithread = TRUE)

# Taxonomic assignment to the species level by exact matching against SILVA and NCBI reference datasets
spec_silva = assignSpecies(seqs, ref2, allowMultiple = FALSE, tryRC = TRUE, verbose = TRUE)
spec_ncbi = assignSpecies(seqs, ref3, allowMultiple = FALSE, tryRC = TRUE, verbose = TRUE)

# Combine species-level taxonomic assignment from 2 reference sources
SVformat = paste("%0",nchar(as.character(ncol(seqtab.nochim))),"d", sep = "")
svid = paste0("SV", sprintf(SVformat, seq(ncol(seqtab.nochim))))

s_silva = as.data.frame(spec_silva, stringsAsFactors = FALSE)
rownames(s_silva) = svid

s_ncbi = as.data.frame(spec_ncbi, stringsAsFactors = FALSE)
rownames(s_ncbi) = svid
s_ncbi$Genus = gsub("\\[|\\]", "", s_ncbi$Genus)

s_merged = cbind(s_ncbi, s_silva)
colnames(s_merged) = c("nGenus","nSpecies","sGenus","sSpecies")
s_merged1 = s_merged[!is.na(s_merged$nSpecies),]	                        # NCBI assignment not empty
colnames(s_merged1)[1:2] = c("Genus","Species")
s_merged2 = s_merged[is.na(s_merged$nSpecies) & !is.na(s_merged$sSpecies),]	# no NCBI assignment but Silva not empty
colnames(s_merged2)[3:4] = c("Genus","Species")
s_merged3 = s_merged[is.na(s_merged$nSpecies) & is.na(s_merged$sSpecies),]	# no NCBI and Silva assignment
colnames(s_merged3)[3:4] = c("Genus","Species")

s_final = rbind(s_merged1[,c("Genus","Species")], s_merged2[,c("Genus","Species")], s_merged3[,c("Genus","Species")])
s_final = s_final[order(row.names(s_final)),]

s_final = as.matrix(s_final)

if("Genus" %in% colnames(taxtab$tax)) {
        gcol = which(colnames(taxtab$tax) == "Genus")
} else { 
	gcol = ncol(taxtab$tax) 
}

matchGenera <- function(gen.tax, gen.binom, split.glyph = "/") {
	if(is.na(gen.tax) || is.na(gen.binom)) { return(FALSE) }
	if((gen.tax == gen.binom) || grepl(paste0("^", gen.binom, "[ _", split.glyph, "]"), gen.tax) || grepl(paste0(split.glyph, gen.binom, "$"), gen.tax)) {
		return(TRUE)
	} else {
		return(FALSE)
	}
}

gen.match = mapply(matchGenera, taxtab$tax[,gcol], s_final[,1])
taxtab$tax = cbind(taxtab$tax, s_final[,2])
colnames(taxtab$tax)[ncol(taxtab$tax)] = "Species"
print(paste(sum(!is.na(s_final[,2])), "out of", nrow(s_final), "were assigned to the species level."))

taxtab$tax[!gen.match,"Species"] = NA
print(paste("Of which", sum(!is.na(taxtab$tax[,"Species"])), "had genera consistent with the input table."))

# Prepare df
df = data.frame(sequence = seqs, abundance = colSums(seqtab.nochim), stringsAsFactors = FALSE)
df$id = svid

df = merge(df, as.data.frame(taxtab), by = "row.names")
rownames(df) = df$id
df = df[order(df$id),2:ncol(df)]

# Construct phylogenetic tree
alignment = AlignSeqs(Biostrings::DNAStringSet(setNames(df$sequence, df$id)), anchor = NA)

# Export alignment
phang.align = phyDat(as(alignment, "matrix"), type = "DNA")
write.phyDat(phang.align, file = "alignment.fasta", format = "fasta")
write.phyDat(phang.align, file = "alignment.aln", format = "phylip")

# Phylogenetic tree construction
# Use system2() to execute commands in R
# Note to change the PATH to raxml and raxmlng accordingly
# Classic RAxML HPC multithreaded
raxml = "/home/patwuch/miniforge3/envs/16S/bin/raxmlHPC-PTHREADS-SSE3"
# RAxML-NG
raxmlng = "/home/patwuch/miniforge3/envs/16S/bin/raxml-ng"

system2(raxml, args = c("-T 2", "-f E", "-p 1234", "-x 5678", "-m GTRCAT", "-N 1", "-s alignment.aln", "-n raxml_tree_GTRCAT", "--redo"))
system2(raxmlng, args = c("--evaluate", "--force", "--seed 1234", "--log progress", "--threads 2", "--msa alignment.fasta", "--model GTR+G", "--tree RAxML_fastTree.raxml_tree_GTRCAT", "--brlen scaled", "--prefix GTRCAT", "--redo"))

# Import tree
# raxml-ng (v0.7.0): GTRCAT.raxml.mlTrees
# raxml-ng (v0.9.0): GTRCAT.raxml.bestTree
raxml_tree = read_tree("GTRCAT.raxml.bestTree")

# Load sample info (sample.meta)
samdf = data.frame(fread(file.path(fastq, "metadata.tsv"), colClasses = "character"))
str(samdf)
head(samdf)
anyDuplicated(samdf$sampleid)
anyNA(samdf$sampleid)


# Set rownames to sample IDs
rownames(samdf) <- samdf$sampleid

# Convert columns to factors with proper ordering
samdf$sampleid <- factor(samdf$sampleid, levels = sort(unique(samdf$sampleid)))
samdf$Group <- as.factor(samdf$Group)
samdf$MainType <- as.factor(samdf$MainType)
samdf$Modifier <- as.factor(samdf$Modifier)

# ----------------------------
# Prepare OTU Table
# ----------------------------
# Clean up seqtab rownames
new_seqtab <- seqtab.nochim
rownames(new_seqtab) <- sub("_S.*", "", rownames(new_seqtab)) 

# Ensure samples are in the same order and only include those with metadata
# Use `intersect` to find the common samples
common_samples <- intersect(rownames(new_seqtab), samdf$sampleid)

new_taxtab <- as.data.frame(taxtab)

length(rownames(new_seqtab))
length(samdf[match(rownames(new_seqtab), samdf$sampleid), "sampleid"])
sum(is.na(match(rownames(new_seqtab), samdf$sampleid)))
setdiff(rownames(new_seqtab), samdf$sampleid)

# Replace column names (taxa) with short IDs from df
colnames(new_seqtab) <- df[match(colnames(new_seqtab), df$sequence), "id"]
# Replace rownames with sample IDs
rownames(new_seqtab) <- samdf[match(rownames(new_seqtab), samdf$sampleid), "sampleid"]
sum(is.na(match(rownames(new_seqtab), samdf$sampleid)))

# ----------------------------
# Prepare Taxonomy Table
# ----------------------------
tax <- as.data.frame(new_taxtab$tax)
# Ensure taxa IDs match OTU table columns
rownames(tax) <- colnames(new_seqtab)

# Convert Family and Genus to character (optional, avoids factor issues)
tax$Family <- as.character(tax$Family)
tax$Genus <- as.character(tax$Genus)

# ----------------------------
# Create Phyloseq Object
# ----------------------------
ps <- phyloseq(
	otu_table(new_seqtab, taxa_are_rows = FALSE),
	tax_table(as.matrix(tax)),
	sample_data(samdf),
	phy_tree(raxml_tree)
)

# Inspect phyloseq object
ps

# ----------------------------
# Save Workspace
# ----------------------------
save.image(file = "image2.RData")

