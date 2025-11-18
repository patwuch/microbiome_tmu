#!/usr/bin/env python3
"""
Prune a SILVA tree and FASTA so that they match exactly.
Handles typical SILVA ID differences (version numbers, taxonomic annotations, coordinates in NTREE).
"""

from ete3 import Tree
from Bio import SeqIO
import re

# --------- USER CONFIG ---------
tree_file = "/home/patwuch/projects/microbiome/reference/SILVA_138.2_SSURef_NR99.accessions.ntree"
fasta_file = "/home/patwuch/projects/microbiome/reference/SILVA_138.2_SSURef_NR99_tax_silva.fasta"
pruned_tree_file = "SILVA_pruned.ntree"
pruned_fasta_file = "SILVA_pruned.fasta"
# -------------------------------

def normalize_fasta_id(fasta_id):
    """Keep only the accession prefix (remove version and suffix)."""
    return fasta_id.split('.')[0]

def normalize_tree_leaf(leaf_name):
    """Extract accession ID from tree leaf name like 'CP010125.4624908.4626464'."""
    # Remove quotes if present
    leaf_name = leaf_name.strip("'\"")
    # Take the part before the first dot (accession)
    return leaf_name.split('.')[0]

# 1. Load tree
print("Loading tree...")
tree = Tree(tree_file, format=1)

# 2. Load FASTA IDs
print("Loading FASTA IDs...")
fasta_ids = set(normalize_fasta_id(rec.id) for rec in SeqIO.parse(fasta_file, "fasta"))

# 3. Prune tree
print("Pruning tree to match FASTA...")
for leaf in tree.get_leaves():
    if normalize_tree_leaf(leaf.name) not in fasta_ids:
        leaf.delete()

print(f"Tree pruned. Leaves remaining: {len(tree)}")

# 4. Prune FASTA to match pruned tree
pruned_tree_leaves = set(normalize_tree_leaf(leaf.name) for leaf in tree.get_leaves())
pruned_records = [rec for rec in SeqIO.parse(fasta_file, "fasta")
                  if normalize_fasta_id(rec.id) in pruned_tree_leaves]

print(f"FASTAs pruned. Sequences remaining: {len(pruned_records)}")

# 5. Save outputs
print("Saving pruned tree and FASTA...")
tree.write(outfile=pruned_tree_file)
SeqIO.write(pruned_records, pruned_fasta_file, "fasta")

print("Done!")
