from Bio import SeqIO
from ete3 import Tree
import pandas as pd
import gzip

# --- Paths ---
sig_taxa_csv = "/home/patwuch/projects/microbiome/experiments/楊/sig_taxa_comp4.csv"
silva_tree_gz = "/home/patwuch/projects/microbiome/reference/SILVA_138.2_SSURef_NR99.accessions.ntree.gz"
silva_fasta = "/home/patwuch/projects/microbiome/reference/SILVA_138.2_SSURef_NR99_tax_silva.fasta"
output_tree = "/home/patwuch/projects/microbiome/experiments/楊/sig_taxa_comp4_genus_fuzzy.tree"
annotation_file = "/home/patwuch/projects/microbiome/experiments/楊/gp_annotations.txt"

# --- Load CSV ---
sig_df = pd.read_csv(sig_taxa_csv)

# Parse CSV: order, family, genus
csv_taxa = []
for taxon_str in sig_df['Taxon']:
    levels = taxon_str.split(";")
    order = levels[3].split("__")[-1] if len(levels) > 3 else ""
    family = levels[4].split("__")[-1] if len(levels) > 4 else ""
    genus = levels[5].split("__")[-1] if len(levels) > 5 else ""
    csv_taxa.append({'order': order, 'family': family, 'genus': genus})

# --- Load SILVA FASTA and build family+genus → accessions mapping ---
fg_to_acc = {}
acc_to_ranks = {}
for record in SeqIO.parse(silva_fasta, "fasta"):
    ranks = record.description.split(None,1)[1].split(";")
    order = ranks[3] if len(ranks) > 3 else ""
    family = ranks[4] if len(ranks) > 4 else ""
    genus = ranks[5] if len(ranks) > 5 else ""
    acc_to_ranks[record.id] = {'order': order, 'family': family, 'genus': genus}
    key = (order, family, genus)
    fg_to_acc.setdefault(key, []).append(record.id)

# --- Efficient matching: CSV taxa → SILVA accessions ---
keep_accessions = set()
for tax in csv_taxa:
    order, family, genus = tax['order'], tax['family'], tax['genus']
    
    # Exact genus match if available
    if genus:
        for key in fg_to_acc:
            if key[0] == order and key[1] == family and genus in key[2]:
                keep_accessions.update(fg_to_acc[key])
    else:
        # Keep all sequences at family level if genus missing
        for key in fg_to_acc:
            if key[0] == order and key[1] == family:
                keep_accessions.update(fg_to_acc[key])

keep_accessions = list(keep_accessions)
print(f"Total SILVA accessions retained: {len(keep_accessions)}")

# --- Load and prune SILVA tree ---
with gzip.open(silva_tree_gz, 'rt') as f:
    tree = Tree(f.read(), format=1)

tree_leaves = set(tree.get_leaf_names())
valid_accessions = [acc for acc in keep_accessions if acc in tree_leaves]

# --- Relabel leaves with genus/family for GraPhlAn ---
leaf_label_map = {}
for leaf in tree:
    if leaf.name not in acc_to_ranks:
        continue
    ranks = acc_to_ranks[leaf.name]
    genus = ranks['genus']
    family = ranks['family']
    order = ranks['order']
    
    # Prefer genus, fallback to family if missing
    new_label = genus if genus else family if family else leaf.name
    
    # Ensure uniqueness
    if new_label in leaf_label_map.values():
        new_label = f"{new_label}_{leaf.name}"
    
    leaf_label_map[leaf.name] = new_label
    leaf.name = new_label
print("Example FASTA IDs:", list(acc_to_ranks.keys())[:5])
print("Example tree leaves:", list(tree.get_leaf_names())[:5])

# --- Prune tree ---
keep_labels = [leaf_label_map[acc] for acc in valid_accessions if acc in leaf_label_map]
if not keep_labels:
    raise ValueError("No matching leaves found for pruning!")
tree.prune(keep_labels, preserve_branch_length=True)
tree.write(outfile=output_tree, format=1)
print(f"Pruned tree saved: {output_tree} (Leaves: {len(keep_labels)})")

# --- Generate GraPhlAn annotation file ---
with open(annotation_file, "w") as f:
    f.write("#FORMAT: annotation\n")
    for leaf_name, tax in zip(keep_labels, valid_accessions):
        # Example: color by significance (can be adapted)
        sig_color = "red"  # You could map this from your CSV values
        f.write(f"{leaf_name}\tclade\tcolor\t{sig_color}\n")

print(f"GraPhlAn annotation file saved: {annotation_file}")
