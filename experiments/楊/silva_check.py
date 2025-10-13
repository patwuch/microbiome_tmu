from ete3 import Tree
from Bio import SeqIO

# --- Paths ---
pruned_tree_path = "/home/patwuch/projects/microbiome/experiments/楊/sig_taxa_comp4_genus_fuzzy.tree"
silva_fasta = "/home/patwuch/projects/microbiome/reference/SILVA_138.2_SSURef_NR99_tax_silva.fasta"
clean_tree_path = "/home/patwuch/projects/microbiome/experiments/楊/sig_taxa_comp4_genus_fuzzy_labels.tree"

# --- Step 1: Build accession -> genus/family mapping ---
acc_to_label = {}
for record in SeqIO.parse(silva_fasta, "fasta"):
    ranks = record.description.split(None,1)[1].split(";")
    family = ranks[4] if len(ranks) > 4 else ""
    genus = ranks[5] if len(ranks) > 5 else ""
    label = genus if genus != "" else family
    acc_to_label[record.id] = label

# --- Step 2: Load pruned tree ---
tree = Tree(pruned_tree_path, format=1)

# --- Step 3: Map tree tips to genus/family using base accession ---
for leaf in tree:
    leaf_base = leaf.name.split(".")[0]  # take only the first part of tip
    mapped = False
    for fasta_acc, label in acc_to_label.items():
        if fasta_acc.startswith(leaf_base):
            leaf.name = label
            mapped = True
            break
    if not mapped:
        # fallback: keep original accession
        leaf.name = leaf_base

# --- Step 4: Clean tip names ---
for leaf in tree:
    leaf.name = leaf.name.replace(" ", "_")           # replace spaces
    leaf.name = leaf.name.encode('ascii','ignore').decode()  # remove non-ASCII
    for char in ["(", ")", "[", "]", "/", "\\", ",", ";", ":"]:
        leaf.name = leaf.name.replace(char, "_")

# --- Step 5: Remove internal node names ---
for node in tree.traverse():
    if not node.is_leaf():
        node.name = ""

# --- Step 6: Write clean GraPhlAn-ready tree ---
tree.write(outfile=clean_tree_path, format=1)
print(f"GraPhlAn-ready tree written to: {clean_tree_path}")

# Optional: inspect first 10 tip labels
print("First 10 tip labels:")
for i, leaf in enumerate(tree.iter_leaves()):
    print(f"{i+1}: {leaf.name}")
    if i >= 9:
        break
