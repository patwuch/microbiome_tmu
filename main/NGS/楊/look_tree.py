from ete3 import Tree

tree_path = "/home/patwuch/projects/microbiome/experiments/楊/sig_taxa_comp4_genus_fuzzy.tree"

# Try loading with ETE3
try:
    tree = Tree(tree_path, format=1)
    print("Tree loaded successfully!")
except Exception as e:
    print("Tree loading error:", e)

# Print first 5 leaves
print("First 5 leaf names:")
for i, leaf in enumerate(tree.iter_leaves()):
    print(f"{i+1}: {leaf.name}")
    if i >= 4:
        break