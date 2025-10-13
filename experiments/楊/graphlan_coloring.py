# -*- coding: utf-8 -*-
"""
Generate a GraPhlAn annotation file with:
- Leaf colors for significant taxa
- Ring (circular heatmap) for significant taxa
- Optional subtree shading and labels
"""

from ete3 import Tree
import pandas as pd

# -------------------------
# INPUT FILES
# -------------------------
pruned_tree_file = "/home/patwuch/projects/microbiome/experiments/楊/sig_taxa_genus.tree"
sig_taxa_csv = "/home/patwuch/projects/microbiome/experiments/楊/sig_taxa_comp4.csv"
annotation_file = "/home/patwuch/projects/microbiome/experiments/楊/gp_annotations.txt"

# -------------------------
# Load pruned tree
# -------------------------
tree = Tree(pruned_tree_file, format=1)

# -------------------------
# Load significant taxa from CSV
# -------------------------
sig_df = pd.read_csv(sig_taxa_csv)
sig_taxa = set()
for taxon_str in sig_df['Taxon']:
    levels = taxon_str.split(";")
    genus = levels[5].split("__")[-1].strip() if len(levels) > 5 else ""
    if genus:
        sig_taxa.add(genus)

# -------------------------
# Annotation generator
# -------------------------
def write_graphlan_annotations(tree, output_file, sig_taxa,
                               leaf_color="red", ring_level=1,
                               ring_width=0.1, ring_height=0.1,
                               ring_alpha=1.0):
    """
    Generates a fully compliant GraPhlAn annotation file.
    """
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("#FORMAT: annotation\n")

        # ---------------------
        # Global ring settings
        # ---------------------
        f.write(f"ring_label\t{ring_level}\tSignificant taxa\n")
        f.write(f"ring_label_font_size\t{ring_level}\t10\n")
        f.write(f"ring_internal_separator_thickness\t{ring_level}\t0.05\n")
        f.write(f"ring_external_separator_thickness\t{ring_level}\t0.05\n")
        f.write(f"ring_separator_color\t{ring_level}\tk\n")  # black

        # ---------------------
        # Leaf coloring and ring entries
        # ---------------------
        # sig_taxa = set of genera from CSV
        # leaf_label_map = leaf_name in tree -> genus (lowercase)
        for leaf in tree.get_leaf_names():
            genus = leaf_label_map.get(leaf, "").lower()
            if genus in sig_taxa:
                # Leaf color
                f.write(f"{leaf}\tcolor\tred\n")
                # Ring entries
                f.write(f"{leaf}ring_color\t1\tred\n")
                f.write(f"{leaf}ring_width\t1\t0.1\n")
                f.write(f"{leaf}ring_height\t1\t0.1\n")
                f.write(f"{leaf}ring_alpha\t1\t1.0\n")
                # Subtree annotation
                f.write(f"[{leaf}]\tannotation_background_color\t#FFCCCC\n")
                f.write(f"[{leaf}]\tannotation\t{genus}\n")
                f.write(f"[{leaf}]\tannotation_font_size\t8\n")


# -------------------------
# Run annotation generation
# -------------------------
write_graphlan_annotations(tree, annotation_file, sig_taxa)
print(f"GraPhlAn annotation file written: {annotation_file}")
