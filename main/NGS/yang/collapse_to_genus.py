#!/usr/bin/env python3
import os
import logging
import sys
from collections import defaultdict
from typing import Dict, Set, Tuple
from Bio import SeqIO
from ete3 import Tree
from tqdm import tqdm

from ete3 import TreeNode
from ete3 import NodeStyle # Needed for the .collapse() method

# Increase recursion limit for potentially deep phylogenetic trees
sys.setrecursionlimit(10000)

# ----------------- User-editable paths -----------------
# Full FASTA file used for taxonomy annotation
SILVA_FASTA = "/home/patwuch/projects/microbiome/reference/SILVA_full.fasta" 
# Full Newick tree file
FULL_TREE_FILE = "/home/patwuch/projects/microbiome/reference/SILVA_full.ntree" 

# Output file paths for the single tree
OUTPUT_TREE_FILE = "full_genus_collapsed.ntree" 
OUTPUT_ANNOTATION_FILE = "full_genus_annotation.txt" 
# -------------------------------------------------------

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# --- Helper Functions (Retained) ---

def parse_taxonomy_string(t: str) -> dict:
    out = {'order':'','family':'','genus':''}
    if not t: return out
    parts = [p.strip() for p in t.split(';') if p.strip()]
    for p in parts:
        if p.startswith('o__'): out['order'] = p[3:].lower()
        elif p.startswith('f__'): out['family'] = p[3:].lower()
        elif p.startswith('g__'): out['genus'] = p[3:].lower()
    if not out['order'] and len(parts)>3: out['order'] = parts[3].lower()
    if not out['family'] and len(parts)>4: out['family'] = parts[4].lower()
    if not out['genus'] and len(parts)>5: out['genus'] = parts[5].lower()
    return out

def index_fasta(fasta_path: str) -> Dict[str, dict]:
    """Indexes FASTA to map sequence IDs to parsed taxonomic ranks."""
    acc_to_ranks = {}
    for rec in SeqIO.parse(fasta_path,'fasta'):
        full_id = rec.id
        stripped_id = full_id.split('.')[0]
        desc = rec.description
        try:
            tax_str = desc.split(' ',1)[1].strip() if ' ' in desc else desc[len(full_id):].strip()
        except IndexError:
            tax_str = desc[len(full_id):].strip()

        ranks = parse_taxonomy_string(tax_str)
        acc_to_ranks[full_id] = ranks
        if stripped_id != full_id: acc_to_ranks[stripped_id] = ranks
            
    logger.info("Indexed %d sequences from FASTA", len(acc_to_ranks))
    return acc_to_ranks

def normalize_leaf_name(name: str) -> str:
    """Normalizes leaf names (e.g., removing quotes and version numbers)."""
    name = name.strip("'\"")
    return name.split('.')[0]

def create_genus_collapsed_tree(tree_file: str, acc_to_ranks: dict, output_tree_file: str):
    """
    Collapse a full phylogenetic tree to genus-level safely.
    Handles polyphyletic genera by splitting into connected clades,
    normalizes leaf names, and supports both GUI and headless ETE3 builds.
    """
    logger.info("Loading full tree: %s", tree_file)
    try:
        tree = Tree(tree_file, format=1)
    except Exception as e:
        logger.error("Could not load tree from %s: %s", tree_file, e)
        return

    # --- Normalize tree leaf names ---
    normalized_tree_leaves = {normalize_leaf_name(leaf.name): leaf for leaf in tree.iter_leaves()}

    # --- Build genus -> normalized leaves mapping ---
    genus_to_leaves = defaultdict(set)
    for leaf_name, leaf_node in normalized_tree_leaves.items():
        ranks = acc_to_ranks.get(leaf_node.name) or acc_to_ranks.get(leaf_name) or {}
        genus = ranks.get("genus", "").lower()
        if genus:
            genus_to_leaves[genus].add(leaf_name)

    total_genera = len(genus_to_leaves)
    logger.info("Found %d genera in tree.", total_genera)
    if not genus_to_leaves:
        logger.warning("No genera found; aborting.")
        return

    used_leaves = set()
    genus_nodes_count = 0
    polyphyletic_count = 0

    # --- helper: collapse a connected clade safely ---
    def collapse_connected_clade(node, genus_name, size, index=None):
        name = genus_name.capitalize() if not index else f"{genus_name.capitalize()}_{index}"
        node.name = name
        node.add_feature("genus_name", name)
        node.add_feature("genus_size", size)
        if hasattr(node, "collapse") and callable(getattr(node, "collapse")):
            node.collapse()
        else:
            for child in list(node.children):
                node.remove_child(child)
        return name

    # --- main genus loop ---
    for genus, leaves in genus_to_leaves.items():
        remaining = [leaf for leaf in leaves if leaf in normalized_tree_leaves and leaf not in used_leaves]
        if not remaining:
            continue

        genus_capitalized = genus.capitalize()

        # --- single-leaf genus ---
        if len(remaining) == 1:
            leaf_node = normalized_tree_leaves[remaining[0]]
            leaf_node.name = genus_capitalized
            leaf_node.add_feature("genus_name", genus_capitalized)
            leaf_node.add_feature("genus_size", 1)
            used_leaves.update(remaining)
            genus_nodes_count += 1
            continue

        # --- multi-leaf genus (handle possible polyphyly) ---
        leaf_nodes = [normalized_tree_leaves[n] for n in remaining]
        unprocessed = leaf_nodes[:]
        cluster_index = 1

        while unprocessed:
            group = [unprocessed.pop(0)]
            added = True
            while added:
                added = False
                for node in unprocessed[:]:
                    try:
                        names = [n.name for n in group] + [node.name]
                        _ = tree.get_common_ancestor(*names)
                        group.append(node)
                        unprocessed.remove(node)
                        added = True
                    except Exception:
                        continue

            try:
                mrca = tree.get_common_ancestor(*group)
            except Exception as e:
                logger.error("Failed to get MRCA for genus '%s' subgroup: %s", genus_capitalized, e)
                continue
            collapse_connected_clade(mrca, genus, len(group),
                                     index=(cluster_index if len(unprocessed) else None))
            used_leaves.update([normalize_leaf_name(x.name) for x in group])
            genus_nodes_count += 1
            if cluster_index > 1:
                polyphyletic_count += 1
            cluster_index += 1

    logger.info(
        "Tree collapsed. Genus nodes processed: %d (polyphyletic genera split into %d extra groups).",
        genus_nodes_count,
        polyphyletic_count,
    )

    # --- write final tree ---
    try:
        tree.write(outfile=output_tree_file, format=1)
        logger.info("Genus-collapsed tree saved: %s", output_tree_file)
    except Exception as e:
        logger.error("Failed to write tree to %s: %s", output_tree_file, e)







def main():
    # Output file path for the single tree
    OUTPUT_TREE_FILE = "full_genus_collapsed.ntree" 

    acc_to_ranks = index_fasta(SILVA_FASTA)

    logger.info("--- Starting Single Full Tree Collapse ---")
    
    # Call the core function once
    create_genus_collapsed_tree(
        FULL_TREE_FILE, acc_to_ranks, 
        OUTPUT_TREE_FILE
    )

    logger.info("--- Processing complete. ---")


if __name__ == "__main__":
    main()