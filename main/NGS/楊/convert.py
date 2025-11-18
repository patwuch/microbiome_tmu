#!/usr/bin/env python3
import os
import glob
import logging
import sys # Needed for stderr and recursion
import io # Needed for in-memory string I/O
from collections import defaultdict
from typing import Dict, Set, Tuple
from Bio import SeqIO, Phylo # Phylo added for PhyloXML conversion
from ete3 import Tree
from ete3 import NodeStyle, TextFace # NodeStyle and TextFace needed for tree annotation
import pandas as pd

# Increase recursion limit for potentially deep phylogenetic trees
sys.setrecursionlimit(2000)

# ----------------- User-editable paths -----------------
SIG_TAXA_DIR = "/home/patwuch/projects/microbiome/experiments/楊/ALDEx2_results"
SILVA_FASTA = "/home/patwuch/projects/microbiome/reference/SILVA_full.fasta"
FULL_TREE_FILE = "/home/patwuch/projects/microbiome/reference/SILVA_full.ntree"

# Output directories for the two-step process
PRUNED_TREES_DIR = "pruned_trees"
CONVERTED_TREES_DIR = "converted_trees"

os.makedirs(PRUNED_TREES_DIR, exist_ok=True)
os.makedirs(CONVERTED_TREES_DIR, exist_ok=True)
# -------------------------------------------------------

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

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

def parse_csv_taxa(csv_path: str) -> set:
    df = pd.read_csv(csv_path)
    if 'Taxon' not in df.columns:
        raise ValueError("CSV must contain a 'Taxon' column")
    
    sig_genera = set()
    for t in df['Taxon']:
        ranks = parse_taxonomy_string(str(t))
        if ranks['genus']:
            sig_genera.add(ranks['genus'])
    logger.info("Parsed %d unique significant genera from %s", len(sig_genera), os.path.basename(csv_path))
    return sig_genera

def index_fasta(fasta_path: str) -> Tuple[Dict[str, dict], Dict[Tuple[str,str], Set[str]]]:
    acc_to_ranks = {}
    rank_index = defaultdict(set)
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
        
        for level in ('order','family','genus'):
            val = ranks.get(level,'')
            if val:
                rank_index[(level,val)].add(full_id)
                if stripped_id != full_id:
                    rank_index[(level,val)].add(stripped_id)
                    
    logger.info("Indexed %d sequences from FASTA", len(acc_to_ranks))
    return acc_to_ranks, rank_index

def get_genus_to_leaves(tree: Tree, acc_to_ranks: Dict[str, dict]) -> Dict[str, Set[str]]:
    """Maps each genus name to the set of leaf names belonging to that genus."""
    genus_map = defaultdict(set)
    for leaf in tree.iter_leaves():
        acc_key = normalize_leaf_name(leaf.name)
        ranks = acc_to_ranks.get(leaf.name) or acc_to_ranks.get(acc_key) or {}
        genus = ranks.get('genus', '').lower()
        if genus:
            genus_map[genus].add(leaf.name)
    return genus_map

def normalize_leaf_name(name: str) -> str:
    name = name.strip("'\"")
    return name.split('.')[0]

def prune_and_collapse_tree(tree_file: str, sig_genera: Set[str], acc_to_ranks: Dict[str, dict],
                           comp_name: str, output_tree_file: str, output_annotation_file: str):
    """
    Prunes the tree to only include leaves that have a genus assignment,
    then collapses all members of the same genus into a single node.
    """
    logger.info("Loading full tree for collapsing...")
    tree = Tree(tree_file, format=1)
    
    genus_to_leaves = get_genus_to_leaves(tree, acc_to_ranks)
    all_leaves_with_genus = set().union(*genus_to_leaves.values())
    
    # 1. Prune the tree: Keep only leaves that have a genus annotation
    tree.prune(all_leaves_with_genus)
    logger.info("Tree pruned to %d leaves with genus annotation.", len(tree.get_leaves()))
     
    # 2. Collect all MRCA nodes to be collapsed
    mrca_nodes_to_collapse = []
    genus_to_mrca = {}
    
    for genus, leaves in genus_to_leaves.items():
        existing_leaves = [leaf for leaf in leaves if tree.search_nodes(name=leaf)]
        
        if not existing_leaves:
            continue
        
        # Find the MRCA for all members of the genus
        mrca = tree.get_common_ancestor(existing_leaves)
        
        # We only collapse nodes that are *not* already terminal (i.e., contain more than one leaf)
        if len(mrca.get_leaves()) > 1:
            mrca_nodes_to_collapse.append(mrca)
            
        mrca.add_face(ete3.TextFace(genus.capitalize(), fsize=14), column=0, position='branch-top') # Example ETE3 annotation
        mrca.name = genus.capitalize()
        mrca.add_feature("is_significant", genus.lower() in sig_genera)
        mrca.add_feature("genus_size", len(existing_leaves))
        genus_to_mrca[genus] = mrca

    # 3. Permanently collapse the nodes (this significantly reduces the leaf count)
    # NOTE: You must import 'NodeStyle' from ete3 for this to work
    for node in mrca_nodes_to_collapse:
        node.collapse()
    
    logger.info("Tree collapsed. New number of leaves (genus nodes): %d", len(tree.get_leaves()))
            
    # 3. Final Write of the modified Newick tree (format=1 is plain Newick)
    tree.write(outfile=output_tree_file, format=1) 
    logger.info("Genus-collapsed tree saved: %s", output_tree_file)

    # 4. Prepare Annotation File (Graphlan format)
    with open(output_annotation_file, 'w') as f:
        f.write("# Graphlan Annotation for Genus-Collapsed Tree\n")
        
        for genus, node in collapsed_nodes.items():
            node_name = node.name
            
            f.write(f"{node_name}\tclade_label\t{node_name}\n")
            
            if node.is_significant:
                f.write(f"{node_name}\tclade_label_color\t#FF0000\n") 
                f.write(f"{node_name}\tclade_marker_color\t#FF0000\n")
                f.write(f"{node_name}\tbranch_color\t#FF0000\n")
            else:
                f.write(f"{node_name}\tbranch_color\t#A9A9A9\n") 

            f.write(f"{node_name}\text_annotations\t({node.genus_size} ASVs)\n")

    logger.info("Annotation file (Graphlan format) saved: %s", output_annotation_file)


def main():
    acc_to_ranks, rank_index = index_fasta(SILVA_FASTA)

    # --- Step 1: Prune, Collapse, and Annotate Trees ---
    logger.info("--- Starting Tree Pruning and Collapse ---")
    tree_files_to_convert = []
    for csv_file in glob.glob(os.path.join(SIG_TAXA_DIR,"sig_taxa_*.csv")):
        comp_name = os.path.basename(csv_file).replace("sig_taxa_","").replace(".csv","")
        sig_genera = parse_csv_taxa(csv_file) 
        
        output_tree_file = os.path.join(PRUNED_TREES_DIR,f"{comp_name}_genus_collapsed.ntree")
        output_annotation_file = os.path.join(PRUNED_TREES_DIR,f"{comp_name}_genus_annotation.txt")
        
        prune_and_collapse_tree(
            FULL_TREE_FILE, sig_genera, acc_to_ranks, comp_name,
            output_tree_file, output_annotation_file
        )
        tree_files_to_convert.append(output_tree_file)

    # -------------------------------------------------------------------
    # --- Step 2: Convert Newick (.ntree) to PhyloXML (.xml) format ---
    # -------------------------------------------------------------------
    logger.info("--- Starting PhyloXML Conversion ---")

    for infile in tree_files_to_convert:
        # NOTE: Using CONVERTED_TREES_DIR as the output directory
        outfile = os.path.join(CONVERTED_TREES_DIR, os.path.basename(infile).replace(".ntree", ".xml"))
        
        print(f"Converting {infile} → {outfile}")
        
        try:
            # 1. Load the ETE Tree (using the Newick file output from Step 1)
            # We use format=1 since the prune_and_collapse_tree function explicitly wrote Newick format.
            t = Tree(infile, format=1) 
            
            # 2. Export the ETE Tree to a Newick string
            # ETE's internal structure may contain more features than standard Newick.
            newick_string = t.write(format=1) 
            
            # 3. Use Biopython's Phylo module to read the Newick string
            handle = io.StringIO(newick_string)
            biopython_tree = Phylo.read(handle, "newick")

            # 4. Use Biopython to write the tree to the final PhyloXML file
            Phylo.write(biopython_tree, outfile, "phyloxml")
            
            print(f"Successfully exported to {outfile}")
            
        except Exception as e:
            print(f"ERROR converting {infile}: {e}", file=sys.stderr)

    logger.info("Conversion batch completed.")


if __name__ == "__main__":
    main()