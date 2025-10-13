#!/usr/bin/env python3
import os
import glob
import logging
from collections import defaultdict
from typing import Dict, Set, Tuple
from Bio import SeqIO
from ete3 import Tree
import pandas as pd

# ----------------- User-editable paths -----------------
SIG_TAXA_DIR = "/home/patwuch/projects/microbiome/experiments/楊/ALDEx2_results"
SILVA_FASTA = "/home/patwuch/projects/microbiome/reference/SILVA_pruned.fasta"
FULL_TREE_FILE = "/home/patwuch/projects/microbiome/reference/SILVA_pruned.ntree"
OUTPUT_DIR = "pruned_trees"
os.makedirs(OUTPUT_DIR, exist_ok=True)
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

def parse_csv_taxa(csv_path: str) -> list:
    df = pd.read_csv(csv_path)
    if 'Taxon' not in df.columns:
        raise ValueError("CSV must contain a 'Taxon' column")
    taxa_list = [parse_taxonomy_string(str(t)) for t in df['Taxon']]
    logger.info("Parsed %d taxa from %s", len(taxa_list), os.path.basename(csv_path))
    return taxa_list

def index_fasta(fasta_path: str) -> Tuple[Dict[str, dict], Dict[Tuple[str,str], Set[str]]]:
    acc_to_ranks = {}
    rank_index = defaultdict(set)
    for rec in SeqIO.parse(fasta_path,'fasta'):
        full_id = rec.id
        stripped_id = full_id.split('.')[0]
        desc = rec.description
        tax_str = desc.split(' ',1)[1].strip() if ' ' in desc else desc[len(full_id):].strip()
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

def normalize_leaf_name(name: str) -> str:
    name = name.strip("'\"")
    return name.split('.')[0]

def annotate_tree(tree_file: str, sig_taxa_list: list, acc_to_ranks: Dict[str, dict],
                  output_tree_file: str, output_annotation_file: str):
    """
    Annotate the tree with significant taxa without pruning.
    Leaves matching taxa in sig_taxa_list are marked for Graphlan.
    """
    tree = Tree(tree_file, format=1)
    annotations = {}  # leaf_name -> {'order','family','genus','highlight'}

    # Create a set of accessions to highlight
    highlight_accs = set()
    for tax in sig_taxa_list:
        genus = tax.get('genus','').lower()
        family = tax.get('family','').lower()
        order = tax.get('order','').lower()
        for level,val in [('genus',genus),('family',family),('order',order)]:
            if val:
                for acc,ranks in acc_to_ranks.items():
                    if ranks.get(level,'').lower() == val:
                        highlight_accs.add(acc)
                break

    # Annotate tree leaves
    for leaf in tree.iter_leaves():
        acc_key = normalize_leaf_name(leaf.name)
        ranks = acc_to_ranks.get(leaf.name) or acc_to_ranks.get(acc_key) or {}
        genus = ranks.get('genus','')
        family = ranks.get('family','')
        order = ranks.get('order','')
        highlight = 'yes' if leaf.name in highlight_accs or acc_key in highlight_accs else 'no'
        annotations[leaf.name] = {'order':order,'family':family,'genus':genus,'highlight':highlight}

    # Write annotated tree (full tree, unpruned)
    tree.write(outfile=output_tree_file, format=5)
    logger.info("Annotated tree saved: %s", output_tree_file)

    # Write annotation table
    with open(output_annotation_file,'w') as f:
        f.write("Leaf\tOrder\tFamily\tGenus\tHighlight\n")
        for leaf_name,tax in annotations.items():
            f.write(f"{leaf_name}\t{tax['order']}\t{tax['family']}\t{tax['genus']}\t{tax['highlight']}\n")
    logger.info("Annotation file saved: %s", output_annotation_file)

def main():
    acc_to_ranks, rank_index = index_fasta(SILVA_FASTA)

    for csv_file in glob.glob(os.path.join(SIG_TAXA_DIR,"sig_taxa_*.csv")):
        comp_name = os.path.basename(csv_file).replace("sig_taxa_","").replace(".csv","")
        taxa_list = parse_csv_taxa(csv_file)
        output_tree_file = os.path.join(OUTPUT_DIR,f"{comp_name}_annotated.tree")
        output_annotation_file = os.path.join(OUTPUT_DIR,f"{comp_name}_annotation.txt")
        annotate_tree(FULL_TREE_FILE, taxa_list, acc_to_ranks, output_tree_file, output_annotation_file)

if __name__ == "__main__":
    main()
