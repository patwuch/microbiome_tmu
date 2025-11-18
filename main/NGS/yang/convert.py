#!/usr/bin/env python3
from ete3 import Tree
from Bio import Phylo
import sys
import glob
import os
import io 
import threading # Required for getting/setting the recursion limit

# ----------------------------------------------------
# ⚠️ FIX: INCREASE RECURSION LIMIT ⚠️
# The default is usually 1000. Increase it for deep trees.

# You can also use a fixed large number, like:
sys.setrecursionlimit(10000) 
# ----------------------------------------------------

in_dir = "/home/patwuch/projects/microbiome/experiments/yang/pruned_trees"
out_dir = "/home/patwuch/projects/microbiome/experiments/yang/converted_trees"
os.makedirs(out_dir, exist_ok=True)

# ... (rest of the script from the previous successful conversion) ...
for infile in glob.glob(os.path.join(in_dir, "*.tree")):
    outfile = os.path.join(out_dir, os.path.basename(infile).replace(".tree", ".xml"))
    
    print(f"Converting {infile} → {outfile}")
    
    try:
        # 1. Load the ETE Tree (will now use the higher recursion limit)
        t = Tree(infile, format=1)
        
        # 2. Export the ETE Tree to a Newick string
        newick_string = t.write(format=1) 
        
        # 3. Use Biopython's Phylo module to read the Newick string
        handle = io.StringIO(newick_string)
        biopython_tree = Phylo.read(handle, "newick")

        # 4. Use Biopython to write the tree to the final PhyloXML file
        Phylo.write(biopython_tree, outfile, "phyloxml")
        
        print(f"Successfully exported to {outfile}")
        
    except Exception as e:
        print(f"ERROR converting {infile}: {e}", file=sys.stderr)

print("Conversion batch completed.")