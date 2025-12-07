# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Generates a phylogenetic tree for seventeen felidae species.

Processes the distance matrix from a CSV file, uses the UPGMA class to build a
phylogenetic tree, and runs tree_plotter.py to generate a visualization of the
tree. The taxon labels are hard-coded.
"""
import csv
import subprocess
from phylo_clustering import UPGMA

def read_csv_to_int_array(filename):
    array = []
    with open(filename, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            int_row = [int(x) for x in row]
            array.append(int_row)
    return array

if __name__ == "__main__":
    labels = ["asiatic_golden_cat", "bobcat", "cheetah", "domestic_cat", 
            "eurasian_lynx", "european_wildcat", "jaguar", "leopard", 
            "lion", "mainland_clouded_leopard", "mainland_leopard_cat", 
            "marbled_cat", "pallas_s_cat", "puma", "sand_cat", 
            "sunda_clouded_leopard", "tiger"]
    
    # Create distance matrix from CSV file
    distance_matrix = read_csv_to_int_array("alignment_table.csv")
    
    # Build the UPGMA tree
    tree_builder = UPGMA(distance_matrix, labels)
    tree_builder.run("upgma_output.txt")
    
    # Generate visualization of tree
    subprocess.run([
        "python",
        "tree_plotter.py",
        "upgma_output.txt",
        "-o",
        "UPGMA_plot.png"
    ])
    