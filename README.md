# ICS675_Project

Attempts to recreate a simplified feline phylogenetic tree using the mitochondria genomes of 17 feline species.

## Making the Alignment Table
Run `python3 alignment_table.py` to make the alignment table for all feline genomes. This outputs a file called [alignment_table.csv](output/alignment_table.csv).

To run in demo mode and generate a 3x3 demo matrix using just the asiatic golden cat, bobcat, and cheetah, run `python3 alignment_table.py --demo`. This will generate output into a file called DEMO_alignment_table.py.

## Making the Phylogenetic Tree
To generate the phylogenetic tree, run `python3 phylogenetic_tree.py`. This will use the UPGMA algorithm found in [phylo_clustering.py](phylo_clustering.py) to define all the nodes of the tree, then use [tree_plotter.py](tree_plotter.py) to generate a visualization of the tree. This will produce two output files, one called [upgma_output.txt](output/upgma_output.txt) which lists all internal nodes and edge lengths in the tree, and [UPGMA_plot.png](output/UPGMA_plot.png), which is an image of the phylogenetic tree.

## Declaration of Generative AI Usage
Generative AI was used in this project to modify the [tree_plotter.py](tree_plotter.py) script to use midpoint rooting instead of oupgroup rooting. It was also used to help add multithreading to the [alignment_table.py](alignment_table.py) script to reduce the runtime.
