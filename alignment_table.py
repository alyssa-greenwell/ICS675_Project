#!/usr/bin/env python
from aligner import Aligner
import argparse
import csv

def read_fasta(path):
    '''Reads a single fasta file and returns the sequence as a string.'''
    seq_lines = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_lines.append(line)

    return "".join(seq_lines)

def get_fastas(cats):
    '''Gets all feline fasta genomes and puts their sequences into an array.
    returns the array of sequences.'''
    cat_sequences = []
    for cat in cats:
        path = "cat_mitochondria_genomes/" + cat
        sequence = read_fasta(path)
        cat_sequences.append(sequence)
    
def make_alignment_table(cats, cats_len):
    '''Creates a table of alignment scores of each feline compared to every
    other feline.'''
    cat_sequences = get_fastas(cats)
    alignment_matrix = [[0] * cats_len for _ in range(cats_len)]
    scoring_matrix = {"A": {"A": 2, "T": -1, "C": -1, "G": -1},
                      "T": {"A": -1, "T": 2, "C": -1, "G": -1},
                      "C": {"A": -1, "T": -1, "C": 2, "G": -1},
                      "G": {"A": -1, "T": -1, "C": -1, "G": 2}}
    
    gap_penalty = -2
    my_aligner = Aligner(scoring_matrix, gap_penalty)
    for i in range(cats_len):
        for j in range(cats_len):
           score = my_aligner.global_align(cat_sequences[i], cat_sequences[j]) 
           alignment_matrix[i][j] = score
           alignment_matrix[j][i] = score
           print(f"Finished alignment for {cats[i]} and {cats[j]}")
    return alignment_matrix


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Set the mode to run the file in"
    )
    
    parser.add_argument(
        "--demo",
        help="Activates demo mode"
    )
    
    args = parser.parse_args()
    
    if (args.demo):
        print("alinment_table running in Demo Mode.")
        cats = ["asiatic_golden_cat", "bobcat", "cheetah"]
        cats_len = 3
        output_file = "DEMO_alignment_table"
    else:
        cats = ["asiatic_golden_cat", "bobcat", "cheetah", "domestic_cat", 
                "eurasian lynx", "european_wildcat", "jaguar", "leopard", 
                "lion", "mainland_clouded_leopard", "mainland_leopard_cat", 
                "marbled_cat", "pallas_s_cat", "puma", "sand_cat", 
                "sunda_clouded_leopard", "tiger"]
        cats_len = 17
        output_file = "alignment_table"
    
    alignment_table = make_alignment_table(cats, cats_len)
    
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(alignment_table)
