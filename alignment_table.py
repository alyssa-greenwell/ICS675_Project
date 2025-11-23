#!/usr/bin/env python
from aligner import Aligner
import argparse
import csv
import numpy as np
from multiprocessing import Pool

NUM_WORKERS = 32
ALIGNER = None
SEQUENCES = None

def init_worker(aligner, sequences):
    global ALIGNER, SEQUENCES
    ALIGNER = aligner
    SEQUENCES = sequences

def compute_pair(pair):
    i, j = pair
    seq1 = SEQUENCES[i]
    seq2 = SEQUENCES[j]
    score = ALIGNER.global_align(seq1, seq2)
    return (i, j, score)

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
        path = "cat_mitochondria_genomes/" + cat + ".fasta"
        sequence = read_fasta(path)
        cat_sequences.append(sequence)
    return cat_sequences
    
def make_alignment_table(cats, cats_len):
    '''Creates a table of alignment scores of each feline compared to every
    other feline.'''
    cat_sequences = get_fastas(cats)
    alignment_matrix = [[0] * cats_len for _ in range(cats_len)]
    score_map = {'A':0, 'C':1, 'G':2, 'T':3}
    # Convert score matrix to array
    scoring_matrix = np.array([
        [0, 1, 1, 1],   # A
        [1, 0, 1, 1],   # C
        [1, 1, 0, 1],   # G
        [1, 1, 1, 0],   # T
    ], dtype=np.int32)
    
    gap_penalty = 1
    my_aligner = Aligner(scoring_matrix, gap_penalty, score_map)
 
    # Create list of unique pairs
    pairs = [(i, j) for i in range(cats_len) for j in range(i, cats_len)]

    with Pool(
            processes=NUM_WORKERS,
            initializer=init_worker,
            initargs=(my_aligner, cat_sequences)
            ) as pool:
        for i, j, score in pool.imap_unordered(compute_pair, pairs):
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
        action="store_true",
        help="Activates demo mode"
    )
    
    args = parser.parse_args()
    
    if (args.demo):
        NUM_WORKERS = 6
        print(f"alignment_table running in Demo Mode. Num workers = {NUM_WORKERS}")
        cats = ["asiatic_golden_cat", "bobcat", "cheetah"]
        cats_len = 3
        output_file = "DEMO_alignment_table.csv"
    else:
        NUM_WORKERS = 32
        cats = ["asiatic_golden_cat", "bobcat", "cheetah", "domestic_cat", 
                "eurasian_lynx", "european_wildcat", "jaguar", "leopard", 
                "lion", "mainland_clouded_leopard", "mainland_leopard_cat", 
                "marbled_cat", "pallas_s_cat", "puma", "sand_cat", 
                "sunda_clouded_leopard", "tiger"]
        cats_len = 17
        output_file = "alignment_table.csv"
    
    alignment_table = make_alignment_table(cats, cats_len)
    
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(alignment_table)
