#!/usr/bin/env python
import numpy as np

class Aligner():
    """
    Implements the Levenshtein distance algorithm to calculate the distance
    between two sequences.
    
    Input:
    score_matrix: 2-d numpy array with the scoring matrix.
    score_map: maps the IUPAC codes (A, T, C, G) to indexes in the scoring matrix.
    gap_penalty: score for a gap in the alignment.
    seq1: the first sequence to align.
    seq2: the second sequence to align.
    
    Outputs the alignment score as an integer.
    """
    def __init__(self, score_matrix, gap_penalty, score_map):
        self.gap_penalty = gap_penalty
        self.score_matrix = score_matrix
        self.score_map = score_map
        return
    
    def global_align(self, seq1, seq2):
        # Initialize global matrix
        x_len = len(seq1) + 1 # matrix width
        y_len = len(seq2) + 1 # matrix height
        global_matrix = np.zeros((y_len, x_len), dtype=np.int32)

        # Create base cases for global matrix
        global_matrix[0, :] = np.arange(x_len) * self.gap_penalty
        global_matrix[:, 0] = np.arange(y_len) * self.gap_penalty
        # Fill in global matrix
        for i in range(1, y_len):
            seq2_char = seq2[i-1]
            for j in range(1, x_len):
                seq1_char = seq1[j-1]
                try:
                    dia_score = self.score_matrix[self.score_map[seq2_char]][self.score_map[seq1_char]]
                except Exception as e:
                    dia_score = 1
                dia = global_matrix[i-1][j-1] + dia_score
                hor = global_matrix[i][j-1] + self.gap_penalty
                ver = global_matrix[i-1][j] + self.gap_penalty
                new_cell = min(hor, ver, dia)
                global_matrix[i, j] = new_cell

        #Traceback
        alignment_score = global_matrix[y_len-1][x_len-1]
       
        return alignment_score
    