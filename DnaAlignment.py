from tkinter import *
from tkinter import messagebox
import numpy as np
from enum import IntEnum
import OutputScreen
DNA_nucleotide=["A","C","G","T"]

class Trace(IntEnum):
   STOP = 0
   LEFT = 1
   UP = 2
   DIAGONAL = 3


def check_input(seq):
    valid=True
    if(seq==""):
        return False
    for i in seq:
        if (i not in DNA_nucleotide):
            return False
    return True


def Local_Alignmenet_DNA(seq1, seq2, match, mismatch, gap):
    align_loc = []
    """Do a local alignment between x and y"""
    # create a zero-filled matrix
    # Generating the empty matrices for storing scores and tracing
    row = len(seq2) + 1
    col = len(seq1) + 1
    matrix = np.zeros((row, col))
    tracing_matrix = np.zeros((row, col))

    # Initialising the variables to find the highest scoring cell
    max_score = -1
    max_index = (-1, -1)

    # Calculating the scores for all cells in the matrix
    for i in range(1, row):
        for j in range(1, col):
            # Calculating the diagonal score (match score)
            match_value = match if seq1[j - 1] == seq2[i - 1] else mismatch
            diagonal_score = matrix[i - 1, j - 1] + match_value

            # Calculating the vertical gap score
            vertical_score = matrix[i - 1, j] + gap

            # Calculating the horizontal gap score
            horizontal_score = matrix[i, j - 1] + gap

            # Taking the highest score
            matrix[i, j] = max(0, diagonal_score, vertical_score, horizontal_score)

            # Tracking where the cell's value is coming from
            if matrix[i, j] == 0:
                tracing_matrix[i, j] = Trace.STOP

            elif matrix[i, j] == horizontal_score:
                tracing_matrix[i, j] = Trace.LEFT

            elif matrix[i, j] == vertical_score:
                tracing_matrix[i, j] = Trace.UP

            elif matrix[i, j] == diagonal_score:
                tracing_matrix[i, j] = Trace.DIAGONAL

                # Tracking the cell with the maximum score
            if matrix[i, j] >= max_score:
                max_index = (i, j)
                max_score = matrix[i, j]

    # Initialising the variables for tracing
    aligned_seq1 = ""
    aligned_seq2 = ""
    current_aligned_seq1 = ""
    current_aligned_seq2 = ""
    (max_i, max_j) = max_index
    Score=matrix[max_i, max_j]

    # Tracing and computing the pathway with the local alignment
    while tracing_matrix[max_i, max_j] != Trace.STOP:
        if tracing_matrix[max_i, max_j] == Trace.DIAGONAL:
            current_aligned_seq1 = seq1[max_j - 1]
            current_aligned_seq2 = seq2[max_i - 1]
            max_i = max_i - 1
            max_j = max_j - 1
            align_loc.append((max_i + 1, max_j + 1))

        elif tracing_matrix[max_i, max_j] == Trace.UP:
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = '-'
            max_i = max_i - 1
            align_loc.append((max_i + 1, max_j + 1))

        elif tracing_matrix[max_i, max_j] == Trace.LEFT:
            current_aligned_seq1 = '-'
            current_aligned_seq2 = seq2[max_j - 1]
            max_j = max_j - 1
            align_loc.append((max_i + 1, max_j + 1))


        aligned_seq1 = aligned_seq1 + current_aligned_seq1
        aligned_seq2 = aligned_seq2 + current_aligned_seq2
    align_loc.append((max_i, max_j ))
    # Reversing the order of the sequences
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]
    OutputScreen.output_information(matrix,align_loc,aligned_seq1,aligned_seq2,seq1,seq2,Score)

def get_argments_DNA(seq1,seq2,match,mismatch,gap,window):
   TOP_SEQ=check_input(seq1.upper())
   BOTTOM_SEQ=check_input(seq2.upper())
   if(TOP_SEQ==True &BOTTOM_SEQ==True):
       window.destroy()
       Local_Alignmenet_DNA(seq1,seq2,match,mismatch,gap)
   else:
       messagebox.showerror("ERROR","Please Enter A Correct DNA Sequence Letters")






