from tkinter import messagebox
from enum import IntEnum
import numpy as np
import OutputScreen
def get_PAM_Matrix():
    try:
        data = open("PAM.txt", 'r')
        dimensions = data.readline()
        try:
            n = int(dimensions)
        except ValueError:
            print("Wrong file")
            exit()
        letters = data.readline()
        letters = letters.replace('\n', '')
        letters_arr = letters.split('\t')
        score_matrix = np.zeros((n, n))
        for i in range(0, n):
            arr = data.readline().split("\t")
            for j in range(0, n):
                score_matrix[i][j] = float(arr[j])
        return  score_matrix
    except FileNotFoundError:
        print("There is no matrix file.")
        exit()

Protein_Letters={"A":0,"R":1,"N":2	,"D":3,"C":4,"Q":5,"E":6,"G":7,"H":8,"I":9,"L":10,"K":11,"M":12,"F":13,"P":14,"S":15,"T":16,"W":17,"Y":18,"V":19}


class Trace(IntEnum):
   STOP = 0
   LEFT = 1
   UP = 2
   DIAGONAL = 3




def check_input(seq):
    valid=True
    if (seq == ""):
        return False
    for i in seq:
        if (i not in Protein_Letters.keys()):
            return False
    return True


def Local_Alignmenet_Protein(seq1, seq2 ,gap):
    align_loc = []
    """Do a local alignment between x and y"""
    # create a zero-filled matrix
    # Generating the empty matrices for storing scores and tracing
    row = len(seq1) + 1
    col = len(seq2) + 1
    PAM_matrix=get_PAM_Matrix()
    matrix = np.zeros((row, col))
    tracing_matrix = np.zeros((row, col))

    # Initialising the variables to find the highest scoring cell
    max_score = -1
    max_index = (-1, -1)

    # Calculating the scores for all cells in the matrix
    for i in range(1, row):
        for j in range(1, col):
            # Calculating the diagonal score (match score)
            diagonal_score = matrix[i - 1, j - 1] + PAM_matrix[Protein_Letters[seq1[i-1]],Protein_Letters[seq2[j-1]]]

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
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = seq2[max_j - 1]
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

def get_argments_Protein(seq1,seq2,gap,window):
   TOP_SEQ=check_input(seq1)
   BOTTOM_SEQ=check_input(seq2)
   if(TOP_SEQ==True &BOTTOM_SEQ==True):
       window.destroy()
       Local_Alignmenet_Protein(seq1,seq2,gap)
   else:
       messagebox.showerror("ERROR","Please Enter A Correct Protein Sequence Letters")






