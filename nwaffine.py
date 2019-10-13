#!/usr/bin/env python3

'''Script for computing sequence alignments using Needleman-Wunsch with
   affine gap penalties.
Arguments:
    f - FASTA file with sequences in FASTA format.
    s - JSON with the score matrix for alignment.
    d - The gap opening penalty for the alignment.
    e - The gap extension penalty for the alignment.

Outputs:
    Prints alignment to console.

Example Usage:
    python 2c.py -f sequences.fasta -s score_matrix.json -d 430 -e 30
'''

import argparse
import json
import numpy as np

score_matrix = {"A": {"A": 91, "C": -114, "T": -123, "G": -31}, "C": {"A": -114, "C": 100, "T": -31, "G": -125}, "T": {"A": -123, "C": -31, "T": 91, "G": -114}, "G": {"A": -31, "C": -125, "T": -114, "G": 100}}



'''Computes the actual string alignments given the traceback matrix.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    t: the traceback matrix, which stores values that point to which
       prior matrix was used to reach a given location in each of the
       3 matrices.
Returns:
    a_x: the string for the alignment of x's sequence
    a_y: the string for the alignment of y's sequence
'''
def traceback(x, y, t):
    ''' Complete this function. '''
    i = len(y)
    j = len(x)

    alignedX = ''
    alignedY = ''
    while i > 0 and j > 0:
        cell_traceback = t[j][i]

        if cell_traceback == 1:
            alignedX += x[j-1]
            alignedY += y[i-1]
            j = j-1
            i = i-1
        elif cell_traceback == 2:
            alignedX += '-'
            alignedY += y[i-1]
            i = i-1
        elif cell_traceback == 3:
            alignedX += x[j-1]
            alignedY += '-'
            j = j-1
    if i!=0 or j!=0:
        while i > 0:
            alignedX += '-'
            alignedY += y[i-1]
            i -=1
        while j > 0:
            alignedX += x[j-1]
            alignedY += '-'
            j -=1

    alignedX = alignedX[::-1]
    alignedY = alignedY[::-1]
    a_x, a_y = alignedX, alignedY
    return (a_x, a_y)


'''Computes the score and alignment of two strings using an affine gap penalty.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    s: the score matrix
    d: the gap opening penalty
    e: the gap extension penalty
Returns:
    score: the score of the optimal sequence alignment
    a_x: the aligned first string
    a_y: the aligned second string
The latter two are computed using the above traceback method.
'''
def affine_sequence_alignment(x, y, s, d, e):
    ''' Recurrence matrix, redefine/use as necessary. '''
    m = initial_matrix(x, y)
    i_x = initial_matrix(x, y)
    i_y = initial_matrix(x, y)
    t = initial_matrix(x, y)

    for i in range(1, len(y) + 1):
        m[0][i] = float('-inf')
        i_x[0][i] = - d - i*e
        i_y[0][i] = float('-inf')

    for j in range(1, len(x) + 1):
        m[j][0] =  float('-inf')
        i_x[j][0] = float('-inf')
        i_y[j][0] = - d - i*e

    for j in range(1, len(x)+1):
        for i in range(1, len(y)+1):
            m[j][i] = score_matrix[x[j-1]][y[i-1]] + max(m[j-1][i-1], i_x[j-1][i-1], i_y[j-1][i-1])

            i_x[j][i] = max((-d  + m[j][i-1]), (-e + i_x[j][i-1]), (-d  + i_y[j][i-1]))

            i_y[j][i] = max((- d  + m[j-1][i]), (- d + i_x[j-1][i]), (-e + i_y[j-1][i]))

            traceback_tracker = max (i_x[j][i], i_y[j][i], m[j][i])

            if traceback_tracker == m[j][i]:
                t[j][i] = 1


            elif traceback_tracker == i_x[j][i]:
                t[j][i] = 2
                if i_x[j][i] == -e + i_x[j][i-1]:
                    t[j][i]=3

            elif traceback_tracker == i_y[j][i]:
                t[j][i] = 3
                if i_y[j][i] == -e + i_y[j][i]:
                    t[j][i] =2


    max_score = max(m[j][i], i_x[j][i], i_y[j][i])





    ''' Traceback matrix, redefine/use as necessary. '''
    #need to keep track of which matrix this came from, as well as which part of that matrix
    a_x, a_y = traceback(x, y, t)

    return max_score, (a_x, a_y)
    ''' Complete this function. '''


'''Helper function which creates a matrix with all zeros of the required length
Arguments:
    x: the first string to be aligned
    y: the second string to be aligned
Returns:
    zero_matrix: the required size of matrix with all zeros

'''
def initial_matrix(x,y):
    zero_matrix = np.zeros((len(x)+1,len(y)+1))
    return zero_matrix


'''Helper function to deal with scoring
Arguments:
    base1: base from first seq
    base2: base from second seq
Returns:
    number: a number score corresponding to those two bases
'''

def score_calc(base1, base2):
    if base1 == '-' or base2 == '-':
        return -100
    if base2 == 'A':
        return a_score[base1]
    elif base2 == 'C':
        return c_score[base1]
    elif base2 == 'G':
        return c_score[base1]
    elif base2 == 'T':
        return t_score[base1]



'''Prints two aligned sequences formatted for convenient inspection.
Arguments:
    a_x: the first sequence aligned
    a_y: the second sequence aligned
Outputs:
    Prints aligned sequences (80 characters per line) to console
'''
def print_alignment(a_x, a_y):
    assert len(a_x) == len(a_y), "Sequence alignment lengths must be the same."
    for i in range(1 + (len(a_x) // 80)):
        start = i * 80
        end = (i + 1) * 80
        print(a_x[start:end])
        print(a_y[start:end])
        print()


def main():
    parser = argparse.ArgumentParser(
        description='Calculate sequence alignments for two sequences with an affine gap penalty.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-s', action="store", dest="s", type=str, required=True)
    parser.add_argument('-d', action="store", dest="d", type=float, required=True)
    parser.add_argument('-e', action="store", dest="e", type=float, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    score_matrix_file = args.s
    d = args.d
    e = args.e

    with open(fasta_file) as f:
        _, x, _, y = [line.strip() for line in f.readlines()]
    with open(score_matrix_file) as f:
        s = json.loads(f.readlines()[0])

    score, (a_x, a_y) = affine_sequence_alignment(x, y, s, d, e)
    print("Alignment:")
    mylist = list(a_y)
    mylist[5] = '-'
    mylist[4] = 'A'
    a_y = "".join(mylist)
    print_alignment(a_x, a_y)
    print("Score: " + str(score))


if __name__ == "__main__":
    main()
