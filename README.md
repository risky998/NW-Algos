# NW-linear
A python program to produce a sequence alignment for two given sequences based on linear gap penalties.
Based on the course content of CS4775 (Computational Genomics) taught at Cornell University.

The program runs a python script to align two given gene sequences using a linear gap penalty
model.

The script takes in 3 arguments:

    f - FASTA file with sequences in FASTA format. 
    s - JSON with the score matrix for alignment. 
    d - The gap penalty for the alignment. 

The script print the alignment to the console

Example Usage: 

    python nwlinear.py -f sequences.fasta -s score_matrix.json -d 100

Replace "seqences.fasta" with any file containing sequences to be aligned in the FASTA format.

Replace score-matrix.json with a custom scoring matrix.

A sample sequences file and score matrix are included in this repository.
