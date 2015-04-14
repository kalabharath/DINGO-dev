#!/usr/bin/env python

"""
Project_Name: main/filters/sequence, File_name: sequence_similarity.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 14/04/15 , Time: 04:52 PM

Globally align the amino acid sequences in smotifs against the target sequence
"""




def SequenceSimilarity(s1_def, s2_def, smotif, threshold):
    """
    return sequence identity for given unique seqs and
    new queried sequences
    """
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
    matrix = matlist.blosum62
    gap_open = -10
    gap_extend = -0.5
    hit = True
    for seq in seq_array:
        alns = pairwise2.align.globalds(seq, qseq, matrix, gap_open, gap_extend)
        top_aln = alns[0]
        seqa, qseqa, score, begin, end = top_aln
        j, k = 0.0, 0.0
        for i in range (0,len(qseqa)):
            if qseqa[i] != '-' and seqa[i] != '-':
                j +=1
            if qseqa[i] == seqa[i]:
                k +=1
        seq_id = (k/j)*100
        #print seq_id
        if seq_id > cutoff:
            return False
    return hit