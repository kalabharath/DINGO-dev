#!/usr/bin/env python

"""
Project_Name: main/filters/sequence, File_name: sequence_similarity.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 14/04/15 , Time: 04:52 PM

Globally align the amino acid sequences in smotifs against the target sequence
"""

import utility.io_util as io


def getSmotifAASeq(ss1, ss2):
    one_letter = {'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q', \
                  'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y', \
                  'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A', \
                  'GLY': 'G', 'PRO': 'P', 'CYS': 'C', 'ASX': 'A', 'GLX': 'G'}
    seq = ''

    for entry in ss1:
        if entry[2] == 'CA':
            aa = entry[1]
            seq = seq + one_letter[aa]
    for entry in ss2:
        if entry[2] == 'CA':
            aa = entry[1]
            seq = seq + one_letter[aa]

    return seq


def SequenceSimilarity(s1_def, s2_def, smotif, threshold):
    """
    return sequence identity for given unique seqs and
    new queried sequences
    """
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist

    matrix = matlist.blosum62
    # matrix = matlist.pam120
    gap_open = -10
    gap_extend = -0.5
    hit = True

    header, aa_seq = io.readAASeq()
    # print s1_def
    # print s2_def
    # print aa_seq
    native_seq = aa_seq[s1_def[4] - 1:s1_def[5]] + aa_seq[s2_def[4] - 1:s2_def[5]]  # -1 to fix residue numbering

    smotif_def = smotif[0][0]
    smotif_ss1 = smotif[0][1]
    smotif_ss2 = smotif[0][2]
    smotif_seq = getSmotifAASeq(smotif_ss1, smotif_ss2)
    # print smotif_def
    # print native_seq, len(native_seq)
    # print smotif_seq, len(smotif_seq)


    alns = pairwise2.align.globalds(native_seq, smotif_seq, matrix, gap_open, gap_extend)
    top_aln = alns[0]

    seqa, qseqa, score, begin, end = top_aln
    j, k = 0.0, 0.0
    for i in range(0, len(qseqa)):
        if qseqa[i] != '-' and seqa[i] != '-':
            j += 1
            if qseqa[i] == seqa[i]:
                k += 1
    # seq_id = (k/j)*100
    seq_id = (k / len(smotif_seq)) * 100

    if seq_id >= threshold and score > 0:
        return seq_id, score, True
    else:
        return seq_id, score, False



