#!/usr/bin/env python

"""
Project_Name: main/filters/sequence, File_name: sequence_similarity.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 14/04/15 , Time: 04:52 PM

Globally align the amino acid sequences in smotifs against the target sequence
"""

def orderSeq(previous_smotif, current_seq, direction):
    """

    :param previous_smotif:
    :param current_seq:
    :param direction:
    :return:
    """
    previous_seq = ''

    for entry in previous_smotif:
        if entry[0] == 'seq_filter':
            seq_filter = entry
            previous_seq = seq_filter[1]

    if direction == 'left':
        concat_seq = current_seq + previous_seq
    else:
        concat_seq = previous_seq + current_seq

    return concat_seq

def getSmotifAASeq(ss1, ss2):
    one_letter = {'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q',
                  'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y',
                  'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A',
                  'GLY': 'G', 'PRO': 'P', 'CYS': 'C', 'ASX': 'D', 'GLX': 'G', 'UNK': 'A'}
    seq = ''
    # TODO try to optimize this for loop
    for entry in ss1:
        if entry[2] == 'CA':
            aa = entry[1]
            seq = seq + one_letter[aa]
    for entry in ss2:
        if entry[2] == 'CA':
            aa = entry[1]
            seq = seq + one_letter[aa]

    return seq

def getSmotifAASeq_v2(sse):

    # TODO delete this or above def

    one_letter = {'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q',
                  'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y',
                  'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A',
                  'GLY': 'G', 'PRO': 'P', 'CYS': 'C', 'ASX': 'D', 'GLX': 'G', 'UNK': 'A'}
    seq = ''

    for entry in sse:
        if entry[2] == 'CA':
            aa = entry[1]
            seq = seq + one_letter[aa]
    return seq


def SequenceSimilarity(s1_def, s2_def, smotif, exp_data):
    """
    return sequence identity for given unique seqs and
    new queried sequences
    """
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist

    matrix = matlist.blosum62
    gap_open = -10
    gap_extend = -0.5

    aa_seq = exp_data['aa_seq']

    # get the target and smotif seq information alone and exclude the loop regions
    native_seq = aa_seq[s1_def[4] - 1:s1_def[5]] + aa_seq[s2_def[4] - 1:s2_def[5]]  # -1 to fix residue numbering
    smotif_seq = getSmotifAASeq(smotif[1], smotif[2])

    # Perform the alignment
    alns = pairwise2.align.globalds(native_seq, smotif_seq, matrix, gap_open, gap_extend)
    # the top alignment is in the first entry of the array
    top_aln = alns[0]

    # Compute teh sequence identity
    seqa, qseqa, score, begin, end = top_aln
    j, k = 0.0, 0.0
    for i in range(0, len(qseqa)):
        if qseqa[i] != '-' and seqa[i] != '-':
            j += 1
            if qseqa[i] == seqa[i]:
                k += 1
    # seq_id = (k/j)*100
    seq_id = (k / len(smotif_seq)) * 100
    return smotif_seq, seq_id, score

def S2SequenceSimilarity(ss_def, smotif, direction, exp_data):
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

    aa_seq = exp_data['aa_seq']

    # change below such that previous and current sses are chosen for the native seq
    # NO, only current sse , based on the direction should be chosen and only one SSE seq is aligned
    native_seq = aa_seq[ss_def[4] - 1:ss_def[5]]

    if direction == 'left':
        smotif_sse = smotif[1]
    else:
        smotif_sse = smotif[2]

    smotif_seq = getSmotifAASeq_v2(smotif_sse)

    # Perform alignment
    alns = pairwise2.align.globalds(native_seq, smotif_seq, matrix, gap_open, gap_extend)
    # the best alignment is in the first entry
    top_aln = alns[0]

    # Calculate the sequence identity
    seqa, qseqa, score, begin, end = top_aln
    j, k = 0.0, 0.0
    for i in range(0, len(qseqa)):
        if qseqa[i] != '-' and seqa[i] != '-':
            j += 1
            if qseqa[i] == seqa[i]:
                k += 1
    # seq_id = (k/j)*100
    seq_id = (k / len(smotif_seq)) * 100
    return smotif_seq, seq_id, score

def getS1SeqIdentity(s1_def, s2_def, smotif, exp_data):
    aa_seq = exp_data['aa_seq']
    # get the target and smotif seq information alone and exclude the loop regions
    native_seq = aa_seq[s1_def[4] - 1:s1_def[5]] + aa_seq[s2_def[4] - 1:s2_def[5]]  # -1 to fix residue numbering
    smotif_seq = getSmotifAASeq(smotif[1], smotif[2])
    k = 0.0
    for i in range(0, len(native_seq)):
        if native_seq[i] == smotif_seq[i]:
            k += 1.0
    seq_id = k/float(len(native_seq))
    return smotif_seq, seq_id

def getSXSeqIdentity(ss_def, smotif, direction, exp_data, psmotif, sse_ordered):
    aa_seq = exp_data['aa_seq']

    # change below such that previous and current sses are chosen for the native seq
    # NO, only current sse , based on the direction should be chosen and only one SSE seq is aligned

    if direction == 'left':
        smotif_sse = smotif[1]
    else:
        smotif_sse = smotif[2]

    smotif_seq = getSmotifAASeq_v2(smotif_sse)

    concat_seq = orderSeq(psmotif, smotif_seq, direction)

    native_sse_seq = ''
    for sse in sse_ordered:
        sse_seq = aa_seq[sse[4] - 1: sse[5]]
        native_sse_seq = native_sse_seq + sse_seq
    k = 0.0
    if len(concat_seq) != len(native_sse_seq):
        print "Something is wrong with extracting sequence information"
    for i in range(0, len(concat_seq)):
        # print native_sse_seq[i], concat_seq[i]
        if native_sse_seq[i] == concat_seq[i]:
            k += 1
    seq_id = (k / float(len(concat_seq))) * 100
    return seq_id, concat_seq

def getGlobalSequenceIdentity(concat_seq, exp_data, sse_ordered):
    aa_seq = exp_data['aa_seq']
    # [['helix', 12, 3, 6, 4, 15], ['helix', 11, 4, 10, 24, 34], ['helix', 7, 10, 8, 45, 51]]
    native_sse_seq = ''
    for sse in sse_ordered:
        sse_seq = aa_seq[sse[4] - 1: sse[5]]
        native_sse_seq = native_sse_seq + sse_seq
    k=0.0
    for i in range(0, len(concat_seq)):
        # print native_sse_seq[i], concat_seq[i]
        if native_sse_seq[i] == concat_seq[i]:
            k += 1
    seq_id = (k / float(len(concat_seq))) * 100
    return seq_id