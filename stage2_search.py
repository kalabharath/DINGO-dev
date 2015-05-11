#!/usr/bin/env python

"""
Project_Name: main, File_name: stage2_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 7/05/15 , Time:10:05 PM

Perform stage 2 in parallel
"""
import utility.stage2_util as uts2
import utility.smotif_util as sm
import utility.io_util as io
import filters.sequence.sequence_similarity as Sfilter
import filters.contacts.contacts_filter  as Cfilter
import filters.rmsd.qcp as qcp



def getfromDB(searched_smotifs, current_ss, direction):

    psmotif = searched_smotifs[-1]

    if direction == 'left':
        previous_ss = psmotif[0]
    else:
        previous_ss = psmotif[1]

    # current_ss, previous_ss
    if direction == 'left':#double check this implementation
        smotif_def = sm.getSmotif(current_ss, previous_ss)
    else:
        smotif_def = sm.getSmotif(previous_ss, current_ss)

    return sm.readSmotifDatabase(smotif_def)


def orderSSE(previous_seq, current_sse):
    ordered_SSE = []

    for sse_array in previous_seq[1]:
        ordered_SSE.append(sse_array)
    ordered_SSE.append(current_sse)
    return ordered_SSE


def SmotifSearch(index_array):

    #print index_array

    psmotif = uts2.getPreviousSmotif(index_array[0])
    current_ss, direction = uts2.getSS2(index_array[1])


    csmotif_data = getfromDB(psmotif[-1], current_ss, direction)
    exp_data = io.readPickle("exp_data.pickle")
    """
    always narrow down to previous sse and current sse and operate on them individually

    """
    dump_log = []

    sse_ordered = orderSSE(psmotif[-1], current_ss)
    for i in range(0, len(csmotif_data)):
    #for i in range(0,1):

        csmotif = csmotif_data[i]

        ##QCP RMSD
        rmsd, transformed_coos = qcp.rmsdQCP(psmotif[0],csmotif, direction)

        if rmsd <= 2.0:

            ## Sequence filter, align native and smotif aa_seq as a measure of sequence similarity = structure similarity
            csse_seq, seq_identity, blosum62_score, bool_sequence_similarity \
            = Sfilter.S2SequenceSimilarity(current_ss, csmotif, direction, exp_data, threshold=40)


            ## Contacts filter,
            no_of_contacts, percent_of_satisfied_contacts \
            = Cfilter.S2ContactPredicition(transformed_coos, sse_ordered, exp_data)


            if bool_sequence_similarity and percent_of_satisfied_contacts > 50.0:
                print 'blosum62 score', blosum62_score, "seq_id", seq_identity, "rmsd=", rmsd, "Contacts", percent_of_satisfied_contacts
                dump_log.append([transformed_coos,['seq_filter', csse_seq, seq_identity, blosum62_score],
                ['contacts_filter', no_of_contacts, percent_of_satisfied_contacts], sse_ordered])

    if len(dump_log) > 0 :
        io.dumpPickle("tx_"+str(index_array[0])+"_"+str(index_array[1])+".pickle",dump_log)

    return True